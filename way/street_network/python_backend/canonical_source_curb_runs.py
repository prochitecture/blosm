"""Canonical ownership of a short run on an identical original source curb."""
from copy import copy
from dataclasses import replace
import numpy as np
from .compounds import rim_label,RimSpan,RimSegment
from .regions import BoundaryLoop
from .coincident_source_curbs import _intervals,source_station,curve_parameter,same_bank


def canonicalize(group,maximum_length,*,maximum_runs=32):
    """Remove only a redundant A/B/A ownership island, not physical geometry.

    Both labels must map continuously to the identical original Path and offset
    throughout the finite island. The surrounding A curve must contain every
    mapped station. The normal analytical sampler can then sample one continuous
    curb without introducing tiny cells at incidental owner changes.
    """
    from .source_curb_coordinates import canonical_native_banks
    group,native_banks=canonical_native_banks(group)
    sources={f.region.owner:f for f in group.footprints};loops=[];records=[]
    for loop in group.loops:
        spans=list(loop.spans)
        for _ in range(maximum_runs-len(records)):
            labels=[rim_label(s,sources) for s in spans]
            tokens=[(s.owner,l.kind,l.index) for s,l in zip(spans,labels)]
            changed=False
            for start,token in enumerate(tokens):
                if token[1]!='curb' or token==tokens[start-1]:continue
                indices=[];end=start
                while tokens[end]==token:
                    indices.append(end);end=(end+1)%len(tokens)
                    if end==start:break
                if end==start:continue
                previous=(start-1)%len(tokens);target=tokens[previous]
                if target!=tokens[end] or target[1]!='curb' or target==token:continue
                length=sum(float(np.linalg.norm(loop.xy[(i+1)%len(spans)]-loop.xy[i])) for i in indices)
                if not 0.<length<=maximum_length:continue
                old_curve=sources[token[0]].layout.boundaries[token[2]]
                new_curve=sources[target[0]].layout.boundaries[target[2]]
                target_rows=_intervals(new_curve,0.,1.,.5)
                mapped=[];source=None;maximum_error=0.
                for index in indices:
                    label=labels[index];rows=_intervals(old_curve,label.first,label.last,.5*(label.first+label.last))
                    match=None
                    for row in rows:
                        first,last=source_station(row,label.first),source_station(row,label.last)
                        if min(first,last)<row['low']-1e-9 or max(first,last)>row['high']+1e-9:continue
                        if max(abs(curve_parameter(row,s)-t) for s,t in ((first,label.first),(last,label.last)))>1e-9:continue
                        for candidate in target_rows:
                            if not same_bank(candidate,row):continue
                            if min(first,last)<candidate['low']-1e-9 or max(first,last)>candidate['high']+1e-9:continue
                            a,b=curve_parameter(candidate,first),curve_parameter(candidate,last)
                            if b<=a:continue
                            error=max(float(np.linalg.norm(new_curve(t)-old_curve(u)))
                                for t,u in ((a,label.first),(b,label.last)))
                            if error>1e-7:continue
                            match=RimSpan(target[0],RimSegment('curb',target[2],a,b))
                            source=row['path'].section.id;maximum_error=max(maximum_error,error)
                            break
                        if match is not None:break
                    if match is None:break
                    mapped.append(match)
                if len(mapped)!=len(indices):continue
                for index,span in zip(indices,mapped):spans[index]=span
                records.append(dict(previous_owner=token[0].text(),canonical_owner=target[0].text(),
                    source_section=source,length=length,endpoint_error=maximum_error,
                    proof='IDENTICAL_ORIGINAL_PATH_OFFSET_AND_FINITE_STATION_DOMAIN'))
                changed=True;break
            if not changed:break
        loops.append(BoundaryLoop(loop.xy,tuple(spans)))
    if not records:return group,()
    physical=group.physical_layout
    if physical is not None:
        physical=copy(physical)
        physical.source_curb_owner_canonicalization=tuple(records)
    return replace(group,loops=tuple(loops),physical_layout=physical),tuple(records)
