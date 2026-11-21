"""Recover analytical curb contacts hidden by microscopic internal mouth caps.

Clipping sampled curves can leave a tiny transverse cap of an already consumed
road between two actual curbs. It is neither a new external profile nor a road
feature. Recover their analytical intersection within the existing sampling
budget, then prove complete source-width containment before retaining it.
"""
from dataclasses import replace
import numpy as np
from .compounds import rim_label
from .regions import BoundaryLoop


def normalize(context,group):
    """Keep every external profile; normalize only bounded internal fragments.

    Unchanged, large, externally owned or unproved cases return the original
    group for its normal complete-exterior construction. No geometry failure
    is accepted by raising a tolerance or omitting an actual incoming port.
    """
    if group.physical_layout is not None:return group
    from .compound_ports import CompoundPortPlanner
    if not CompoundPortPlanner.incomplete_ports(group):return group
    from .sampling import AnalyticalBoundarySampler
    sources={f.region.owner:f for f in group.footprints}
    if any(getattr(f.layout,'absorbed_transitions',()) for f in group.footprints):return group
    claimed={claim.section_id for claim in group.claims}
    limit=min(.001,.2*context.network.settings.curve_tolerance)
    sampler=AnalyticalBoundarySampler(tolerance=context.network.settings.curve_tolerance,
        spacing=context.network.settings.spacing)
    loops=[];records=[];removed_points=[]
    for loop_index,loop in enumerate(group.loops):
        labels=[rim_label(span,sources) for span in loop.spans]
        tokens=[(span.owner,label.index) if label.kind=='mouth' else None
                for span,label in zip(loop.spans,labels)]
        remove=set()
        for start,token in enumerate(tokens):
            if token is None or token==tokens[start-1]:continue
            indices=[];at=start
            while tokens[at]==token:
                indices.append(at);at=(at+1)%len(tokens)
                if at==start:return group
            approach=sources[token[0]].layout.approaches[token[1]]
            if (getattr(approach.path,'absorbed_sites',())
                    or getattr(approach.path,'mouth_piece',None) is not None
                    or not all(piece.section.id in claimed for piece in approach.path.pieces)):
                continue
            first,last=labels[start].first,labels[indices[-1]].last
            if abs(first-1.)<=1e-8 and abs(last)<=1e-8:continue
            length=sum(float(np.linalg.norm(loop.xy[(i+1)%len(loop.xy)]-loop.xy[i])) for i in indices)
            previous=(start-1)%len(tokens)
            if not 0.<length<=limit or labels[previous].kind!='curb' or labels[at].kind!='curb':continue
            remove.update(indices)
            removed_points.append((loop_index,loop.xy[[*indices,at]]))
            records.append(dict(owner=token[0].text(),approach=token[1],section=approach.source.id,
                source_sections=[piece.section.id for piece in approach.path.pieces],
                source_parameter_interval=(float(first),float(last)),fragment_length=length,
                fragment_limit=limit,
                analytical_curve_tolerance=context.network.settings.curve_tolerance))
        keep=[i for i in range(len(loop.spans)) if i not in remove]
        loops.append(BoundaryLoop(loop.xy[keep],tuple(loop.spans[i] for i in keep)))
    if not records:return group
    filtered=replace(group,loops=tuple(loops))
    try:
        sampled=sampler.build(filtered)
        if CompoundPortPlanner.incomplete_ports(sampled):return group
        if len(sampled.loops)!=len(group.loops):return group
        for record,(loop_index,points) in zip(records,removed_points):
            boundary=sampled.loops[loop_index]
            if boundary.hole!=group.loops[loop_index].hole:return group
            first=boundary.xy;delta=np.roll(first,-1,axis=0)-first
            squared=np.sum(delta*delta,axis=1)
            if np.any(squared<=1e-16):return group
            offset=points[:,None]-first
            t=np.clip(np.sum(offset*delta,axis=2)/squared,0.,1.)
            movement=float(np.linalg.norm(offset-t[...,None]*delta,axis=2).min(axis=1).max())
            if movement>context.network.settings.curve_tolerance:return group
            record['removed_cap_to_refined_boundary_error']=movement
        from .conditions import BoundaryConditions
        layouts={owner:f.layout for owner,f in sources.items()}
        conditions=[BoundaryConditions().build(sampled,loop,layouts) for loop in sampled.loops]
        from .source_strip_containment import certify
        proof=certify(context,sampled,claimed)
    except ValueError:return group
    from .intersections import JunctionLayout,Mouth,half_intervals
    from .model import Key
    ports=[port for boundary in conditions for port in boundary.ports]
    approaches=tuple(layouts[port.owner].approaches[port.approach_index] for port in ports)
    distances=tuple(layouts[port.owner].outward_stations[port.approach_index] for port in ports)
    radial=half_intervals(context.network.settings);mouths={}
    for layout in layouts.values():
        for key,mouth in layout.mouths.items():
            previous=mouths.get(key)
            if previous is None or (mouth.station>previous.station if key[1]==0 else mouth.station<previous.station):
                mouths[key]=mouth
    for sid in claimed:
        mouths[sid,0]=Mouth(context[Key('path',sid)].length,2*radial)
        mouths[sid,1]=Mouth(0.,2*radial)
    physical=JunctionLayout(approaches,tuple(c for layout in layouts.values() for c in layout.corners),(),
        distances,mouths,len(sampled.loops[0].xy),radial,np.zeros(2),[0.]*len(approaches),mesh_strategy='boundary_chart')
    physical.source_region_union=True
    physical.source_curb_owner_canonicalization=getattr(sampler,'source_curb_canonicalizations',())
    physical.internal_profile_fragments=dict(fragments=records,complete_source_width=proof,
        policy='ANALYTICAL_CURB_CONTACT_REFINEMENT',external_profiles_preserved=True)
    return replace(sampled,physical_layout=physical)
