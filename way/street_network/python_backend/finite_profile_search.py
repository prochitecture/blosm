"""Bounded complete-mouth search on existing finite original source stations."""
import numpy as np
from .finite_containment import segment_on_polygon_side,separated_segment
from .intersections import IntersectionProvider
from .model import End,Key
from .routes import Route


class FiniteProfileUnavailable(ValueError):
    """The existing complete-mouth family exhausts its finite source interval."""


def _stations(context,approach,partitions):
    if not getattr(partitions,'complete_source_ranges',True):
        raise ValueError('Finite profile search needs complete certified original station arrays')
    path=approach.path
    if (not isinstance(path,Route) or getattr(path,'absorbed_sites',())
            or getattr(path,'mouth_piece',None) is not None
            or getattr(path,'deferred_profile_reserve',None) is not None):
        raise ValueError('Finite profile search requires an unmodified original source route')
    result=[]
    for i,piece in enumerate(path.pieces):
        if (piece.path is not context[Key('path',piece.section.id)]
                or abs(piece.section.width-approach.width)>1e-8
                or abs(piece.section.offset*(1. if piece.forward else -1.)-approach.offset)>1e-8):
            raise ValueError('Finite profile search changed original Path identity or dimensions')
        values=np.asarray(partitions[piece.section.id],float)
        if (values.ndim!=1 or len(values)<2 or not np.isfinite(values).all()
                or abs(values[0])>1e-9 or abs(values[-1]-piece.path.length)>1e-8
                or np.any(np.diff(values)<=0.)):
            raise ValueError('Finite profile search needs complete certified original station arrays')
        if not piece.forward:values=piece.path.length-values[::-1]
        result.extend(float(path.stations[i]+s) for s in values)
    return np.unique(result)


def advance(planner,approaches,outward,missing,ports,sources,partitions,members,*,max_candidates=4096,
            full_source_exclusion=False):
    """Advance only failed mouths after the ordinary lead search exhausts.

    Every candidate is an already certified original-Path station. The entire
    actual transverse profile must lie outside every other current source
    operand, using continuous segment/polygon containment. The caller must
    rebuild its full union and repeat profile, root, width and ownership guards.
    """
    if type(max_candidates) is not int or not 1<=max_candidates<=4096:
        raise ValueError('Finite profile search needs at most 4096 candidates')
    if type(full_source_exclusion) is not bool:
        raise ValueError('Finite profile exclusion needs an explicit boolean policy')
    if len(approaches)!=len(outward) or len(ports)!=len(approaches):
        raise ValueError('Finite profile search lost its original physical port inventory')
    result=list(outward);reports=[];tested=0;members=set(members)
    regions=[(s.region.owner,s.region.xy,s.region.xy.min(0),s.region.xy.max(0),None) for s in sources]
    source_error=min(.001,.05*planner.context.network.settings.curve_tolerance)
    clearance=2.*source_error+1e-7
    cells=0
    if full_source_exclusion:
        # These are exclusion witnesses only. They never join the physical
        # pavement union or claim an unowned source interval. Two source chord
        # allowances cover a later compulsory re-partition of the same curb.
        for port,a in zip(ports,approaches):
            for original in getattr(a,'original_approaches',(a,)):
                at=_stations(planner.context,original,partitions)
                right=original.point(at,original.offset+.5*original.width)
                left=original.point(at,original.offset-.5*original.width)
                for i in range(len(at)-1):
                    polygon=np.array((right[i],right[i+1],left[i+1],left[i]))
                    regions.append((port.region.owner,polygon,polygon.min(0),polygon.max(0),clearance));cells+=1
    for index in missing:
        a=approaches[index];path=a.path
        if not isinstance(path,Route) or hasattr(a,'original_approaches'):
            raise ValueError('Finite profile search requires an unmodified original source route')
        candidate=planner.topology.trace(End(a.source.id,a.end))
        visits=tuple((v.section_id,v.forward) for v in candidate.traversals)
        actual=tuple((p.section.id,p.forward) for p in path.pieces)
        if (actual!=visits[:len(actual)] or not actual
                or any(p.path is not planner.context[Key('path',p.section.id)] for p in path.pieces)):
            raise ValueError('Finite profile search changed original source incidence or Path identity')
        complete=len(actual)==len(visits)
        if (path.terminal_free and (not complete or candidate.terminal_scope or candidate.terminal_selection
                                   or candidate.terminal_site is not None)):
            raise ValueError('Finite profile search cannot relabel a finite scope or physical owner as free')
        neighbor=candidate.terminal_site if complete else None
        if neighbor is not None and neighbor not in members:
            layout=planner.layouts.get(neighbor)
            if layout is None or not hasattr(layout,'mouths'):
                raise ValueError('Finite profile search needs available original neighbor ownership')
        values=_stations(planner.context,a,partitions)
        values=values[values>outward[index]+1e-9];selected=None
        selected_regions=[r for r in regions if r[0]!=ports[index].region.owner]
        bottoms=np.array([r[2]-(r[4] or 1e-7) for r in selected_regions])
        tops=np.array([r[3]+(r[4] or 1e-7) for r in selected_regions])
        for proposed in values:
            station=IntersectionProvider.mouth_station(a,proposed,planner.context.network.settings)
            if not IntersectionProvider.mouth_fits(a,station):continue
            if tested>=max_candidates:
                raise ValueError('Finite profile search exceeds its bounded candidate budget')
            tested+=1
            profile=a.point(station,np.array((a.offset-.5*a.width,a.offset+.5*a.width)))
            low,high=profile.min(0),profile.max(0)
            near=(np.flatnonzero(np.all(high>=bottoms,axis=1)&np.all(low<=tops,axis=1))
                  if selected_regions else ())
            if any(not (separated_segment(*profile,selected_regions[i][1],clearance=selected_regions[i][4])
                        if selected_regions[i][4] is not None else
                        segment_on_polygon_side(*profile,selected_regions[i][1],inside=False,tolerance=1e-7))
                   for i in near):continue
            # The same directed original mouths used by normal source-union
            # closure decide whether an available neighboring owner must join.
            trims=path.trims(station);neighbors=set()
            for key,layout in planner.layouts.items():
                if key in members or not hasattr(layout,'mouths'):continue
                for (sid,end),at in trims.items():
                    opposite=layout.mouths.get((sid,1-end))
                    if opposite is not None and (at>opposite.station+1e-7 if end==0 else at<opposite.station-1e-7):
                        neighbors.add(key)
            if neighbors:
                from .compound_preparation import CompoundExpansionRequired
                raise CompoundExpansionRequired(neighbors)
            selected=station;break
        if selected is None:
            raise FiniteProfileUnavailable('Finite profile search has no complete mouth before its original finite boundary')
        result[index]=selected
        reports.append(dict(policy='FINITE_ORIGINAL_CERTIFIED_STATION_PROFILE_SEARCH',
            section=a.source.id,end=a.end,previous_station=float(outward[index]),station=float(selected),
            finite_route_length=a.length,candidates_tested=tested,
            complete_transverse_profile='continuous_outside_all_other_source_operands',
            numerical_boundary_band=1e-7,original_sources_profiles_and_corners_unchanged=True,
            endpoint_extrapolation=False,terminal_scope=candidate.terminal_scope,
            terminal_selection=candidate.terminal_selection,full_finite_source_exclusion=full_source_exclusion,
            exclusion_only_source_cells=cells,source_chord_error=source_error,
            future_source_chord_repartition_clearance=clearance if full_source_exclusion else None))
    return result,reports
