"""Keep a complete original mouth beyond its owner's last source overlap.

A transverse mouth can be wholly exposed while another, consumed source strip
crosses the same outgoing road farther away. Certify the whole finite remainder,
not just that one cross-section, before assigning its pavement to the road.
"""
import numpy as np

from .finite_containment import segment_on_polygon_side
from .finite_profile_search import _stations, advance
from .mesh import cross
from .model import End
from .regions import _classify
from .routes import Route


def _positive_overlap(cell, polygon, *, tolerance=1e-7):
    # All intersections partition each edge, so even a short excursion into
    # either interior is found. Boundary-only contacts keep the old mouth.
    if any(not segment_on_polygon_side(a,b,polygon,inside=False,tolerance=tolerance)
           for a,b in zip(cell,np.roll(cell,-1,axis=0))):
        return True
    selected=polygon[np.all(polygon>=cell.min(0),axis=1)
                     &np.all(polygon<=cell.max(0),axis=1)]
    if len(selected):
        inside,boundary=_classify(selected,cell,tolerance)
        if np.any(inside & (boundary<0)):return True
    return False


def advance_past_overlap(planner,approaches,outward,ports,sources,partitions,members,*,
                         maximum_cells=20000,maximum_pairs=2000000,max_candidates=4096):
    """Advance only actual original mouths with a positive remote overlap.

    Source cells use the same complete certified station partitions as the
    physical source union. Every original Path, incidence, width and full
    transverse profile stays unchanged. The ordinary finite-profile search
    applies its original neighbor, scope and finite-boundary guards to the new
    cut. The caller must rebuild the union and repeat every ownership proof.
    An incursion inside the two existing chord-error allowances does not prove
    a positive re-entry and cannot create a new ownership request.
    """
    if (type(maximum_cells) is not int or not 1<=maximum_cells<=20000
            or type(maximum_pairs) is not int or not 1<=maximum_pairs<=2000000
            or type(max_candidates) is not int or not 0<=max_candidates<=4096):
        raise ValueError('Source-mouth clearance needs bounded finite cells, pairs and candidates')
    if len(approaches)!=len(outward) or len(ports)!=len(approaches):
        raise ValueError('Source-mouth clearance lost its original port inventory')
    regions=[(s.region.owner,s.region.xy,s.region.xy.min(0),s.region.xy.max(0)) for s in sources]
    source_error=min(.001,.05*planner.context.network.settings.curve_tolerance)
    overlap_band=2.*source_error+1e-7
    result=list(outward);failed=[];witnesses={};cells=pairs=0
    for index,(a,current,port) in enumerate(zip(approaches,outward,ports)):
        if (not isinstance(a.path,Route) or hasattr(a,'original_approaches')
                or getattr(a.path,'absorbed_sites',())
                or getattr(a.path,'mouth_piece',None) is not None
                or getattr(a.path,'deferred_profile_reserve',None) is not None):
            continue
        values=_stations(planner.context,a,partitions)
        values=np.r_[current,values[values>current+1e-9]]
        if len(values)<2:continue
        cells+=len(values)-1
        if cells>maximum_cells:
            raise ValueError('Source-mouth clearance exceeds its finite cell budget')
        high=a.point(values,a.offset+.5*a.width)
        low=a.point(values,a.offset-.5*a.width)
        selected=[r for r in regions if r[0]!=port.region.owner]
        bottoms=np.array([r[2] for r in selected]);tops=np.array([r[3] for r in selected])
        last=None;owners=set();overlaps=0
        for i in range(len(values)-1):
            cell=np.array((high[i],high[i+1],low[i+1],low[i]))
            edge=np.roll(cell,-1,axis=0)-cell
            turns=cross(edge,np.roll(edge,-1,axis=0))
            if not (np.all(turns>1e-12) or np.all(turns<-1e-12)):
                raise ValueError('Source-mouth clearance needs unfolded finite source cells')
            near=(np.flatnonzero(np.all(cell.max(0)>=bottoms,axis=1)
                                &np.all(cell.min(0)<=tops,axis=1)) if selected else ())
            for other in near:
                pairs+=1
                if pairs>maximum_pairs:
                    raise ValueError('Source-mouth clearance exceeds its finite pair budget')
                if _positive_overlap(cell,selected[other][1],tolerance=overlap_band):
                    last=float(values[i+1]);owners.add(selected[other][0]);overlaps+=1
        if last is not None:
            from .intersections import IntersectionProvider
            if not IntersectionProvider.mouth_fits(a,last):
                route=planner.topology.trace(End(a.source.id,a.end))
                if (route.terminal_site in planner.layouts and route.terminal_site not in members
                        and not route.terminal_scope and not route.terminal_selection
                        and tuple((p.section.id,p.forward) for p in a.path.pieces)
                            ==tuple((v.section_id,v.forward) for v in route.traversals)):
                    from .compound_preparation import CompoundExpansionRequired
                    proof=dict(policy='ACTUAL_OWNED_STRIP_REACHES_ORIGINAL_FINITE_NEIGHBOR',
                        section=a.source.id,end=a.end,previous_station=float(current),
                        last_overlapping_cell_end=last,finite_route_length=a.length,
                        overlap_operands=[owner.text() for owner in sorted(owners)],
                        overlapping_pairs=overlaps,certified_remainder_cells=len(values)-1,
                        positive_overlap_source_error_band=overlap_band,
                        previous_members=[key.text() for key in sorted(members)],
                        original_neighbor=route.terminal_site.text(),endpoint_extrapolation=False,
                        terminal_scope=False,terminal_selection=False,proof_only_future_operands=False)
                    records=getattr(planner,'source_mouth_expansions',{})
                    records[tuple(sorted(members)),route.terminal_site,a.source.id,a.end]=proof
                    planner.source_mouth_expansions=records
                    failure=CompoundExpansionRequired((route.terminal_site,))
                    failure.source_mouth_clearance=proof
                    raise failure
            failed.append(index);result[index]=last
            witnesses[index]=dict(previous_station=float(current),last_overlapping_cell_end=last,
                overlap_operands=[owner.text() for owner in sorted(owners)],overlapping_pairs=overlaps,
                certified_remainder_cells=len(values)-1,finite_route_length=a.length)
    if not failed:return list(outward),[]
    if not max_candidates:
        raise ValueError('Source-mouth clearance exceeds its finite candidate budget')
    result,reports=advance(planner,approaches,result,failed,ports,sources,partitions,members,
                           max_candidates=max_candidates)
    for index,report in zip(failed,reports):
        report.update(witnesses[index],policy='COMPLETE_ORIGINAL_MOUTH_PAST_LAST_OWNER_SOURCE_OVERLAP',
            complete_remaining_strip='continuous_convex_cell_edges_and_operand_vertex_exclusion',
            boundary_only_contacts_preserved=True,positive_overlap_source_error_band=overlap_band,
            total_cells_tested=cells,total_pairs_tested=pairs)
    return result,reports
