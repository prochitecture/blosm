"""Finite source admission for a whole member beside a truly free member."""
import math
import numpy as np
from .bundles import CenterlineRange, certify_disjoint
from .cycle_dimensional_retry import certify_layout
from .finite_containment import separated_segment
from .mesh import cross
from .model import End, Key, fingerprint
from .regions import _classify
from .bundle_axis import AxialChartUnavailable


class WhollyOwnedBundleUnavailable(AxialChartUnavailable):
    """The explicit owned/free preparation contract is not available."""


def candidate(context,bundle,consumed,members):
    """An original graph policy only; final geometry must still certify."""
    if len(bundle.members)!=2 or not all(m.complete and len(m.path)==1 for m in bundle.members):
        return None
    visits=tuple(m.path[0] for m in bundle.members)
    included=[v.section_id in consumed for v in visits]
    if sum(included)!=1:return None
    role=included.index(True);other=1-role
    owned,free=(context.network.section_index[visits[i].section_id] for i in (role,other))
    if owned.scope_ends or owned.selection_ends or free.scope_ends or free.selection_ends:return None
    incident={end:site for site in (*context.network.intersections,*context.network.connections) for end in site.ends}
    roots=(End(owned.id,0),End(owned.id,1))
    if any(end not in incident or len(incident[end].ends)<3 or incident[end].key not in members for end in roots):
        return None
    if any(end.section_id==free.id for end in incident):return None
    paths=tuple(context[Key('path',v.section_id)] for v in visits)
    separation=certify_disjoint(*((CenterlineRange.build(path,0.,path.length),) for path in paths))
    return dict(bundle=bundle.id,owned_role=role,owned_section=owned.id,free_section=free.id,
        original_owner_sites=[incident[end].key.text() for end in roots],
        complete_original_members=True,true_free_source_ends=True,source_aliases_created=0,
        original_centerline_separation=separation)


def certify(context,group,*,maximum_cells=8192,store=True):
    """Keep the entire free ribbon outside the OUTER region, including holes."""
    consumed={claim.section_id for claim in group.claims}
    candidates=[record for bundle in context.network.bundles
                if (record:=candidate(context,bundle,consumed,set(group.members))) is not None]
    if not candidates:return group
    if group.physical_layout is None:
        raise WhollyOwnedBundleUnavailable('Owned/free source certification needs a prepared physical layout')
    if type(maximum_cells) is not int or maximum_cells<1:
        raise ValueError('Invalid owned/free source cell allocation')
    exterior=[loop.xy for loop in group.loops if not loop.hole]
    if len(exterior)!=1:raise ValueError('Owned/free source preparation needs one actual physical exterior')
    outer=exterior[0]
    complete_width=certify_layout(context,group,consumed)
    rows=0;reports=[]
    for record in candidates:
        path=context[Key('path',record['free_section'])];section=path.section
        epsilon=min(.00025,.05*context.network.settings.curve_tolerance)
        if path.exact:
            stations=np.asarray(path.knots);chord_error=0.
        else:
            k=float(path.curvature_bound);reach=abs(section.offset)+.5*section.width
            if reach*k>=1.:raise ValueError('Free source ribbon needs a positive offset curvature certificate')
            curvature=k/(1.-reach*k);speed=1.+reach*k
            spacing=min(.5,context.network.settings.spacing)
            if curvature>1e-12:
                spacing=min(spacing,math.sqrt(8.*epsilon/curvature)/speed,math.pi/(16.*curvature*speed))
            if path.length/spacing>maximum_cells-rows:
                raise ValueError('Owned/free source ribbon exceeds its finite width-cell budget')
            stations=path.samples(0.,path.length,spacing,tolerance=epsilon,max_samples=maximum_cells-rows+1)
            chord_error=epsilon
        rows+=len(stations)-1
        if rows>maximum_cells:raise ValueError('Owned/free source ribbon exceeds its finite width-cell budget')
        high=path.point(stations,section.offset+.5*section.width)
        low=path.point(stations,section.offset-.5*section.width)
        cells=np.stack((high[:-1],high[1:],low[1:],low[:-1]),axis=1)
        edge=np.roll(cells,-1,axis=1)-cells
        if np.min(cross(edge,np.roll(edge,-1,axis=1)))<=1e-12:
            raise ValueError('Owned/free source ribbon needs unfolded complete source cells')
        outer_allowance=context.network.settings.curve_tolerance
        clearance=outer_allowance+chord_error+1e-7
        for cell in cells:
            if any(not separated_segment(a,b,outer,clearance=clearance)
                   for a,b in zip(cell,np.roll(cell,-1,axis=0))):
                raise ValueError('Wholly owned source region overlaps or approaches the complete free counterpart ribbon')
            # A tiny owned component enclosed in a source cell must not evade
            # an edge-only test. Holes are deliberately not subtracted here.
            selected=np.flatnonzero(np.all(outer>=cell.min(axis=0),axis=1)&np.all(outer<=cell.max(axis=0),axis=1))
            if len(selected):
                inside,near=_classify(outer[selected],cell,1e-10)
                if np.any(inside|(near>=0)):
                    raise ValueError('The complete free counterpart ribbon encloses an owned outer boundary')
        reports.append(dict(record,free_width_cells=len(cells),source_chord_error=chord_error,
            free_source_width=section.width,free_source_offset=section.offset,free_source_length=path.length,
            free_source_fingerprint=fingerprint(section),
            actual_outer_curve_allowance=outer_allowance,numeric_allowance=1e-7,
            certified_positive_clearance=clearance,outer_region_including_holes_excluded=True,
            proof='continuous_convex_full_width_cells_outside_outer_polygon_and_boundary_capsules'))
    certificate=dict(candidates=reports,total_free_width_cells=rows,
        complete_consumed_and_prefix_width=complete_width,opposite_source_remains_independent=True,
        acquired_opposite_owners=0,physical_source_connections_created=0)
    if not store:return certificate
    group.physical_layout.wholly_owned_free_counterparts=certificate
    group.physical_layout.wholly_owned_free_source_dependencies=tuple(sorted({row['free_section'] for row in reports}))
    return group
