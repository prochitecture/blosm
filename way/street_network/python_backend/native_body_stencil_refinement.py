"""Failed-only conforming discretization of one proven smooth native body.

The caller retains every original Path/profile/owned interval and every old
mesh vertex. Only new computational points may move; their heights are always
reauthenticated against the actual finite original body map.
"""
import numpy as np
from .model import Key,fingerprint
from .source_body_heights import prepare_steps
from .updates import drain_steps
from .source_body_references import BodyGuideLimits
from .slope_projection import corner_matrices,ProjectionIterationBudgetError
from .refinement import QuadRefiner
from .interior_quad_refinement import InteriorQuadRefinement
from .quality_correction import QuadShapeCorrection
from .entrance_height_recovery import recover_steps as recover_height_steps
from .finite_containment import contains_segment,_capsule_intervals
from .four_side_charts import quality
from .publication_slopes import PublicationSlopeGuard
from .junction_bank_transitions import certify_body_domains
from .source_bank_reconciliation import owned_ranges
from .source_curb_coordinates import CurbCoordinates
from .source_identity import matches_current_section

POLICY = '1-conforming-native-body-stencils-with-fixed-original-guides'

class NativeBodyStencilUnavailable(ValueError):pass

def _bad(mesh,body,fixed,limit):
    """Find faces whose fixed heights make a corner gradient infeasible."""
    known = fixed.copy()
    known[body.core_nodes] = True
    heights = mesh.xyz[:, 2].copy()
    heights[body.core_nodes] = body.core_targets
    matrices = corner_matrices(mesh)
    weights = np.array((1, 2, 4, 8))
    masks = known[mesh.quads] @ weights
    bounds = np.zeros(len(mesh.quads))
    for mask in np.unique(masks):
        if mask == 0:
            continue
        indices = np.flatnonzero(masks == mask)
        fixed_columns = (mask & weights) != 0
        for start in range(0, len(indices), 512):
            rows = indices[start:start+512]
            corners = matrices[rows]
            fixed_heights = heights[mesh.quads[rows][:, fixed_columns]]
            residual = np.einsum('fcij,fj->fci', corners[..., fixed_columns], fixed_heights)
            if not fixed_columns.all():
                free = corners[..., ~fixed_columns]
                # Project onto the orthogonal complement of the gradients
                # attainable by changing the free vertices of this quad.
                without_fixed = np.linalg.pinv(free, rcond=1e-12) @ -residual[..., None]
                residual += (free @ without_fixed)[..., 0]
            bounds[rows] = np.linalg.norm(residual, axis=2).max(axis=1)
    return np.flatnonzero(bounds > limit+1e-9).astype(np.int32), float(np.max(bounds, initial=0.))

def recover(builder,context,group,mesh,refs,conditions,field,body,limit,max_cells,*,previous_failure):
    return drain_steps(recover_steps(builder,context,group,mesh,refs,conditions,field,body,limit,max_cells,
                                     previous_failure=previous_failure))


def recover_steps(builder,context,group,mesh,refs,conditions,field,body,limit,max_cells,*,previous_failure):
    """Return a locally qualified remesh to the unchanged publication path."""
    if type(previous_failure) is not ProjectionIterationBudgetError:
        raise NativeBodyStencilUnavailable('Native body refinement requires the actual typed default iteration refusal')
    if type(max_cells) is not int or not len(mesh.quads)<=max_cells<=50000:
        raise NativeBodyStencilUnavailable('Native body refinement lacks the original finite cell allocation')
    maximum_passes=min(3,builder.max_refinements)
    if maximum_passes<1:raise NativeBodyStencilUnavailable('Native body refinement has no configured refinement allocation')
    if not body.ledger['complete_supported_scope'] or body.ledger['unsupported']:
        raise NativeBodyStencilUnavailable('Native body refinement cannot omit an unsupported original body domain')
    for reference in body.references:
        sid=reference.section_id
        if (reference.path is not context[Key('path',sid)] or reference.profile is not context[Key('profile',sid)]
                or not matches_current_section(context, reference.path.section)):
            raise NativeBodyStencilUnavailable('Native body refinement changed its exact original dependency identity')
    network_token=fingerprint(context.network)
    original_mesh=mesh;original_core=body.core_nodes.copy();original_targets=body.core_targets.copy()
    original_fixed=np.zeros(len(mesh.xyz),bool);original_fixed[list(refs)]=True
    original_rim=mesh.xyz[original_fixed].copy();current=body.seed;refs=dict(refs)
    mapping=np.arange(len(mesh.xyz),dtype=np.int32);fixed=original_fixed.copy()
    selected,lower=_bad(current,body,fixed,limit)
    if not 1<=len(selected)<=512:
        raise NativeBodyStencilUnavailable('Failed original body solve has no bounded locally infeasible stencil witness')
    hard=np.zeros(len(current.xyz),bool);hard[body.core_nodes]=True
    nodes=np.unique(current.quads[selected])
    if not hard[nodes].all() or fixed[nodes].any():
        raise NativeBodyStencilUnavailable('Native stencil refinement cannot replace a rim or mixed-source anchor conflict')
    references=[]
    for reference in body.references:
        station,q,valid,proof=reference.inverse(current.xyz[nodes,:2],max_cells=max_cells)
        if valid.all() and np.max(abs(reference.profile.datum+reference.profile.value(q)-current.xyz[nodes,2]))<=1e-8:
            references.append(reference)
    if len(references)!=1:
        raise NativeBodyStencilUnavailable('Hard stencils lack one unique actual finite original body map')
    reference=references[0];path=reference.path;profile=reference.profile
    ranges=owned_ranges(context,group).get(reference.section_id,())
    if not any(a<=reference.first<reference.last<=b for a,b in ranges):
        raise NativeBodyStencilUnavailable('Native body refinement leaves actual finite original ownership')
    if path.exact or profile.half:
        raise NativeBodyStencilUnavailable('This bounded stencil template requires a smooth original full crown')
    gradient=4.*profile.rise/path.section.width
    regularity=path.curvature_bound*(abs(path.section.offset)+.5*path.section.width)
    if gradient>limit or regularity>=1.-1e-10:
        raise NativeBodyStencilUnavailable('Original continuous body map lacks the unchanged gradient/regularity certificate')
    source=next((area for area in group.area_footprints if area.region.owner==Key('cycle_strip',reference.section_id)),None)
    if source is None:raise NativeBodyStencilUnavailable('Native body map has no actual source area operand')
    source_xy=source.region.xy.copy();coordinate_sources=[]
    for curve in source.layout.boundaries:
        coordinates=getattr(curve,'source_coordinates',None)
        if (type(coordinates) is not CurbCoordinates or coordinates.path is not path or coordinates.approach is not None):
            raise NativeBodyStencilUnavailable('Native body stencil lacks actual single-Path typed source bank coordinates')
        coordinates.validate(curve.points,context);coordinate_sources.append((coordinates,curve))
    budgets=[curve.source_coordinates.maximum_chord_error for curve in source.layout.boundaries
        if getattr(curve,'source_coordinates',None) is not None]
    if not budgets or any(value is None for value in budgets):
        raise NativeBodyStencilUnavailable('Native body source lost its continuous chord-error declaration')
    budget=max(budgets)
    for fi in selected:
        polygon=current.xyz[current.quads[fi],:2]
        if any(not contains_segment(P,Q,source.region.xy,tolerance=1e-10)
                or _capsule_intervals(P,Q-P,source.region.xy,budget)
                for P,Q in zip(polygon,np.roll(polygon,-1,axis=0))):
            raise NativeBodyStencilUnavailable('A complete native stencil leaves its finite original strip/chord interior')
    details=dict(policy=POLICY,previous_original_body_failure=str(previous_failure),
        source_section=reference.section_id,continuous_native_gradient_bound=gradient,
        normal_offset_curvature_product=regularity,continuous_finite_stencil_containment=True,
        original_native_chord_budget=budget,original_quads=len(mesh.quads),original_body_targets=len(original_core),
        initial_infeasible_faces=selected.tolist(),initial_local_lower_bound=lower,
        maximum_passes=maximum_passes,maximum_candidates=512,maximum_cells=max_cells,passes=[])

    def rebuild(refined,remap,family,extra):
        nonlocal current,refs,fixed,body,mapping,selected,lower
        old=np.zeros(len(refined.xyz),bool);old[remap]=True;added=np.flatnonzero(~old)
        refined.xyz[remap]=current.xyz
        refined.xyz[added,2]=reference.values(refined.xyz[added,:2])
        mapping=remap[mapping];refs={int(remap[node]):ref for node,ref in refs.items()}
        fixed=np.zeros(len(refined.xyz),bool);fixed[list(refs)]=True
        if not np.array_equal(refined.xyz[mapping[np.flatnonzero(original_fixed)]],original_rim):
            raise NativeBodyStencilUnavailable('Native refinement moved an original physical rim sample')
        body=yield from prepare_steps(context,group,refined,refs,limits=BodyGuideLimits(maximum_mesh_cells=max_cells))
        certify_body_domains(body,conditions)
        targets={int(node):float(value) for node,value in zip(body.core_nodes,body.core_targets)}
        if any(int(node) not in targets or abs(targets[int(node)]-value)>1e-12
                for node,value in zip(mapping[original_core],original_targets)):
            raise NativeBodyStencilUnavailable('Native refinement lost or changed an original required body target')
        selected,lower=_bad(refined,body,fixed,limit);current=refined
        details['passes'].append(dict(family=family,quads=len(current.quads),added_vertices=len(added),
            remaining_infeasible_faces=len(selected),maximum_local_lower_bound=lower,**extra))

    for _ in range(maximum_passes):
        if not len(selected):break
        if len(selected)>512:raise NativeBodyStencilUnavailable('Native stencil count exceeds unchanged 512-candidate allocation')
        refined,remap=yield from QuadRefiner(context.patterns,max_cells=max_cells).refine_steps(current,selected,lambda xy:np.zeros(len(xy)))
        yield from rebuild(refined,remap,'existing_conforming_1_to_3',{})
    if len(selected):raise NativeBodyStencilUnavailable('Existing three refinement passes leave infeasible native body stencils')
    lengths=np.linalg.norm(np.roll(current.xyz[current.quads],-1,axis=1)-current.xyz[current.quads],axis=2)
    quality_bad=np.flatnonzero(lengths.max(axis=1)/lengths.min(axis=1)>8.)
    if len(quality_bad):
        if len(details['passes'])>=maximum_passes or len(quality_bad)>512:
            raise NativeBodyStencilUnavailable('Native quality refinement exceeds the same three-pass/512-candidate allocation')
        previous=current,refs,fixed,body,mapping
        refined,remap,split=InteriorQuadRefinement(context.patterns,max_cells=max_cells)._refine(current,quality_bad)
        yield from rebuild(refined,remap,'existing_anisotropic_long_edge',dict(split_interior_edges=split))
        if len(selected):
            from .native_stencil_alignment import align
            # The alternate topology can make a smooth profile infeasible.
            # Return to the last slope-feasible mesh before adjusting its new vertices.
            current,refs,fixed,body,mapping=previous
            details['rejected_long_edge_transition']=details['passes'].pop()
            current,alignment=align(current,mapping,fixed,body,reference,limit)
            body=yield from prepare_steps(context,group,current,refs,limits=BodyGuideLimits(maximum_mesh_cells=max_cells))
            certify_body_domains(body,conditions)
            selected,lower=_bad(current,body,fixed,limit)
            if len(selected):
                raise NativeBodyStencilUnavailable('Native vertex alignment leaves an infeasible original body stencil')
            details['passes'].append(dict(family='align_added_vertices_on_original_profile',
                quads=len(current.quads),maximum_local_lower_bound=lower,**alignment))
    # Preserve the already qualified local sequence: finish the ordinary
    # body solve on this topology before adjusting only introduced XY points.
    current,pre_shape_reconciliation=yield from recover_height_steps(builder,body.seed,field,fixed,conditions,limit,max_cells,
        at_corners=False,additional_slabs=body.height_slabs,body_guides=body)
    planar=current.with_positions(np.c_[current.xyz[:,:2],np.zeros(len(current.xyz))])
    corrected,correction=yield from QuadShapeCorrection(max_cells=max_cells).improve_steps(planar,fixed_vertices=mapping)
    if not correction['qualified']:raise NativeBodyStencilUnavailable('Fixed-original shape correction did not meet full quad quality')
    added=np.ones(len(current.xyz),bool);added[mapping]=False
    corrected.xyz[mapping,2]=current.xyz[mapping,2]
    corrected.xyz[added,2]=reference.values(corrected.xyz[added,:2])
    if (not np.array_equal(corrected.xyz[mapping,:2],original_mesh.xyz[:,:2])
            or not np.array_equal(corrected.xyz[mapping[np.flatnonzero(original_fixed)]],original_rim)):
        raise NativeBodyStencilUnavailable('Shape correction moved an original guide or physical profile')
    body=yield from prepare_steps(context,group,corrected,refs,limits=BodyGuideLimits(maximum_mesh_cells=max_cells))
    certify_body_domains(body,conditions)
    targets={int(node):float(value) for node,value in zip(body.core_nodes,body.core_targets)}
    if any(int(node) not in targets or abs(targets[int(node)]-value)>1e-12
            for node,value in zip(mapping[original_core],original_targets)):
        raise NativeBodyStencilUnavailable('Shape correction changed an original required body target')
    selected,lower=_bad(corrected,body,fixed,limit)
    if len(selected):raise NativeBodyStencilUnavailable('Shape correction reintroduces a hard native body conflict')
    final,reconciliation=yield from recover_height_steps(builder,body.seed,field,fixed,conditions,limit,max_cells,
        at_corners=False,additional_slabs=body.height_slabs,body_guides=body)
    physical=quality(final);posted=quality(PublicationSlopeGuard.physical_buffer(final,field.datum,
        context.network.settings.exaggeration or 1.))
    if (max(physical['maximum_edge_ratio'],posted['maximum_edge_ratio'])>8.
            or min(physical['minimum_corner_sine'],posted['minimum_corner_sine'])<.03):
        raise NativeBodyStencilUnavailable('Refined original body failed physical/stored quad quality')
    if (not np.array_equal(final.xyz[mapping,:2],original_mesh.xyz[:,:2])
            or not np.array_equal(final.xyz[mapping[np.flatnonzero(original_fixed)]],original_rim)):
        raise NativeBodyStencilUnavailable('Refined body solve moved original geometry or physical profiles')
    error=float(np.max(abs(final.xyz[mapping[original_core],2]-original_targets)))
    if error>1e-8 or body.height_slabs.violation(final.xyz[:,2])>1e-8:
        raise NativeBodyStencilUnavailable('Refined body solve lost an original native height constraint')
    body.assert_unchanged(final)
    if fingerprint(context.network)!=network_token or not np.array_equal(source.region.xy,source_xy):
        raise NativeBodyStencilUnavailable('Original source buffers changed during a finite native body refinement')
    for coordinates,curve in coordinate_sources:coordinates.validate(curve.points,context)
    details.update(pre_shape_height_reconciliation=pre_shape_reconciliation,
        shape_correction=correction,physical_quality=physical,stored_quality=posted,
        original_XY_and_rim_XYZ_exact=True,original_body_targets_retained=True,original_body_target_error=error,
        final_quads=len(final.quads),final_local_lower_bound=lower,only_introduced_native_points_resampled=True,
        no_profile_radius_geometry_or_limit_changes=True,normal_publication_still_required=True,
        height_reconciliation=reconciliation)
    return final,refs,fixed,body,details
