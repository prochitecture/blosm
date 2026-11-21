"""Failed-only micrometre normalization of a proven source-cycle exterior."""
from .source_rim_normalization import normalize
from .source_strip_containment import certify


def recover(parent,group,context,original):
    layout=group.physical_layout
    if (layout is None or not hasattr(layout,'source_cycle_proofs')
            or not group.claims or len(group.loops)!=1 or group.loops[0].hole
            or getattr(layout,'source_rim_normalization',None)):
        raise ValueError('Contact normalization requires an unnormalized source-proven cycle exterior')
    sources={claim.section_id for claim in group.claims}
    before=certify(context,group,sources)
    normalized=normalize(group,maximum_projection=1e-6)
    prepared=certify(context,normalized,sources)
    # The existing normalization marks its immutable provenance. Repeated
    # normalization is ineligible, so another failed contact cannot recurse.
    result=parent.build(normalized,context)
    after=certify(context,result.partition.group,sources)
    result.metadata['source_width_visibility']=after
    result.metadata['source_contact_normalization']=dict(
        previous_contact_failure=str(original),original_width_proof=before,
        normalized_width_proof=prepared,final_width_proof=after,
        normalization=normalized.physical_layout.source_rim_normalization,
        maximum_added_polyline_projection=1e-6,original_profiles_exact=True,
        original_analytical_primitives_retained=True,contact_tolerance_unchanged=True)
    return result
