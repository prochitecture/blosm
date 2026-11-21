"""Physical height reconciliation only for a new typed finite free closure."""
from .fixed_grass_heights import boundary_mask, reconcile as _reconcile


def reconcile(mesh, *, datum, ceiling, exaggeration=1., physical_mesh=None,
              max_cells=50000):
    return _reconcile(mesh, datum=datum, ceiling=ceiling, exaggeration=exaggeration,
        physical_mesh=physical_mesh, max_cells=max_cells,
        policy='fixed_boundary_finite_free_grass')


def physical_separator(context, plan):
    """Rebuild the same finite plan in a temporary original-physical view."""
    from .bundle_partial_caps import dependencies
    from .bundle_merged_plan import physical_context
    from .bundle_caps import finite_plan
    sections = set(plan.bundle.sections) | dependencies(context.network, plan.bundle)
    sites = tuple(site.key for site in context.network.intersections
                  if any(end.section_id in sections for end in site.ends))
    physical = physical_context(context, required_owners=sites)
    candidate = finite_plan(physical, plan.bundle)
    if (candidate is None or type(candidate) is not type(plan)
            or candidate.margin_factor != plan.margin_factor):
        raise ValueError('Finite free physical reconstruction changed its accepted finite plan')
    part = candidate.separator(physical)
    if 'finite_free_reconciliation' not in part.metadata:
        raise ValueError('Finite free physical reconstruction lost its typed closure')
    return part.mesh
