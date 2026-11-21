# Complete source-strip ownership

An owned station interval and a contained centreline do not prove ownership of
the road's full width. A compound can retain both complete incoming profiles
yet clip a carriageway curb between them. Partial or whole Bundle claims now
require an independent complete source-strip proof before they can authorize
removing original road bodies.

`source_strip_containment.certify` works in physical metres on the unchanged
original paths. Smooth offset curbs use their complete curvature/speed bounds
to choose chords with at most 0.25 mm source error. Exact paths retain their
original miter knots. Consecutive transverse rows form convex source cells;
folded or excessive cell inventories are rejected.

`certify_intervals` preserves each exact owned station range on that same
original Path. Overlapping or exactly touching ranges merge; disjoint endpoint
prefixes of a returning Section do not authorize its middle. Only a real whole
consumed claim adds `[0,L]`. The wrapper batches at most one range per Section
through the unchanged proof and charges every batch to one 8,192-cell budget.
Its reports retain every certified interval explicitly. A single merged range
per Section uses the original call and report unchanged.
The shared `SOURCE_RANGE_POLICY` invalidates PLAN and the projected layout,
surface and seam signatures. Raw Path/Profile/physical-footprint resources keep
their original signatures; stale derived geometry is reevaluated normally.

Every cell edge is checked continuously against the single physical outer
boundary and every known hole. Boundary intersections split a segment into
intervals of constant inside/outside classification. A wrong-side interval
must be covered entirely by the permitted finite-edge/disk boundary capsules.
A midpoint close to a curb cannot excuse a larger excursion elsewhere.

Each proof invocation owns a bounded prepared point classifier for its unchanged
exterior and holes. Conservative y bins select possible ray-crossing edges and
boundary candidates; the original arithmetic, epsilon and nearest-edge tie order
remain unchanged. Continuous segment intersections, finite boundary capsules and
the width intervals still use their original algorithms. Other finite-containment
callers keep the original classifier unless they explicitly supply one.

The cache pins polygon buffers so fresh views share one preparation. It holds at
most sixteen polygon records and 8 MiB of prepared work arrays, with no eviction;
an index has at most 20,000 edges, 1,024 bins and 131,072 entries. Unsupported or
pathological indices permanently use the baseline for that proof. A new proof
after a boundary edit constructs a new cache; no cache survives globally or in a
published product.

The source-curb chord error and physical boundary band total at most 2 mm.
Known hole vertices are also tested against the swept source interior, so a
small hole cannot disappear between boundary samples or behind an artificial
sampling cross-section. Outer-boundary vertices receive the same test: a large
outside bay with a very narrow entrance cannot hide wholly inside one source
cell. Returning Exact sources retain their real central holes rather than
being treated as filled source polygons.

Partial Bundle certification still requires complete station ownership and
the actual exposed original half profiles. Whole Bundle claims additionally
include declared interior Sections in the width proof. Bundle centreline
separation, source-cycle construction, separator closure, mesh quality, UVs,
drainage and publication checks remain separate requirements.

The ownership-plan task declares Bundle source dependencies because interval
permission depends on member order, completeness and original terminals even
when the resolved crown values are unchanged. Surface, mouth and adjacent-seam
signatures retain the relevant local Bundle source digests. An unrelated Bundle
edit does not change a local junction's geometry signature.

The focused tests include missing midspan width, tiny interior holes, a bay with
a sub-millimetre neck, returning Exact geometry and disjoint prefixes, aggregate
cell-budget exhaustion, exact overlap/touch merging, and cached/fresh rejection of
an invalid Bundle member-order edit that leaves every crown fingerprint equal.
