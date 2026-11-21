# Finite wide free caps

A short, wide authored free carriageway can be wholly owned by its original
junction while two other ground approaches remain external. A two-curb through
pool can fail when those external roads turn significantly or have different
widths. The full original pad footprint still needs physical support.

`finite_free_caps.py` provides a separate bounded policy for this configuration.
It requires three original incidences, a uniquely short straight true free pad,
unbundled full profiles, common datum and style radii. Scope boundaries,
selection cuts, strict-radius failures and absorbed dimensional reserves remain
outside this policy. The source Sections and original endpoint identities are
unchanged.

The two real incoming roads keep complete original profiles at finite stations
on their source paths. Their original exposed outer fillet is solved and
retained. The free pad gains a separate zero-height curb side on its exact
original terminal plane. Its lateral ends extend far enough to contain the
complete original pad corners and allow rounded connections to the external
roads. The straight cap remains a curb boundary; it imposes no independent
incoming carriageway profile.

Two quintic curb sides connect those real mouths to the terminal curb with
matching tangents. Their endpoint tangents follow the cap line, and their
complete Bernstein curvature certificates respect the fitted style radius.
Preferred radii are tried first, down to the unchanged 0.5 m style minimum.
The requested radius, attempts, actual radius and certified bounds are reported.
The cap-side polynomial on the pad's outer side stays within the original
terminal plane. No source centreline or source endpoint is extrapolated.

The physical preparation unions these rounded boundaries with both original
finite incoming source strips. This retains a curved source curb wherever a
new cap connector lies inside the original carriageway. Original source-curb
samples and circle contacts use a common partition. Bounded prepared curb
pieces preserve all retained physical corners and source approximation budgets.

Before consolidation into one physical exterior, both original incoming
profiles must survive as separate, contiguous, complete parameter runs. An
exposed internal preparatory cross-section is rejected. The resulting physical
owner retains those original two mouths and the complete curb chains between
them. The complete original pad rectangle must remain contained. Neither a
partial profile nor a disconnected/holey physical exterior qualifies.

The ordinary quad edge-ratio limit of eight, original incoming-profile slope
limits, UV mapping and float32 publication checks remain active. Finite source
pools additionally require exact-rim drainage and entrance-join reconciliation.
The interior field can reconcile the former pad crown; all real incoming
profiles and outer zero-height curbs stay fixed.

`test_python_street_finite_free_caps.py` checks the complete mesh, preservation
of original source records and profiles, the exact zero-height terminal plane,
curvature certificates, reversed and rigidly transformed source geometry,
scope/selection/strict rejection, UVs, stored float32 quads/slopes and incoming
joins. Focused qualification does not replace native whole-import acceptance.

## Partial Bundle closure on a published physical curb

`finite_free_closure.py` handles a different configuration: one member of a
partial Bundle has a genuine free end, while the opposite member's rounded
terminal curb continues along an independent original road. The existing walk
to a real free or grade-scope end runs first. A typed `PhysicalJunctionStop`
alone admits the additional policy; selection cuts and unrelated errors do not.

The continuation keeps every actual road-bank edge up to the stopping junction.
The declared bank side may then follow exactly one published junction-curb
chain. That chain comes from its accepted physical product and ends at the next
actual mouth. It does not connect, traverse, consume, or reclassify the next
road. A unique transverse crossing of the original free profile's terminal
plane selects the existing inside curb vertex, strictly before that next mouth.
The small inward offset is bounded by the one published crossing chord. Tangent,
coplanar, missing or multiple crossings reject the closure.

A new computational exterior grass connector joins this vertex to the unchanged
free member's original inner peak. It creates no source incidence or physical
mouth. If the complete grass rim has odd parity, this connector alone permits
one collinear sample; no pavement edge is subdivided. Both the published curb
prefix and the independent roads retain their original whole edges.

After assembly, this new typed closure uses fixed-boundary height projection,
boundary drainage and the stored-buffer slope guard. All XY positions, quads,
actual curbs, incoming profile heights and connector heights stay fixed. The
rounded grass slope limit is ten percent; grass corners remain positive with
edge ratio at most eight. Non-unit exaggeration reconstructs a temporary normal
physical plan and requires an identical chart. In particular, E0 never infers
physical heights from a flat display mesh. Height queries use the reconciled
grass, including after publication geometry is released.

`test_python_street_finite_free_closure.py` checks plane uniqueness, original
vertex selection, the single computational parity sample, fixed-boundary
feasibility and physical reconstruction at E0/E10. Its original airport-control
fixture retains the independent next road and the matching whole pavement/grass
edges. The established partial-cap tests explicitly forbid entering the new
policy when the original terminal walk succeeds. Source/task fingerprints also
include the extra stopping junction's actual physical curb.
