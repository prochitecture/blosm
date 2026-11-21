# Finite ground scope ports

Two ground approaches can share an original junction yet remain physically
overlapping all the way to two distinct recorded scope endpoints. A circular
internal gore then cannot fit inside the retained source lengths. Extending
either ground road past its source endpoint would invent supported geometry.
Treating the endpoints as free ends would erase the scope handoff.

`finite_scope_ports.py` provides a bounded physical-port policy for this case.
It currently accepts two straight, unbundled, full-profile ground approaches
with a shared original junction and separately recorded far scope incidences.
Their original finite corridors must overlap continuously in both directions.
Selection cuts, unrecorded/free endpoints, Bundle profiles, absorbed seams and
curved tails do not qualify. A strict requested internal radius remains strict;
the automatic policy applies only to a failed style-radius corner.

The source Sections, source points, original junction and both scope incidences
remain unchanged. The exact union of the two finite source strips supplies the
distal XY outline. The two original outer curb fillets retain their circular
contacts and radii. Their curb leads end at the original outer scope corners.
No coordinate outside the finite union is added to the distal outline, and the
new port rejects longitudinal extrapolation.

Overlapping full crowns cannot both remain independent boundary conditions.
The new physical scope mouth therefore has an explicit **reconciled complete
crown**. This is a new derived boundary profile, not two claimed original full
profiles. Its datum is the common original datum and its ceiling is the highest
original incoming crown. Both halves are monotone, with zero height at the
outer curb joins and one complete crest. Quintic pieces have stationary jets at
the actual finite-outline bends, giving a common height and gradient at those
bends. The ordinary, continuing incoming road keeps its complete original
profile and source-path mouth.

Each source alias records its original junction end, its distinct scope end,
the original scope station and the complete consumed source interval. These
records are published in `finite_scope_ports` metadata alongside the mode and
finite outline. They do not relabel an elevated incidence or generate elevated
road geometry at ground level.

Before surface construction, one simple physical exterior must contain every
complete consumed finite source-strip footprint. Straight source segments are
split at every physical-boundary intersection; every resulting open interval
is classified. This continuous proof detects exits shorter than any fixed
sampling interval. Original terminal points may lie on the explicit open scope
port. The new mouth sampler retains every finite-outline bend and its complete
crest.

The ordinary direct quad charts, fixed-rim quality correction, overlap checks,
UV mapping and stored-coordinate publication checks remain active. Finite
source pools additionally require drainage and entrance-join reconciliation
even when the initial field already meets the physical slope limit. Their rim
XYZ stays fixed, the quad edge-ratio limit remains eight, and the entrance
longitudinal slope limit remains 0.03. A focused passing construction is not a
whole-import coverage claim.

`test_python_street_finite_scope_ports.py` checks complete mesh preparation,
source and scope identity, outer fillets, complete profiles, UVs, stored float32
quad quality/slopes, reversed directions, rigid transforms, resolution changes,
rejection of unsupported boundary meanings and continuous containment of a
segment crossing a sub-sample notch.

## Open scope spur beside a continuing ground road

`finite_scope_spurs.py` handles a different finite configuration. A narrow
straight ground spur has a recorded far grade-scope endpoint, while a wider
ground road continues through the original shared junction. The two facing
curbs have only a small gap at the scope plane. Even the smallest permitted
style fillet would need a contact beyond the retained ground spur. The omitted
elevated continuation supplies neither a tangent extension nor pavement.

The continuing wide road keeps a complete original profile farther along its
real, finite ground route. Its original facing curb returns to the spur's
exact scope plane. The finite closing edge across the gap belongs to an
explicit `OPEN_SCOPE_BOUNDARY`; it does not claim to be an internal circular
gore or a closed dead end. The spur's scope mouth extends along that same cut
plane to the wide curb. Its original outer crown half and crest remain exact;
its inner crown half is reconciled to a zero jet at the wide-curb endpoint.
The reconciled half cannot exceed the original spur's maximum crossfall or
height. The other original outer corners retain their circular fillets.

The complete original spur footprint must remain contained. The wide ground
prefix contributes its original finite strip to the physical union, so a
negative auxiliary lead from another source road cannot silently cut inside
that strip. Shared original curb portions use one common chord partition.
Every exposed incoming mouth must remain complete, and an internal support
cross-section may not become a physical curb. Beyond the spur's scope plane,
the resulting surface is restricted to the original continuing ground strip.
Smooth source-curb chords retain a certified sub-millimetre sampling error and
the ordinary published curve tolerance; this is not endpoint extrapolation.

This narrow policy currently requires a three-incidence ground junction,
unbundled full profiles, style radii, one straight scope spur, a continuing
road within 30 degrees of that spur, and an approximately opposite third
ground approach. Selection cuts, true free ends and strict-radius failures do
not qualify. Unsupported Exact miter policies remain rejected. Source records,
scope aliases and the wide route's complete-profile stations are reported.

`test_python_street_finite_scope_spurs.py` checks the full mesh with a compatible
source seam, exact outer crown preservation, the continuing original profile,
scope identity, source reversals, reflection/rotation, a straight Exact case,
UVs, stored float32 quads/slopes/joins, local updates after buffer release and
transactional rejection when the source boundary meaning changes.

## Composing an adjacent original junction

An already owned open scope pool remains a physical neighbor when an adjacent
source junction must defer. `finite_scope_compounds` composes the two original
three-way sites only when they retain both incidences of one complete ordinary
ground connecting Section. The source rotation determines exterior order;
sorting displaced mouth tangents cannot replace that incidence walk.

If the adjacent continuing road still cannot make a permitted finite fillet
with the original scope spur, the same recorded scope plane can close the
remaining ground gap between the two actual facing curbs. The new outgoing
profile explicitly reconciles both halves, retains the original scope crest
position and bounded rise, and has zero jets at the outer curb joins. The
original scope Section, width, terminal plane and source alias remain recorded.
The two continuing roads retain complete original profiles farther along their
finite routes, before any opposite source-junction claim.

The complete original connecting strip, scope strip and any wholly consumed
continuing Sections require continuous full-width containment. Finite original
prefix strips provide the remaining ground support. Every non-through source
polygon must stay on the recorded ground side of the scope plane; only the two
known original continuing strips may cross it. Internal source cross-sections
cannot become exposed curbs. Bundle, selection, absorbed-transition, strict
radius and ambiguous-neighbor cases require a separate policy.

Two computational cuts on that scope plane separate the retained continuing
tails from the shared pool. Positive fixed-rim leaf atlases can use bounded
local collars around otherwise unrepairable unequal corner edges, preserving
all physical boundary samples and their profile parameters. The assembled
quad, stored-coordinate, height, slope, join and drainage guards remain
mandatory. `test_python_street_finite_scope_compounds.py` checks the original
incidences, complete strip proof, real opposite mouth, scope identity, source
storage reversal, reflection/rotation and unsupported-domain rejection.

The derived union can inherit sub-millimetre clipping pieces from a nearly
straight original fillet. These are computational samples. Composed scope
curbs retain the original union polyline, source spans, source curve references
and original analytical corner primitives, and may resample interior curb
points within an additional 1 micrometre projection allowance. The measured
projection error is reported against the saved **polyline**; analytical
equivalence is not claimed. All full profiles, scope endpoints and visible
corners remain fixed, and complete source-width certification runs again.
