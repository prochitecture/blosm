# Compound junctions and short source prefixes

This document describes the Python backend's additional preparation branch for
connected junctions that cannot first be constructed as independent surfaces.
It does not change the original OSM connectivity, Section records, widths or
profile rises. Current real-network qualification is recorded separately in
[the Beijing and Bergamo report](STREETS_BEIJING_BERGAMO_FIXES.md).

## Deferred preparation

An independent circular corner may exhaust a short Section before reaching the
neighboring junction. `DeferredCompoundFootprint` records the original
approaches, required physical neighbors and diagnostic text. It contains no
provisional polygon or mesh. `OwnershipPlan` resolves these records before any
road or junction surface can be published.

`DeferredCompoundPlanner` follows explicit opposite source incidences. It is
bounded to eight member junctions and eight connected expansion passes. A
candidate's final exterior must retain complete external profiles, be one simple
boundary, and contain the complete consumed source centerlines with certified
clearance. Every consumed internal route must have both original incidences.
Different elevation datums, internal dimensional seams and consumed Bundle
bodies still require other policies; this branch rejects them.

If an external profile already overlaps the opposite source claim of another
physical junction, that ownership is resolved before the provisional boundary
is tested. No provisional faces are generated or accepted during expansion.

## Cyclic order at distinct origins

Final approach tangents cannot determine cyclic order when their mouths start at
different junctions. In particular, two nearly parallel Bundle arms can swap and
connect their outer curbs instead of the facing inner curbs.

Each original junction provides a cyclic incidence order. The planner walks
that order and crosses each wholly consumed internal route to its opposite
incidence. The resulting external cycle preserves source adjacency without
inventing a common junction center. All external incidences and both directions
of every internal route must occur exactly once. Additional internal cycles are
rejected by this branch: a real hole must not silently become a filled disk.

`ExteriorSpecification.ordered_ends` carries the derived ordering to the contact
solver. Ordinary common-origin intersections retain tangent-angle ordering.

## Short source prefixes

A very short Section between a physical junction and an explicit dimensional
seam may not support the tangent imposed by the next Section. `OwnedPrefix`
allows a narrowly defined alternative after the ordinary path construction
fails:

1. Certify the original prefix curve using its unforced endpoint tangent.
2. Use the explicitly connected next Section's actual frame and complete profile
   as the candidate exterior approach.
3. Preserve aliases for the consumed prefix and the visible external mouth.
4. Prove that the entire prefix and absorbed seam lie inside the final junction.
5. Reject publication if any part of the prefix remains an independent visible
   road interval.

This does not approve an invalid visible road or modify source control points.
An external Bundle half-profile normally belongs to a complete original member
terminal. The separate contiguous-interval policy below also permits explicitly
owned internal extents. The same Section/Street IDs continue to identify updates
and furniture anchors.

## Contiguous internal Bundle intervals

The end of a compound can lie partway along an original Bundle. A local source
cycle can contain one complete Section of a lateral carriageway while the rest
of that member, the opposite member and most of the separator remain visible.
Consuming the entire Bundle to avoid that case would discard real source areas.

`partial_bundle_claims.py` records explicit, bounded alternatives. Each requires
two complete original lateral members and contiguous consumed source intervals:

- An internal interval retains both adjacent source ends and half profiles.
- A paired HEAD or TAIL interval includes the corresponding original terminal
  junction of each member and retains both surviving half-profile cuts.
- An asymmetric terminal interval can retain an untouched member's original
  terminal half profile while consuming the corresponding other member's end.
  Both original terminal owners must be in the compound; a retained Section
  whose opposite terminal is also owned must instead be explicitly claimed.
- An entirely consumed pair requires all four original terminal owners. It is
  admitted only by the physical source-cycle boundary policy, which computes
  remaining holes from the union of actual pavement; a single exterior disk may
  not simply fill its separator. Original curve separation is still certified.

Original member traversals identify the surviving ends. Each must belong to an
actual junction in this compound and retain its complete, unchanged half
profile. Disconnected intervals, missing corresponding terminal ownership,
absorbed profile changes and missing original incidences are rejected.

The physical compound must then prove complete coverage of every claimed
Section, across its complete width and not just its centerline. Convex source
cells must lie inside the actual exterior and outside every real hole, including
holes contained entirely in a cell. The existing source-curve approximation
budget is retained. Its actual sampled exterior must expose both adjacent half profiles;
an interval cannot hide those ports or silently consume another member Section.
The output records the original Bundle/member/Section IDs and cut ends in
`partial_bundle_claims` metadata. No source points or connectivity are changed.

Bundle membership and member order are declared dependencies of ownership,
mouth and seam products. Reordering a member can change containment without
changing a Section's crown; unchanged profile fingerprints must not permit a
stale accepted plan. The cache tests compare fresh and cached rejection and
verify that an unrelated Bundle does not invalidate this local signature.

This permission does not approve a separator. The established visible-run
builder must still resolve the remaining member intervals, use accepted physical
junction curb chains for their end caps, certify centerline separation and pass
the final pavement/grass overlap tests. An intervening road mouth can split a
median curb into separate runs; it is never bridged with an invented straight
edge. Whole-Bundle closed source cycles retain their separate physical-hole
policy.

## Paired median end caps

A cap between converging inner curbs may need slightly more than a semicircle.
This is explicitly distinguished from an ordinary turn fillet. The compound
planner identifies adjacent complete terminals of the same Bundle, at the same
Bundle end, with half-profile peaks facing each other.

Only that paired-terminal primitive permits a clockwise circular sweep up to
210 degrees. Contacts must lie on the finite source curves and remain within
`2 × maximum carriageway width + distance between mouth origins`. Radius,
tangency, complete profile, simple-boundary and all surface qualification checks
remain mandatory. Explicit radii remain strict. Ordinary corners retain their
180-degree sector limit. The fitted-radius report identifies a
`PAIRED_MEDIAN_CAP` when this primitive reduces a style radius.

## Persistent construction overlap

A longer pair of approaches can remain too close for even the smallest permitted
fillet up to the next physical junction. `overlapping_approaches.py` supplies a
bounded candidacy test, with at most 4,096 centerline samples per route.

For a chord of the first route, let `d0` and `d1` be endpoint distances to the
second sampled route, `L` its chord length, and `e0`, `e1` conservative curve
approximation errors. The whole chord is sufficiently close only when

```text
max(d0, d1) + L/2 + e0 + e1 < (width0 + width1)/2 + clearance
```

The bound follows from the distance function being 1-Lipschitz. For physical
overlap, clearance is zero. For a required circular fillet, it is twice the
smallest permitted radius. This latter condition can include a narrow positive
gap; it authorizes another planning attempt, never filling an unowned island.

After a member is wholly owned, its original centerline can provide the next
source hop to another exterior approach. The same certificate may then identify
one further real terminal junction. This does not infer a connection from
proximity. Selection/scope boundaries, incompatible profiles, unbounded work and
invented terminal junctions are rejected.

## Finite contact searches through owned source nodes

A short internal Section can end before the fillet between two roads reaches
its physical contact. `compound_contact_routes.py` can continue the contact
search through an already owned original junction. This is a piecewise search
support, not a smooth road seam: original centerlines and offset frames remain
unchanged at a bent join.

The next original incidence must be unique, turn at most 20 degrees, and retain
the same width, oriented profile, offset and interpolation mode. Its Street ID
can differ because an ordinary cross street can split a physical through route
into separate Streets. Coincident coordinates do not supply incidence. Scope
and selection ends, ambiguous forks, changed dimensions and absorbed profiles
stop this branch. Searches use at most 32 source pieces and stop after entering
the first exterior Section.

Contacts must lie away from frame discontinuities. Every source alias maps back
to its original station, and the complete external mouth moves beyond retained
contacts on that source. Patch simplicity, source-width containment, ownership,
unexposed preparation seams and later surface checks remain mandatory.

A curved near-through pair can also need the other contact-normal side from
the one suggested by its tangents at the original node. After ordinary STYLE
fitting fails, `compound_corner_contacts.py` tries that side only within 15
degrees of a through alignment, at the same declared source node and using the
unchanged permitted radius sequence. Contacts remain within twice the larger
width; analytical radius, tangency, sector and negative-extension checks remain
in force. Explicit STRICT overrides do not enable this secondary branch.

This secondary candidate separately audits its facing root halves: the first
road's left half and the second road's right half. Other halves belong to other
local sectors and still need the complete source-union proof. This does not
relax the ordinary whole-cap audit. Fitted-radius metadata records the branch,
source route pieces and owned junctions used by the search.

## Exposed roots and partial paired terminals

An internal corner is not retried merely because its first circular contact
fails. The unchanged source-footprint union must first prove that the local
corner is hidden. Only a typed exposed-root result authorizes the bounded
alternative contact search. This ordering prevents an unnecessary internal
fillet from adding pavement or removing a real separator hole.

When a consumed interval reaches the HEAD or TAIL of one complete Bundle
member, its paired terminal can still belong to a neighboring original
intersection. `bundle_partial_terminal_expansion.py` acquires that intersection
only when an explicit original source route joins it to the compound. This
keeps the matching half profiles under actual physical ownership. A same-XY
endpoint is insufficient. Scope cuts, disconnected intervals, missing owners
and noncontiguous partial ownership remain errors; the ordinary complete-width,
contact and partial-profile checks run again after expansion.

A terminal closure can cross an intervening real physical junction. The search
then acquires **every** intermediate owner, rather than jumping over it. This
extension uses at most four original route hops and 64 search states, excludes
scope/selection handoffs, requires one shortest owner closure, and bounds total
original-path length by the larger of four member widths or two declared corner
radii plus one member width. It authorizes another full construction attempt;
the existing eight-owner budget and complete source/profile checks still apply.

Interior physical cuts needed for overlapping multi-Section Bundle bodies are
evaluation intervals, not new authored Sections or OSM traffic nodes. A
`SourceMouth` with `original_incidence=False` is stored in the separate derived
port registry. It may not replace the original endpoint's `resolve(End)` binding.
The new paired side-junction construction remains disabled until its complete
source, profile, surface and body checks pass.

## Publication precision and immutable profile slopes

The publication guard still accepts only the actual float32 mesh under the
original physical slope limit. Its correction target must also be attainable
without changing prescribed boundary heights. For each quad corner, let
$G$ map its four heights to the two-dimensional gradient and let $G_f$
retain the columns belonging to movable vertices. The component of $Gz$
orthogonal to the range of $G_f$ is immutable. Its norm provides a local
lower bound on every feasible slope target.

`PublicationSlopeGuard.fixed_gradient_floor` computes the largest such bound.
The established roundoff target remains first. If its global safety margin
would push the target below the fixed-gradient floor $f$, the guard instead
tries $f+0.75(L-f)$, where $L$ is the unchanged acceptance limit. This is
a search target, not a relaxed acceptance tolerance. Actual stored slopes,
unchanged XY/topology/profile rows, drainage and incoming longitudinal slopes
must still pass. A fixed boundary that itself exceeds the published limit is
rejected.

An inner numerical target may also exhaust its iteration budget even though the
last iterate already meets the actual surface contract. The default height
solver still raises `ProjectionIterationBudgetError`; it carries an explicitly
unaccepted candidate. The publication guard alone may certify that candidate
using the unchanged physical and actual float32 slope limits, height box and
spill floors, drainage, incoming joins, exact XY, exact boundary and topology.
Failure of any check remains a failure. This does not enlarge a slope tolerance
or accept a failed height solve elsewhere. Metadata distinguishes such a fully
checked candidate from convergence to the stricter optimization target.

Inside the coupled spill-floor loop, the same typed iteration-budget exception
may supply an intermediate physically feasible height iterate. It must preserve
XY, topology and fixed XYZ, respect the current lower spill floors and crown
ceiling, and pass the existing physical slope tolerance. It is not a published
result: the remaining bounded drainage rounds, incoming-profile joins and actual
stored-coordinate checks still have to pass. An intermediate basin or float32
slope error therefore cannot disappear from the final acceptance test merely
because the physical height iterate was useful to the next round.

If every established ordinary chart/refinement attempt also fails an entrance
join, the late coupled recovery can add explicit directional height slabs.
At fixed XY, a bilinear corner gradient is linear in its four vertex heights:
$g=A z$. For each complete original incoming profile edge, both endpoint
gradients are constrained by $|n^T A z|\leq j$, where $n$ is the profile's
unit normal. Convex interpolation then bounds that normal gradient along the
whole edge. The construction uses $j=0.95\,j_{max}$ to reserve numerical room;
the original three-percent acceptance limit stays unchanged.

`HeightSlabs` projects only the free height components. For a row $a$, the
local update is

$$
z' = z +
\frac{\operatorname{clip}(a z,-j,j)-a z}{\lVert a_{free}\rVert^2}
a_{free}.
$$

A zero free row that violates its bound is infeasible and rejects the chart.
Dykstra corrections combine these slabs with the existing physical slope balls,
crown box and monotone spill floors. Row count, colors, iterations and memory
remain bounded. The default solver receives no slabs and retains its previous
arithmetic. The recovery keeps every fixed XYZ and all XY/topology exact;
drainage, actual entrance checks and stored-coordinate publication checks still
follow. It does not make an impossible profile or unsupported footprint valid.

## Finite Bundle body pieces

`MergedRunBody.build_between` reuses the maintained direct-quad body construction
between two prepared common transverse fibres. The four inputs are the head
and tail half profile of each of the two members, each ordered outer curb to
crest. A fibre can contain two separated peaks or one shared crest. The same
core therefore supports two merged ends, two separated ends, and a change
between them with an automatically terminating separator.

`CommonFibreCuts` checks complete profiles, actual curb/crown endpoints, finite
source-curve domains, ordering and the prescribed dynamic height field. It
never chooses a physical owner from coordinates. A crest shared by the two
already-declared halves must be one identical authoritative sample; the core
uses a shared assembly key. All supplied profile vertices stay exact. Generated
grass and pavement share complete edges, including the profile peaks. Separate
pieces use separate key namespaces.

When all established open-island choices fail, `OpenFibreIsland` can retain
the complete original banks and add only collinear samples on their existing
chords. Its two half-strips use cached direct 1:3 transitions and a finite
width-aware row partition. Protected bank vertices and the true closing tip
are fixed. The exported bank samples and keys are then the pavement's exact
whole-edge interface; the grass is not independently snapped onto a completed
pavement mesh.

`FixedRimBodyRecovery` is a corresponding failed-only pavement chart family.
It retains the existing semantic curb/profile/crest topology and every physical
rim point. Bounded fixed-rim shape-energy seeds precede the existing direct
quad rotations and anisotropic interior refinement. This can improve the
computational layout without copying every tiny grass-bank interval across the
full carriageway width. The original dynamic height field is sampled again;
complete boundary heights stay prescribed. Physical and stored 8/0.03 pavement
grading, crown bounds, physical/longitudinal slopes, drainage and overlap are
checked independently before the result is usable.

For a few inverted computational seed cells, `quad_seed_untangling.positive_seeds`
may move one interior vertex while keeping the complete rim and connectivity
fixed. For that single variable vertex, every corner Jacobian in its complete
incident one-ring is an affine function of its XY displacement. Intersecting
those positive half-planes gives a bounded feasible seed region. At most four
nonpositive cells and four candidate seeds are considered, with displacement
limited by the original one-ring's extent. The complete double-precision and
stored meshes must be positive. This is a seed certificate, never a final
quality certificate: all pavement grading, source, height, join and overlap
checks still follow.

Stored-coordinate cell ratio remains at most eight, pavement corner sine at
least 0.03, physical slope at most 0.1, and carriageway longitudinal slope at
most 0.03. Grass retains its established positive-corner and ratio-eight policy;
an original narrow grass tip can make the pavement's angle threshold impossible
without changing its physical banks. Actual grass angles remain visible in the
quality report. A boundary-spill
audit rejects interior depressions. The old one-merged-end entry point retains
its established interface and assembly keys.

The multi-Section `MergedBundlePlan` prepares every actual side junction before
meshing, splits only evaluated source intervals, transports exact physical
mouths to the common-fibre body cuts, and accounts for every original source
interval. Junction areas are included in the height query once and are not
published again by the body. Focused whole-Bundle source coverage, aliases,
UVs, queries, overlap and release checks pass. Fresh native and progressive
publication remain the separate complete-import acceptance gate.

`TerminalMergedBundlePlan` handles a resolved paired terminal whose finite
overlap ends before the next physical side junction. Original Section
cross-sections delimit the replacement prefix. They are computational planes,
with no new End, opposite alias or physical owner. The complete remainder is
rebuilt from the normal original Path/Profile road builder, and its separator
follows actual intervening junction curbs. Whole profile edges and the grass
interface are welded exactly. Every original source interval is accounted for
by the prefix, ordinary remainder or an already accepted physical junction.
The ordinary builder's optional `source_interval` accepts only a positive subset
of its normal visible extent and retains absolute source-station UVs.

A finite terminal bridge may serve an exact terminal-incident original prefix
or suffix of its complete run. This requires explicit ordered source visits,
unchanged membership/direction, actual consecutive incidences and the same
original terminal. Its progress proof still covers every original run Path;
the final assembly separately proves ownership of all remaining intervals.

The established run coordinate frame stays first. When it cannot support all
finite mouths, `bundle_mouth_frames` can intersect two actual prepared mouth
normal lines to choose an angular coordinate pole. All complete original
offsets must certify positive progress in that chart. This changes coordinates,
not source geometry, physical mouths or connectivity. The resulting meshes
still pass all profile, source, stored grading, height and overlap guards.

Only an explicit combined-cell failure admits `BoundedMergedRunBody`'s two-body
family. One strictly merged common fibre divides computation; both halves use
the existing body builder with its eight-thousand-cell limit. The shared crown
and every transverse edge are exact. The complete leaf, including its
connectors, retains the sixteen-thousand-cell limit. Physical rims, source
intervals and genuine grass holes stay unchanged. Rounded grass caps also use
the existing failed-only publication slope correction when float32 storage
violates their ten-percent physical limit; all XY, quads and boundary heights
remain fixed.

### Physical query geometry and inspection exaggeration

Compound junction products retain authoritative `physical_pavement` and
`physical_grass` mesh views from before display exaggeration. At exaggeration
one these share the published geometry buffers. Other display settings preserve
the physical coordinates separately, including exaggeration zero, where the
display mesh alone cannot recover elevations. Bundle preparation therefore does
not infer grass heights from a pavement-only query or invert an erased Z value.

`physical_context` supplies temporary product/Part views over these exact meshes
and the existing physical rims. Source geometry and published Parts are not
mutated. The meshes remain construction resources: `release_geometry` clears
both physical views along with published buffers, while the compact attachment
height query remains available. Any later construction needing the meshes must
use the normal hydration path. `required_physical_owners` resolves only the
Bundle's original incidence owners in the current plan and its explicit paired
preparations. Unrelated released junctions are skipped. A missing required mesh
is rebuilt through the normal junction builder in a temporary physical-height
context; accepted cached Products are not mutated. This also reconstructs
standalone band junctions at physical heights rather than inverting erased E0
coordinates. Regression tests compare both pavement and grass byte-for-byte
across exaggeration zero, one and ten, exercise release, and compare a normal
released/edit/restore sequence with fresh builds.

### Finite terminal cuts and consumed terminal links

A paired terminal is selected at the actual cut required by its original
corners. Overlap at a nearer provisional cut is insufficient. If that required
cut has separated original carriageways, `cycle_terminal_ports` may retain the
two original full mouths instead. Its certificate inverts all four analytical
offset curbs inside their finite original ranges, proves positive axial progress
and regular offsets, and requires a gap larger than both source sampling
allowances. The successful merged-cut path remains first; uncertain or touching
cuts do not establish separation. No source End, half profile or gap height is
invented.

When a failed local STYLE contact faces a fully consumed full-profile road, a
circle may reach the opposite terminal member of the same declared Bundle
through that road's actual far junction. `bundle_terminal_contact_walks` admits
only the immediately adjacent original bank at an owned junction, equal actual
widths, unchanged matching positive half-profile peaks and the original bounded
STYLE radius sequence. The circle closes through the consumed road's original
centreline. Each contacted member keeps its own Path and Profile; a uniform
half-profile Route never replaces the full-profile source road. Complete source
width, mouth, internal-seam and subsequent mesh checks remain mandatory.

An unresolved corner with one consumed arm can be omitted only if the complete
consumed facing curb and the external arm's actual mouth-bounded owned prefix
stay continuously in the pavement and have no positive-length coincidence with
any exterior or hole boundary. Mere overlap near the original node, an unrelated
hole, or a far exposed lead does not qualify this proof. The established
two-consumed-arm hidden-corner policy remains separate.

## Validation and entry points

The mandatory mesh checks still apply after preparation: positive quad
Jacobians, edge ratio at most eight at publication precision, profile and height
constraints, longitudinal slopes, source seams and forbidden surface overlap.
A qualified footprint is not a qualified mesh.

Relevant modules are `compound_preparation.py`, `owned_prefixes.py`,
`overlapping_approaches.py`, `partial_bundle_claims.py`, `compound_coverage.py`, `junction_transitions.py`
and `ownership.py` under `way/street_network/python_backend`.

Focused tests are `test_python_street_deferred_compounds.py`,
`test_python_street_owned_prefixes.py`, `test_python_street_compound_order.py`
and `test_python_street_overlap_candidates.py`, plus
`test_python_street_partial_bundle_claims.py`. Run them with the existing
ownership, source-support, compound-port, corner-radius and decomposition tests.
These numerical tests do not require starting Blender. Actual publication still
needs the isolated Blender integration tests described in the main report.
