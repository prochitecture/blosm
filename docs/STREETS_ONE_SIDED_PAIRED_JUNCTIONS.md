# One-sided junctions in overlapping paired pavement

The merged-entrance provider prepares these source junctions for multi-Section
Bundle runs before physical ownership or surface construction. The standalone
`OneSidedPairedSpecification.candidates` entry remains empty: the paired source
inventory discovers each original incidence once, together with its terminal
and shared original sampling. There is no second competing candidate provider.

An authored three-arm junction on one carriageway can overlap the other
carriageway even though that opposite Section has no source endpoint at the
junction. The original three-arm footprint and original source topology remain
authoritative. A nearest-node replacement or a new authored incidence would
change that topology and is not used here.

## Physical preparation

The caller supplies the original junction and prepared local footprint, its
two explicitly adjacent Bundle visits, the near-side role and the opposite
original Bundle visit. Finite transverse charts evaluate the original Paths;
they do not refit their polynomials. Original curbs are sampled with a bounded
source error and one shared original station partition where evaluations of
the opposite source overlap.

The footprint is the actual union of the local junction and three clipped
original full-width strips: its two near-side leads and the intervening
opposite strip. Two outer fibres expose the paired physical profiles. The
original near-side root profiles and finite strip root caps must be internal;
an exposed root cap rejects preparation. The stem keeps its original profile.

A proposed merged fibre requires actual pavement overlap. A bounded search
can find the first merge farther along the finite original source, followed
by the existing quarter-metre clearance. Any intervening grass gap remains an
explicit union hole. A source gap is never replaced with added pavement.

Continuous certificates retain the original local footprint and every clipped
full-width source polygon, including exclusion of holes contained entirely
within a source strip. The independent full-strip certificate's combined
source chord and boundary budget remains at most 2 mm. Initial derived curb
resampling permits at most 1 micrometre of
additional measured polyline projection; its error is deducted from that
combined budget. Original profiles, polylines, source spans and analytical
primitives remain available as provenance.

The normal surface builder selects `OneSidedQuadChart` only for explicitly
prepared groups with the complete source certificate. Its first attempt uses
the established 1 micrometre sampling. A typed `QuadChartGradingError` can then
select `provenance_rim_sampling`; source or topology errors cannot enable that
retry. The retry spends only the remaining authored curve tolerance after the
old analytical chord errors and earlier projection have been deducted. Its
construction projection is capped at 0.4 mm, and all projection is also
deducted from the separate 2 mm full-strip certificate. Repeated spending is
rejected. The ledger reports old local analytical chord bounds separately
from the original source-strip chord bound: they are different references.
All existing physical profiles, real hole tips and sharp source-rim vertices
remain exact. The old polylines, source spans and analytical primitives are
retained. Normal surfaces require coupled height and drainage reconciliation
and a final physical/stored pavement quality check after height correction.

## Original and derived aliases

`SourceMouth.original_incidence` defaults to `True`. A paired profile contains
one real near-side source alias and one opposite interior-cut alias marked
`False`. Both retain original global Path stations and half-profile parameter
maps. The false alias's direction is an orientation token, not an authored End.

`SeparatorPortIndex.resolve(End)` indexes only real source aliases. Interior
cuts use `resolve_derived(owner, section_id, direction)` and a separate registry.
`CompoundProduct.peaks` and `port_points` likewise retain only original
incidences. `derived_peaks` and `derived_port_points` include the physical owner
in their keys. Thus an interior cut cannot shadow a remote original endpoint
on the same Section.

`SourceEvaluationInterval(section_id, first, last, forward)` restricts metric
evaluation of an original Path. It creates no new Section. A specification's
entrances and opposite cut keys are ordered as the near member's preceding
and following visits. Complete original Bundle membership remains separate;
the local entrance deliberately has no fabricated `BundleRun`.

## Shared original curb sampling

Overlapping finite prefixes must share the same original Path object and the
same global station partition before any physical mouth is qualified.
`shared_source_parameters` pools the actual source stations;
`canonical_reprepare` builds fresh one-sided preparations from that pool.
Additional terminal entrances must also be rebuilt through
`MergedEntranceSpecification.build(source_parameters=...)`, rerunning their
original corners, source union and guards. The final overlap certificate
compares station arrays exactly. Accepted products are never patched or
snapped to a newly sampled body.

Supplied-knot sampling keeps already certified spans. It subdivides only spans
that exceed the existing original offset-curvature bound, retains Exact
control stations and rejects the unchanged 20,000-point allocation limit.
Repeated preparation is idempotent; it does not force two fresh intervals
between every existing sample. Numerical duplicate sampling stations retain
an existing representative before preparation. Ordinary entrance construction
without a supplied partition keeps its earlier sampling path.

`paired_junction_preparation.prepare(context, bundle, terminal)` provides the
explicit preparation handoff. It returns a fresh terminal specification,
one-sided specifications, pooled original metric stations and a joint phase
certificate. The terminal also carries `paired_sources` and
`source_phase_proof`. Discovery walks original consecutive member incidences;
it chooses the opposite original visit through a continuously certified common
frame and a unique complete source interval. It requires actual pavement
overlap beyond the chord uncertainty. Missing incidence, an interrupted inner
curb, ambiguous opposite intervals and unresolved near-contact remain explicit
errors. No nearest-coordinate connection is created. `MergedEntrancePreparation`
calls this API only for multi-Section runs; established single-Section
preparation remains unchanged. The returned terminal and side groups enter one
fresh ownership plan, followed by the normal junction surface factories.
Source-phase-aware local signatures invalidate related physical mouths when
the opposite member or shared source partition changes.

## Computational charts and acceptance

The optional paired cut policy joins existing merged-profile peak vertices.
With one real separator hole, it joins those vertices to the corresponding
existing hole tips instead. Original near-source normals can add bounded
computational bands, avoiding opposed duplicate cut families. Profile vertices
are shared after assembly; their geometry, source aliases and boundary rules
are unchanged.

Planar acceptance still requires positive direct quads, pavement edge ratio at
most 8 and corner sine at least 0.03, including actual stored float32 geometry.
Grass keeps its separate positive-corner and ratio-8 contract. Surface
qualification requires the original elevation envelope, physical slope limits,
exact source boundary heights and shared seams, zero interior minima, entrance
joins at most 3 percent, finite UV maps and no interior overlaps. A qualified
planar source union alone does not qualify a surface or the original Bundle
body.

The original-source regression exercises Exact and Smooth Paths, reversed
authored directions, rotation and translation, complete half-profile maps,
remote original endpoint ownership, the separate derived product registry,
real separator-hole preservation and computational cut containment.

The focused actual B20 preparation finds all four one-sided junctions and
rebuilds their shared original curb phases with the paired terminal. Five
normal Products qualify independently with 1,147 / 1,090 / 1,184 / 1,317 /
5,486 pavement quads for the terminal and four side groups. The last group
retains its real separator hole and 434 grass quads. The first three side
groups pass the strict sampling branch. Only the last group needs the typed
retry: its measured added projection is 0.3983 mm, and its combined original
local analytical bound is 3.5014 mm within the authored 5 mm tolerance. The
independent complete-strip certificate remains 2 mm. All five pass stored
pavement ratio/sine, source alias and peak bijections, exact physical boundary
heights, grass seams, finite UVs, zero interior minima, joins below 3 percent
and overlap checks. The complete source-owned Bundle plan also passes focused
qualification with 14,970 pavement and 2,430 grass quads, including both rounded
external caps, all original Section intervals and query-only junction meshes.
These results permit normal-provider integration; they do not certify a fresh
whole-city import or progressive publication.

## Complete mouths and remote source overlap

An exposed full transverse profile does not prove that the rest of its original
road is outside the owner's pavement. A different consumed source strip can
cross that outgoing road farther along its finite original route.
`source_mouth_clearance` checks every remaining full-width source cell against
the other actual prepared operands, including an operand contained entirely
inside a cell. Boundary-only contacts keep the old mouth. A positive overlap
advances only the affected original mouth beyond its last overlapping cell;
the complete source union, aliases, profiles, roots, neighbor ownership and
full-width proofs then run again. No source Path, radius or authored incidence
changes. Only overlap deeper than the two existing source-chord allowances
plus the unchanged numerical boundary band can request movement. An exhausted
ordinary full-profile search can reuse this actual-owned-overlap proof only
to request its recorded finite original neighbor; an alternative station or a
proof-only future strip cannot authorize ownership. Scope and selection ends
remain failures. The search shares the existing 4,096 finite candidate budget, with
at most 20,000 cells and 2,000,000 broad-phase cell/operand pairs per pass.

The captured Bergamo Bundle source case has an already complete road mouth,
followed by a separate overlap with an original consumed crossing strip.
Replaying those actual operands advances the original mouth from 1.0152 m to
41.2639 m along its unchanged 199.0768 m route. A fresh current source rebuild
already places that mouth at 48.7020 m, so the added check leaves this successful
result unchanged. Its three actual holes retain byte-exact boundaries, and an
independent 2 mm complete-width certificate covers 2,338 source cells including
the extended original prefix. Surface and complete Bundle-body qualification
remain separate requirements.

## Source-cycle roots discovered during ownership

Some complete Bundle runs become paired entrances only while ordinary source
cycles are resolved, rather than in the earlier merged-entrance task. The
`source_companions` preparation runs after those source groups are built and
before any surface Product. It matches a root's **complete ordered original
run** against the declared Bundle, discovers only genuine intermediate source
incidences, and canonically prepares an actually overlapping opposite member.
The existing root is immutable. Exact shared Path identities and station arrays
must agree; a phase change requiring another root construction remains an
error, rather than permission to snap a mouth.

An explicit `PreparedSourceRun` records the actual root group, one complete
run, real side groups, joint source-phase proof and original source IDs. It
does not invent a crossing record or opposite Section endpoint. The ownership
plan installs the actual side junction once and exposes the resource through
`bundle_sources`. A one-run consumer must still cover every original member
Section exactly once and every selected source interval completely. Ordinary
two-run resources retain their existing default and coverage checks.

The source-plan task declares original Path/Profile dependencies. Each derived
side layout/surface signature also includes the opposite source paths and
profiles absent from its original local footprint. Moving or widening that
opposite member therefore invalidates the side's shared geometry. No cached
road mesh supplies an undeclared input. Existing early paired preparations are
excluded to prevent duplicate ownership.

The current Beijing source-task replay creates 43 source groups over 100
original junction sites and one new complete-run resource. Its actual side
junction has one owner. Fourteen focused source, phase, input-dependency and
terminal-prefix regressions pass. This source-stage witness does not qualify
the complete Bundle surface, fresh import or progressive publication.

## Finite merged terminal and separated remainder

`bundle_terminal_unions.TerminalMergedBundlePlan` handles a different case:
an actual paired terminal followed by an ordinary separated run. After the
existing Bundle plan families decline, it resolves both true original terminal
Ends through the physical port registry and checks the complete run identity.
Only the finite overlapping prefix is generated as one merged surface. A
private evaluation plane cuts each unchanged original source beyond its last
overlap, with a complete exact carriageway profile on both sides. It stops
before the next real side junction or incompatible/dimensional seam.

The remaining road pieces use the ordinary carriageway builder over explicit
subsets of their accepted visible source intervals. Real junction caps and
source roads retain their own owners. The assembled prefix and remainder
share complete profile and grass edges; no evaluation cut becomes a physical
End, junction, SourceMouth or registry entry. Original absolute longitudinal
UV stations are retained. Every original source interval is accounted for,
and final source containment, stored shape, height/drainage, longitudinal
joins and global overlaps remain mandatory. Only a qualified **complete**
replacement allows member roads to be suppressed.

The Bundle task declares the actual attachment layouts, crossroad Path/Profile
inputs and physical-port readiness needed by this construction. A focused
normal task entry point currently qualifies 20,529 pavement and 3,169 grass
quads in the real terminal case; fresh whole-import qualification is pending.
