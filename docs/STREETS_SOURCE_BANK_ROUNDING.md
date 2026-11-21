# Finite rounding between original source banks

An original-source union can expose a reflex switch from one carriageway bank
to another. Both banks can have zero height but different crossfalls. Joining
their existing height rules at that point cannot produce a smooth surface.
A small circular infill can replace the exposed switch without changing the
incoming widths, complete profile planes, source nodes or logical junctions.

The geometric primitive is in
[`source_curb_rounding.py`](../way/street_network/python_backend/source_curb_rounding.py).
The source-cycle PLAN path now uses it through the bounded discovery and
admission stages described below. A successful local repair is not a complete
city-import qualification.

## Original parameters and orientation

`FiniteSourceBank` retains an actual native `Path`, its declared bank side and
one or more ordered, disjoint finite station intervals. The Path must be the
context's actual original object. Disjoint ranges never become one hull.
The caller must establish ownership of those ranges and identify the actual
exposed bank switch from the canonical boundary. Nearby XY positions alone
are insufficient. The
[prepared source-coordinate ledger](STREETS_SOURCE_CURB_PARAMETERS.md) provides
original Path, bank and station ancestry without a nearest-source search.

`fit_banks(context, arriving, leaving, radius, mode, crossing_guess)` orients
the first native bank as the right bank of its Approach and the second as the
left bank. This follows the existing counterclockwise physical boundary.
These Approach directions parameterize original curves; they do not declare
that the remote original Ends form a new intersection. No extended Route can
acquire another Section implicitly.

The exposed boundary orientation matters. Rounding the same physical banks in
the opposite order can add a completely hidden circle while leaving the actual
reflex switch intact. A candidate must therefore prove full arc exposure after
the fresh physical union; successful tangent contact alone is insufficient.

## Existing construction and mandatory gates

The helper uses `CornerRadiusPolicy` and the existing circular contact solver.
STYLE radii retain the existing finite fitting attempts and physical floor;
STRICT radii do not shrink. A candidate is accepted only if both its tangent
contacts and original-bank intersection occupy the same declared interval on
each original Path. Negative or extrapolated contacts are refused.

`patch()` closes the infill along those actual source banks to their original
intersection. Its closing sides are internal preparation edges, not new curbs
or invented centreline bridges. When adjacent source polygons share a bank,
the caller must supply the same native sampling partitions to both footprints
and patch closures. Independently sampled versions can otherwise leave small
false holes or expose an internal edge.

The caller must still verify all of the following:

- Original logical members, claims, complete incoming profiles and mouth planes.
- Full exposed circular intervals, actual endpoint tangency and equal or better
  approximation to the intended analytical curves.
- Shared source sampling, finite source support and the existing continuous
  complete-width proof, including original holes and all owned prefixes.
- Source-labelled curb values and gradients, the original height envelope,
  physical and stored mesh quality, slopes, drainage and incoming joins.
- UVs, self-overlap, physical height queries and final normal-task publication.

No successful local constructor result bypasses those normal task gates.

## Source-cycle discovery and PLAN admission

[`source_bank_reconciliation.py`](../way/street_network/python_backend/source_bank_reconciliation.py)
scans adjacent owners on the existing physical boundary. A bank witness must
come from a validated original `CurbCoordinates` ledger or an exact finite
analytical offset map. It retains the actual native Path, authored bank side,
station, disjoint owned intervals and original junction incidences.
Coordinates alone do not establish any of those identities.

The first supported families are an exposed bank-to-bank switch, an original
circular corner clipped before its exit bank, and the reversed corner-entry
case. The corner families require the surviving original tangent anchor.
Tangents come from the declared analytical circle and actual native bank;
polygon chord angles do not trigger a repair. Polynomial corner envelopes
have a different ancestry contract and are not silently treated as circles.

A radius comes from an unambiguous original junction sector. A walk across
several original sectors is eligible only when every intervening incidence is
fully consumed by this same physical owner. Compatible STRICT values remain
strict. Several STYLE sectors use their smallest preferred maximum. A nearby
bank, an unowned interval, an ambiguous original walk or conflicting STRICT
values cannot authorize another candidate.

[`source_bank_preparation.py`](../way/street_network/python_backend/source_bank_preparation.py)
inserts the fitted native contact and crossing stations into the original
shared sampling partitions. It rebuilds affected original strips and ports,
retaining every old source vertex, before adding the local infill. Both raw
and canonical boundaries receive a fresh complete-width proof. The complete
arc must remain exposed, with no interior interval gaps. Original holes,
logical members, Section claims, complete profile planes and mouth positions
must remain unchanged. Existing future-strip checks run without permission to
advance any mouth.

The generated footprint carries a typed `PreparedBankFillet`, derived from
`GeneratedJunctionArea`. Its local computational Corner supplies curve and
frame evaluation only. It is never appended to `physical_layout.corners`,
which retains the original logical Corners. Body-contact bounds skip this
formal generated-area record; ordinary foreign-area influence checks still
include its actual physical area. Source-height conditions use the explicit
native-bank profile witnesses rather than interpreting computational Ends as
new logical intersections.

After ordinary owners are finalized,
[`source_bank_neighbors.py`](../way/street_network/python_backend/source_bank_neighbors.py)
checks every generated footprint against the eventual prepared owner
inventory. A raw Deferred layout is not a substitute for its finished owner.
The guard rejects stolen neighbor mouths, positive added-area overlap and a
missing prepared owner. This is a final PLAN gate, separate from each local
constructor's checks. `SOURCE_BANK_POLICY` and the local bank records enter
PLAN and surface signatures; unrelated guard inventory counters do not become
local mesh dependencies.

The source-cycle end hook runs bank reconciliation first and retained-circle
reconciliation second. The latter is documented by its typed implementation in
`retained_arc_fillets.py`. Both retain the original logical Corner inventory.

## Verification and recorded real case

The maintained
[`test_python_street_source_bank_rounding.py`](../tests/test_python_street_source_bank_rounding.py)
tests bank orientation under reversal and rigid motion, tangent contact,
unchanged authored radius, identity validation, invalid finite inputs and
disjoint ranges that must not bridge an unowned middle, and preferred-radius
fitting that never shrinks a STRICT radius.
The focused primitive tests pass. The additional maintained
[`test_python_street_bank_reconciliation.py`](../tests/test_python_street_bank_reconciliation.py)
checks discovery under reversal and rigid motion, exact Path identity,
disjoint intervals, strict radii, missing ancestry, typed body-contact bounds
and final prepared-neighbor ownership. In particular, a raw distant mouth
cannot hide a nearer finalized neighboring owner.

The captured Beijing 1132/1156 owner requires two exposed original-bank infills.
Both have the existing 1 m STYLE radius, inside the unchanged original owned
intervals. Regenerating the original shared source grids retains all old source
vertices and complete profile planes. The general `fit_banks` API reproduces
the independently qualified circles exactly: centres, tangent stations,
radii, start angles and sweeps.

The actual normal Product and chart constructors then pass with 2,600 pavement
quads. An independent audit measures stored ratio/sine 7.6871/0.04795, stored
physical slope 0.172819 below the source-derived 0.176000 limit, incoming
longitudinal joins 0.013001 below 0.03, zero drainage fill and exact incoming
XY/profile heights. Full 2 mm source-width proof uses 604 continuous source
cells; UVs, self-overlap and bilinear physical queries also pass. Source height
and gradient fit has residual below 1.74e-10.

Evidence is under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114/`:
`original_bank1132_product01_report.json`,
`original_bank1132_product_independent02_report.json` and
`source_bank1132_general_api01_report.json`. Those initial reports use explicit
source-boundary bindings. The later `maintained_bank_cases_report.json`
records the maintained selector and normal Product path for the same actual
owner, including its captured prepared-neighbor inventory. It produces 2,600
quads with the original Corners unchanged. The Source125 captured failure
context cannot qualify the final neighbor gate because it contains Deferred
neighbors rather than a completed PLAN. Its local geometry and normal Product
checks are separate evidence, not an import qualification.

## Contacts through a consumed connector with another profile

Genuine Beijing owner 185 exposes another distinction: its original 98→97
corner is clipped by a wider fully owned bank. A uniform contact Route cannot
be extended through a changed profile or width. The observed physical
boundary continues over fully owned Section 87 to the actual incoming Section 90
bank at the other real end of fully consumed Section 97. A native 98/90
bank-intersection closure would require a negative station on 90 and is
correctly refused.

[`source_owned_connector_fillets.py`](../way/street_network/python_backend/source_owned_connector_fillets.py)
instead retains an explicit original graph walk and uses the existing
`SourceContactWalk`/displaced-patch construction. The logical radius sector is
98/97 at its actual original site; the actual contact profiles are 98/90, with
their own widths and half-crown settings. The 97 centreline supplies only the
internal closure. No uniform-width97→90 Route or invented bank extension is
created. Discovery stops at an exposed mouth, unsupported curve, unrelated
bank or ambiguous original walk. The original sector and the through-site
sector both contribute their STYLE/STRICT restrictions.

The preferred 2 m circle needs more incoming lead than the captured owner's
existing cuts. It is therefore a typed **reservation**, not a mutation of a
finished mouth. A bounded source-cycle retry inserts the patch and actual
native contact stations before the normal port search. That search recomputes
the complete profile leads, extent and future-strip guards. Logical Corners
remain unchanged. Only after this preparation does `bind_contacts()` seal the
two `NativeCurbContact` records to the final, disjoint owned ranges. An unsealed
reservation cannot authorize a source-height condition. The eventual prepared
neighbor guard still runs on the final PLAN.

The initial generic reservation evidence is
`genuine185_successor_retry_report.json`: exact preferred 2 m radius, full
exposed arc, hidden actual 97 centreline closure, one original hole, 3,280
continuous source cells within 2 mm, and a passed captured prepared-owner guard.
The normal construction adds its existing 0.25 m profile lead, giving 98/90 cuts
of 15.44618 m and 8.50090 m. It removes both the clipped original corner switch
and the later 87/90 switch, so the separate 4 m infill is no longer proposed.
This report is an
artifact construction experiment. Fresh maintained source-cycle, source-height
Product and complete native-import results must be reported separately.

The subsequent `genuine185_evaluation_batches_report.json` records fresh
maintained source-cycle construction with original evaluation-batch provenance,
without a source patch or manually selected contacts. Source conditions then
pass in `genuine185_conditions147_report.json`: incoming profiles and planes
are exact, height jumps are zero, and the boundary-jet fit residual is
4.02e-11. `genuine185_normal_product_report.json` uses that immutable numerical
source snapshot and calls the normal `CompoundProduct.build(...,
source_heights=True)` and chart constructors without supplied conditions, cuts,
mesh or chart. It produces 28,625 pavement quads and 108 grass quads.

The separate `genuine185_normal_product_independent_report.json` verifies
stored edge ratio/sine 7.84808/0.050135, stored physical slope 0.128992 below
the source-derived 0.130370 limit, incoming longitudinal slope 0.011716 below
0.03, and zero drainage fill or isolated minima. All 204 incoming XY/height
profile samples are exact. Source-width, UV, self/grass-overlap and physical
query checks pass, including unchanged queries after geometry release. The
normal Product takes about 330 s and peaks at 321 MiB resident in this local
run; the independent audit takes about 18 s. These close the local constructor
gate, with current geometry and no supplied solution. They do not constitute
a fresh native coordinator Task or complete Beijing import qualification.

The portable
[`test_python_street_owned_connector_fillets.py`](../tests/test_python_street_owned_connector_fillets.py)
covers native/contact versus logical-sector identity, rigid motion and reversed
source directions, consumed/scope boundaries, through-site STRICT radii,
duplicate reservations and refusal to publish against old short mouths.
Fresh complete native imports and normal publication remain the authoritative
qualification; manually selected historical contacts are not that evidence.

## A certified core followed by an original circle and bank

The Beijing 691 investigation identified a separate finite geometry
family. Fresh original Task outputs from `runtime26_target_beijing122` reproduce
a reflex switch from a certified inner-width quintic to an actual original
circle. The existing 4 m minimum-curvature certificate and both original curves
remain unchanged. Fitting an auxiliary circle to those two finite exposed curves
fails: its circle-side contact lies beyond the true original arc endpoint.
Neither a negative native station nor extrapolating that circle is permitted.

The observed boundary has an explicit next member: that original circle's real
Section 350 bank, inside its fully owned native interval. An artifact prototype
fits a 2 m auxiliary circle between the finite certified core and that actual
native bank. Its internal closure follows the existing physical boundary chain,
including the original circle. It inserts the two true contact parameters into
the existing source grids, retaining every original source vertex. It creates
no Approach, original End, centreline, logical Corner, changed source profile or
moved mouth. The native contact is at 11.932124756 m on side -1, within the real
owned interval [0, 19.438628022] m.

`genuine691_chain168_source_conditions_report.json` records the current bounded
prototype: the full new 2 m arc and original hole remain exposed, the maintained
retained-circle construction also passes, all complete profiles and planes are
exact, and 276 continuous source cells satisfy the 2 mm width proof. Full future
strips pass without advancing a mouth. The actual core endpoint has height
0.08 m and zero gradient; the native contact has its actual zero edge height and
source gradient. The generated arc is unanchored between those real source
conditions. Boundary-field residual is 1.48e-11, and maximum slope 0.081884 is
below the source-derived 0.1 limit.

This prototype also demonstrates why a generated curve must retain its declared
source sampling budget and compulsory source knots. Its original polygon uses
a 0.25 mm chord budget; a later 5 mm sampler must not coarsen that reference
before another union. The earlier `genuine691_chain_topology167_report.json`
correctly rejects that prototype mistake. Run 168 preserves the original knots
and budget; no topology tolerance or guard is relaxed.

The maintained implementation is
[`source_parent_chain_fillets.py`](../way/street_network/python_backend/source_parent_chain_fillets.py)
with source refinement and admission in
[`source_parent_chain_preparation.py`](../way/street_network/python_backend/source_parent_chain_preparation.py).
Discovery scans only consecutive physical source owners. It requires a real
`DimensionalCompositeCurb` certified quintic, a producer-declared original-circle
interval, and that circle's actual next native bank on a claimed source strip.
Both parent Corners must retain their actual original junction incidences and
STYLE/STRICT authority. A nearby bank or an undeclared circle cannot substitute
for these identities. A STRICT constraint cannot select a smaller STYLE radius.

For core point `P(u)` and native-bank point `B(v)`, the contact solve equates the
two right-normal offset loci:

```text
Cp(u) = P(u) + r right(P'(u) / |P'(u)|)
Cb(v) = B(v) + r right(B'(v) / |B'(v)|)
Cp(u) - Cb(v) = 0
```

Parameters stay strictly inside the actual exposed finite parent intervals.
The initial bounded grid has 33 by 33 pairs; at most 32 Newton iterations use
finite bounds and a five-step decreasing line search. Each radius candidate
also needs continuous regularity certificates: `r κcore < 1` and
`(|native offset| + r) κpath < 1`. A native Exact interval crossing a miter is
excluded from this branch. There are at most 16 discovered chains and 17
radius attempts per chain. No extrapolation, source reassignment or hidden
global geometric search is available to the solver.
The native bank direction comes from its original Smooth frame or the exact
affine-miter derivative on one Exact segment. Differencing translated bank
positions is unnecessary and can produce a numerical tangent floor. Bank
points, residual limits and original controls stay unchanged.

The new `PreparedParentChainFillet` is a `GeneratedJunctionArea`, not a logical
Corner. It records its actual refined core contact, sealed native contact,
original circle and complete existing physical closure ancestry. The global
original Corners and Approaches remain identical. Every raw source point is
retained; only true source contact parameters are added. Whole-arc exposure and
hidden internal closure are checked both before and after the subsequent
retained-circle construction. Continuous source width, original hole identity,
complete profile planes, unchanged future strips and eventual prepared owners
remain independent mandatory guards.

`genuine691_maintained_chain170_report.json` records the maintained generic
selector and preparation, with no manually selected pair or injected geometry.
It passes in about 16 s, reproduces the original bounded prototype's exact
2 m contact geometry, and retains the other two original-circle fillets. The
portable
[`test_python_street_parent_chain_fillets.py`](../tests/test_python_street_parent_chain_fillets.py)
covers finite contacts under rigid motion and reversed source storage,
regularity and Exact-miter refusal, compulsory sampling knots and chord budgets,
and rejection of clipped arcs or exposed internal closures. Twenty-one joint
parent-chain, source-circle, bank and connector tests pass in
`parent_chain_tests172_report.json`.

The formal source-height dispatch uses
[`source_parent_chain_conditions.py`](../way/street_network/python_backend/source_parent_chain_conditions.py).
It validates the exact refined core object/parameter, producer circle authority,
actual native Path/side/owned ranges/contact knot, both original radius policies,
and the prepared arc's original knots. The core condition comes from the
existing source rule; the native condition comes from its own original profile
and physical covector. A generated arc has zero along-arc endpoint height
derivatives, so incompatible tangential jets or an outward-facing source jet
are rejected. Their norms cannot be substituted by clamping the vector.
The exact existing radius-candidate sequence is checked again. The new arc
remains unanchored between the real source conditions.

`parent_chain_conditions180_report.json` checks the actual maintained
descriptor and rejects seven changed-contact/authority cases. The five focused
parent-chain tests, including endpoint-jet refusal, pass under both the bundled
Python and Python 3.10 in `parent_chain_tests181_py13_report.json` and
`parent_chain_tests181_py10_report.json` respectively. Earlier sampling and
integration reports retain their original code versions.

The normal maintained boundary dispatch and `BoundaryJetField` construction
also pass on freshly regenerated original inputs in
`genuine691_normal_conditions182_report.json`: height jumps and incoming
profile/plane errors are zero, field residual is 1.48e-11, and boundary slope is
0.081884 below 0.1. The earlier normal Product attempt 176 then fails its planar
quad grading; the source conditions are not used to hide that failure.

The [finite owned-Bundle normal chart](STREETS_OWNED_BUNDLE_NORMAL_CHARTS.md)
documents the subsequent bounded computational alternative. Fresh normal
Product 189 and independent audit 190 then pass with 5,941 pavement and 142
grass quads, exact incoming profiles, unchanged controls and source widths,
stored slope 0.086350 below 0.1, and no drainage fill. The independent original
body scope is complete and contains zero eligible hard nodes in this junction
support. No saved geometry, chart, mesh, cuts or conditions are supplied to
normal Product construction. These local gates do not accept the eventual
prepared-neighbor inventory, native coordinator Task or complete Beijing import.
