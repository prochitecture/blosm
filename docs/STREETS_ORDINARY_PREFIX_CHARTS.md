# Ordinary original-source prefix charts

`way/street_network/python_backend/ordinary_source_prefix_atlas.py` supplies a
bounded, failed-only computational chart for an ordinary original-root source
union. It is independent of the broader consumed-source cycle atlas and does
not expand `source_cycle_atlas.eligible()`.

The protected `CompoundQuadChart._ordinary_source_prefix()` hook runs only after
the existing assembled-grading, corner seed and ordinary sampling/recovery
paths have failed. Successful existing charts keep their original behavior.
The retry subclass disables that hook, so a rejected candidate cannot recurse.

## Admission and source ownership

`OrdinarySourcePrefixScope.build()` creates a typed scope from the current
original network, the existing original-root-union admission and an independently
rerun `certify_layout()` proof. Its requirements are:

- Exactly one original `Intersection`, with the complete, distinct original
  directed `End` inventory and no consumed Section claim.
- The existing `FINITE_ORIGINAL_SINGLE_JUNCTION_SOURCE_UNION` witness, retaining
  original Path/Profile incidence, complete source-grid proof and strict radius
  inventory.
- Actual original Section and Path object identities for every Route piece;
  full, centered incoming profiles equal to the original oriented profiles.
- No lateral Bundle membership, partial Bundle claim, absorbed site, prepared
  mouth extension or finite-scope port.
- Finite, disjoint owned station intervals from a fresh complete-width proof,
  agreeing with the current stored proof. Disjoint prefixes of a returning
  Section remain separate; their unowned middle is never sampled as a hull.

The scope retains exact physical boundary coordinates, corner radii, incoming
mouth stations and the original upstream objects. The computational cut policy
rejects reuse after any of those inputs change. This is a computational chart
policy, not permission to move a mouth, fit a radius, add an incidence or consume
another source.

## Finite candidates and unchanged guards

Rank original sources by their **owned** interval length:

\[
  L_i=\sum_{[a,b]\in\mathcal R_i}(b-a).
\]

Try the three longest sources individually, then their combined axes. Repeat
that same finite list with cadence two only if cadence one fails. Duplicate axis
sets are removed, giving at most eight attempts. The all-pavement-interval mode
includes every physical pavement intersection along a computational normal; it
does not grant ownership to any additional longitudinal interval.

Each candidate uses the maintained `PhysicalHoleCuts`, `SourceCycleEnergyChart`,
direct quad construction and existing bounded energy correction. Source normals
retain the actual metric eight-width span limit, at most 64 shared bands per
original Section and 128 computational cuts. The 20,000 boundary-sample and
50,000 compound-cell budgets remain unchanged. Candidate admission reruns the
complete source-width proof and checks physical and published float32 quad
quality (edge ratio at most eight; corner sine at least 0.03). No new optimizer,
surface-height policy or numerical tolerance is introduced.

### Optional original-height admission

`CompoundSurfaceBuilder(source_heights=True)` supplies an optional
`SourceHeightChartAdmission` callable. Only this new failed-only prefix family
uses it, after its normal planar checks. The default geometry-only constructor
supplies no callable and keeps the original height-independent candidate order.

The callable prepares the actual configured boundary conditions and original
body targets without solving heights. For each authoritative entrance midpoint
row \(a\cdot z\), immutable profile/body heights are retained and all remaining
heights are allowed anywhere in the original datum/crown box \([d,c]\). Its
necessary feasible interval is

\[
 I_{\min}=\sum_{i\in K}a_i z_i+
          \sum_{i\notin K}\min(a_i d,a_i c),\qquad
 I_{\max}=\sum_{i\in K}a_i z_i+
          \sum_{i\notin K}\max(a_i d,a_i c).
\]

The row's necessary slope floor is
\(\max(0,I_{\min},-I_{\max})\). Arithmetic is expanded outward before testing
the builder's configured entrance limit. If the interval cannot meet that
limit, the next existing axis/cadence candidate is tried. The callable neither
fixes a free value nor releases a mandatory value. Passing this necessary test
does not establish feasibility: the existing height, drainage, publication and
actual entrance checks remain mandatory.

`OriginalBodyReference.inverse()` retains every already-positive normal-root
certificate. When its old coarse bound fails, it may use a tighter bound over
the same complete station bracket and the same global curvature bound. For a
unit-speed original curve with \(\lVert C''\rVert\le\kappa\), bracket length
\(L\) and query point \(p\),

\[
 \sup_s\lVert p-C(s)\rVert\le
 \max(\lVert p-C(s_0)\rVert,\lVert p-C(s_1)\rVert)
       +\kappa L^2/8.
\]

The extra term bounds deviation from the endpoint chord. The same positive
normal-root Jacobian and endpoint sign tests, finite domain and physical inverse
residual checks still apply. A failed uniqueness certificate remains fatal.
This changes neither the original Path nor its curvature bound or bracket.

### Declared source pieces in the boundary field

A native source bank can cross an original Path metric knot or a real Route
seam inside one physical boundary chord. Its exact jet may need separate
polynomials on either side of that declared source boundary. Native rules can
expose `polynomial_breaks()` and `sample_piece(local_t, first, last)`; the latter
retains true one-sided endpoint jets for exact source corners.

`BoundaryJetField` compiles these internal prediction pieces with the existing
degrees 5, 10 and 20 and the unchanged `1e-8` residual requirement. Original
physical edge origins, chords, closest-edge projections and blending weights
remain unchanged. Only polynomial evaluation chooses a declared piece. Total
prediction branches must fit the existing 20,000 cap. Adaptive arc-length lookup
table entries do not define pieces. Rules without declared interior breaks keep
their original fitting and evaluation arithmetic.

The actual singleton's field has 392 physical edges and 395 prediction branches,
with maximum qualified fit residual `9.30143e-11`. This repairs the two native
knot spans that failed the single-polynomial approximation in immutable source02;
it does not relax an elevation, slope or surface-quality guard. Focused tests in
`test_python_street_boundary_jet_pieces.py` verify unchanged physical weighting,
one-sided jets, bounded declarations and rejection of unresolved fits. The
existing boundary-polynomial reference tests verify the unchanged ordinary path.

## Dependency and cache contract

`signature(context, group)` projects the policy token, original Junction and
directed Ends, exact owned station ranges, physical boundary/mouth/radius token,
and digests for every original Path/Profile used by any approach Route piece.
The original Junction digest is included. Surface Task declarations must cover
that same upstream closure. The signature does not rerun geometric proofs on
cache lookup; construction always reruns them before generating cuts.

## Validation scope

Focused tests are in
`tests/test_python_street_ordinary_source_prefix_atlas.py`. They cover directed
source reversal, owned-length ranking, disjoint/stale ranges, wrong original
graphs, original profile/Path identity, finite attempt order, recursion rejection
and exact preservation of a previously successful normal chart.

Generated qualification reports and diagnostic meshes belong under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114`, outside the
repository. A supplied-chart height experiment is not a normal-chart or Task
qualification. Full physical height, drainage, entrance, body-row, publication
and normal coordinator Task gates remain independent requirements.

### Native singleton chart qualification, 2026-11-16

The actual original-source singleton at Bergamo site 731 was constructed through
the normal `CompoundQuadChart.build()` dispatch. No mesh, computational cut
policy or boundary conditions were injected. The first candidate used original
Section 366 and cadence one, generating 3,081 pavement quads. Physical edge ratio
was 7.99762046 and minimum corner sine 0.03439495; after float32 publication these
were 7.99985416 and 0.03413164. Exact physical loops, incoming mouth stations,
corner radii, raw ownership and complete owned-width proof were preserved.

The coherent immutable source contained 1,519 files. All source hashes and the
input network remained unchanged during the replay. Evidence is in
`source351_normal_prefix_chart01_report.json` and
`frozen_ordinary_prefix_chart_source01/source_manifest.json` under the task's
artifact directory. The normal chart took 112.70 seconds in this concurrent,
resource-guarded diagnostic; this is not a controlled performance measurement.

The subsequent normal `CompoundProduct.build(source_heights=True)` still failed
the independent physical height gate. An additional supplied-chart diagnostic
proved why merely rerunning a height solver cannot fix this first chart: one
incoming midpoint gradient row has four immutable vertices (two actual profile
vertices and two exact original body targets), giving longitudinal slope
0.03624643, above the configured 0.03 limit. Global fixed-gradient feasibility
was satisfactory, so that result concerns entrance alignment on this particular
chart. Mandatory body targets and incoming profiles must be preserved while a
different bounded chart or conforming refinement is investigated. This is a
normal **chart** qualification, not a normal Product, Task or full-city import
qualification.
