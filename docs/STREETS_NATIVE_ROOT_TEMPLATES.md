# Original source root templates

This is a failed-only construction family for one original ordinary junction.
The established circular and zero-curvature root families run first. Successful
existing constructions keep their original curves and arithmetic. This family
does not add a mesh solver, raise sampling/cell budgets, move centerline
controls, replace profiles, or override an explicitly STRICT radius.

## Why a contained source root may still need repair

A full-width proof establishes that each native street strip is contained in
the **closed** pavement domain. Its transverse root cross-section can still be
part of the exterior. That is an artificial preparation cap, rather than a
carriageway curb. Its interior profile vertices can also have positive height.
The root audit therefore independently requires every original transverse root
to be internal. Width containment is necessary but does not replace this test.

Use the actual prepared original ports once as physical operands. Rebuilding a
second approximation of the same native bank with another station grid can
create an invalid planar arrangement before a proposed curb is added. An
invalid fixed arrangement must never be treated as a rejected new candidate.
The reference validates the fixed operand union independently, retains every
other original local corner, and removes only the donor being replaced.

## Mixed native-profile reference under qualification

The full-profile reference remains unchanged. A separate diagnostic reference
admits a genuine Bundle half profile on an unrelated original approach when
**both donor banks still belong to full profiles and have their original zero
curb heights**. It authenticates every original Approach, Route, Path and
membership-resolved ProfileProvider result. The required exposed root must be
an actual zero-height half of a non-donor approach; an already internal half
cannot be invented to satisfy this condition.

The same native-curvature core, finite bell families, displacement bounds and
normal field/mesh solver are reused. Actual original port polygons replace the
synthetic proof strips only through their own typed reference and independent
fixed-operand union. All other original area operands remain authoritative.
The signature includes the real oriented profile halves and every relevant
Bundle/member Section dependency. Source-cycle replay uses that exact reference
to authenticate retained-circle contact sampling.

Local normal construction and an independent 6,636-quad mesh audit pass for the
four-arm Karl-Marx-Allee case, including 68 exact profile samples, 31 exact body
targets and zero source-curb height error. Five authored authentication methods
and three existing full-profile regression methods pass. A fresh native import
gets past this root and exposes a different source seam, so this is not yet a
qualified complete native implementation. See the current
[city acceptance report](STREETS_BEIJING_BERGAMO_FIXES.md).

An exposed chain consisting of the donor approaches' **opposite** root halves,
with no exposed non-donor root, lies outside this contract. Do not force an
internal half into `required_roots`, flatten a Bundle profile or borrow the
full-only reference to accept such a case. It needs its own actual-bank-chain
construction and complete area, profile and final-owner proofs.

## Geometric endpoint curvature

The original quintic continuation has zero second endpoint derivatives. At a
curved source bank, its zero-curvature continuation can enter the carriageway
at second order. A bell with zero endpoint position/first/second jets changes
only third and higher orders and cannot correct this local defect.

Let a native centerline have signed curvature `k`, with its unit right normal.
At signed lateral distance `d`, its regular offset bank has curvature

```text
k_bank = k / (1 + d k),       1 + d k > 0.
```

Route reversal changes the signed curvature and lateral orientation explicitly.
The first root endpoint is traversed toward the original junction; the second
is traversed outward. No nearest-bank search determines these signs. True
control knots or Route seams requiring an undeclared one-sided jet are refused.

For the existing endpoint velocity `v`, use geometric normal acceleration

```text
a = k_bank |v|² left(v / |v|).
```

The quintic Bézier controls retain `P0`, `P1`, `P4`, and `P5`, and set

```text
P2 = 2 P1 - P0 + a0 / 20,
P3 = 2 P4 - P5 + a1 / 20.
```

Thus positions and endpoint first jets stay fixed. This matches geometric
curvature (G2); it does not claim that the source bank's external parameterization
has the same complete parametric second jet. The displacement of this new base
relative to the old base has its own Bézier convex-hull bound. It is reported
separately from the bell displacement cap.
The report also supplies a conservative total bound equal to the base
convex-hull displacement bound plus the actual bell peak. The bell cap alone
is not a bound on displacement from the old zero-curvature base.

## Bounded bells and unchanged certification

The existing degree-six bell remains the first bell. Two finite alternatives
are normalized degree-nine Bernstein bases:

```text
b_i(t) = t^i (1-t)^(9-i) / [(i/9)^i (1-i/9)^(9-i)],   i in {3, 6}.
```

They peak at `1/3` and `2/3`, are nonnegative, have maximum one, and have zero
position/first/second endpoint jets. Exact polynomial degree elevation retains
the quintic base before adding `distance * b_i(t) * outward`.

The physical bell displacement is bounded by

```text
max(0.25 * donor_min_width,
    0.5 * max(0, actual_required_root_width - donor_min_width)).
```

An exact root projection supplies the required amplitude, followed by the
existing small chord-error margins. There is one compulsory candidate formed
from actual certified contact planes, inside the existing seventeen-reach
allocation, at most eight amplitudes per radius, and the established finite
STYLE radius sequence. No fitted new source plane is introduced.

The existing positive whole-chord speed projection and Bernstein curvature
certificate apply to the actual polynomial degree. Curvature subdivision stays
bounded by 128 intervals. Sampling uses the existing speed and second-derivative
chord-error bounds. Radius, full exposure, native widths, original profiles,
holes, stored-coordinate quality, slope, drainage, and entrance checks remain
mandatory.

## Contact sampling without adding mesh cuts

When a new tangent/curvature-matched curve is outside the intended native bank
but the coarse sampled bank overlaps it near the contact, refine only the
owned side `[max(0, contact - 0.5 m), contact]`. The external contact-to-mouth
lead retains its existing sampling. Existing source station samples are
included; all new points are evaluated from their exact original Path/Route
parameters and validated by `CurbCoordinates`.

The local analytical chord target is one tenth of the unchanged arrangement
weld distance. This is a refinement target, not a tolerance for accepting
penetration. The complete planar union and full-curve/root exposure tests run
again. Source allocation stays at 20,000 segments; mesh allocation stays at
50,000 cells. Refinement can fail within those limits.

These internal proof knots usually disappear inside the accepted curb patches.
They need not become final physical rim vertices. Keeping generated
approximation polygons bit-identical is not a user requirement: intended
analytical curves, actual controls/seams/contacts/profile planes, equal or
better approximation, and all final geometric proofs are authoritative.

## Independently necessary recipient planes

A width transition may already need to advance its narrow recipient's mouth
before any root construction. `StagedRecipientWitness.capture()` calls the
existing recipient validator on that actual pre-root group. It proves that the
old profile is genuinely incomplete and the normal contact-plus-reserve plane
is complete, with the original source, radius, and certified dimensional curve.

Final admission revalidates that pre-root proof and requires:

- the same original graph, canonical Paths and Profiles;
- the exact original previous plane and already certified current plane;
- the same dimensional corner geometry;
- no root-induced advance of any complete original plane;
- unchanged bounds for the other donor;
- complete final original profiles and a fresh full-width proof.

The ordinary constructor then performs its normal final mouth/STRICT-radius,
analytical sampling, and source-width checks. It creates its own clearance
metadata. Staged records appear under
`ordinary_complete_mouth_clearance.staged_dimensional_recipient_lead_advances`.
They are not treated as new root-donor permissions.

## Code, signatures, and qualification

The implementation is in:

- `source_root_curves.py`: curvature-matched quintic and bounded bells;
- `source_root_ports.py`: actual original port ledgers and contact refinement;
- `source_root_witnesses.py`: independently revalidated staged recipient planes;
- `source_root_templates.py`: typed failed-only scope, construction, and signature.

`SOURCE_ROOT_POLICY` and `construction_signature(context, group)` are exported
by `source_root_templates.py`. The signature includes local curve controls,
native curvature/contact/radius parameters, actual Junction/Path/Profile dependency
digests, and staged contact/plane evidence. Final/global neighbor counters are
excluded. Runtime integration includes this policy in preparation, layout,
projected-surface, and transition cache signatures. Complete PLAN qualification
remains a separate gate from the local constructor proofs.

Portable tests live in `tests/test_python_street_source_root_curves.py` and
`tests/test_python_street_source_root_templates.py`. They construct authored
source networks directly and require no scene, checkpoint, or artifact pickle.
Run them with the repository's numerical test environment and NumPy; Blender
registration is not required. Factory-profile rules still apply to any actual
Blender application launch.

Local preparation and Product proofs do not qualify a whole imported network.
The complete eventual PLAN owner inventory is required for final prepared-
neighbor acceptance, followed by genuine native import tests. Development
reports and scenes belong under the workspace artifacts directory and are not
runtime prerequisites.

### Original full-context qualification

Candidate 34 (`0.1.0.dev93`) passes the fresh original 1,255-Section Bergamo
dependency closure for actual Surface 1947, including complete PLAN preparation
and all 5,211 selected Tasks. No upstream products, task fingerprints, replacement
groups, meshes or builders are injected. Independent final-owner validation
reconstructs 617 eventual physical owners and 904 area operands from 665 original
junctions. The one typed generated root is included in the final guard, and its
published `source_root_prepared_neighbor_guard` matches that complete inventory.

The normal 5,359-quad Product retains 42 exact original body targets and all 51
original incoming profile samples. Physical/stored maximum slopes stay below
0.049702 and entrance slopes below 0.012249, with no interior minima. Independent
field/query and stored-quality checks pass. The same immutable source passes
310 combined numerical tests and native Karl-Marx-Allee/lifecycle regressions.

This qualifies the actual target and complete owner inventory for the root
construction. All remaining mesh Tasks, complete native Bergamo/Beijing imports
and original-city progressive behavior remain separate gates in
[the acceptance report](STREETS_BEIJING_BERGAMO_FIXES.md).
