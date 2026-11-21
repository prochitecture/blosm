# Finite owned-Bundle normal charts

This is a late computational chart family for a source compound with one or
more real holes. It changes no centreline, source profile, rounded corner,
physical ownership or incoming profile plane. Implementation:

- [`owned_bundle_normal_atlas.py`](../way/street_network/python_backend/owned_bundle_normal_atlas.py)
- The late physical-hole recovery sequence in
  [`partition.py`](../way/street_network/python_backend/partition.py).
- Portable contracts in
  [`test_python_street_owned_bundle_normal_atlas.py`](../tests/test_python_street_owned_bundle_normal_atlas.py).

## Why this scope is separate

The existing [paired-prefix chart](STREETS_PAIRED_PREFIX_CHARTS.md) requires two
typed terminal `PartialBundleClaim` records with nonempty completely consumed
member prefixes. A physical compound can also own a finite incoming portion of
one member without consuming an entire Section on that side. That does not
authorize constructing a missing partial claim.

`OwnedBundleNormalScope` instead reads the actual prepared native ranges. Each
range is certified by a fresh continuous source-width proof, including explicit
hole exclusion, with the existing maximum 2 mm error. The range mapping must
equal `source_bank_reconciliation.owned_ranges(context, group)` exactly.
Separate intervals remain separate: their unowned gap is never replaced by a
hull. Path objects and their source Section objects retain their actual
context identities.

There must be exactly one original complete two-member Bundle intersecting the
compound's claims, and both actual member roles must occur in the proven owned
ranges. All other original claims stay in the same physical compound. This
scope does not add a Bundle, change membership, invent a connection or claim
another part of the network.

## Computational cuts

After all earlier physical-hole, paired-prefix and original-corridor choices
fail, the solver tries each original role, at most twice in total. It reuses
`PhysicalHoleCuts` and `SourceCycleEnergyChart`.

For a selected native centreline `C(s)` and its actual normal `N(s)`, the
computational cross-section is

```text
L(s, t) = C(s) + t N(s),    s in one certified owned native interval.
```

Intersections with the unchanged physical rims define the actual pavement
intervals. Both sides of each original member remain physical inputs; selecting
a role chooses computational axes only. Unrelated source normals retain their
existing rule for selecting the interval containing that centreline.

The existing finite budgets remain:

- Two role candidates and at most 16 original compound members.
- At most 64 bands per original Section, shared across its disjoint intervals.
- The existing 128 computational-cut and eight-hole limits.
- Cadence multiplier 2 and an eight-source-width cap on unpaired cut spans.
- The caller's unchanged mesh-cell and quality-correction budgets.

Every cut ends on an unchanged physical rim. Complete incoming profiles remain
prescribed data. The ordinary quad algorithms and bounded energy correction
operate inside those fixed boundaries; no triangle-first construction is
introduced.

## Admission and height semantics

Both pavement and grass must pass the existing physical and stored-coordinate
quad-quality checks and total cell budget. Continuous source-width proof is
recomputed after chart construction. Original members, claims and native ranges
must be unchanged. Any typed parent-chain fillet also repeats its whole-arc and
hidden-closure check after the chart's sampling.

When the Surface caller supplies `SourceHeightChartAdmission`, the same normal
admission callback checks necessary entrance feasibility against the original
profile/body data before accepting a role. It assigns no final height solution
and relaxes no solver guard. Final source heights, complete profile planes,
longitudinal and surface slopes, drainage, stored buffers, UVs, overlaps and
physical queries still require normal Product construction and independent
validation.

The cache policy is `SOURCE_OWNED_BUNDLE_NORMAL_POLICY`.
`signature(context, group)` includes the exact range ledger, actual Bundle,
ordered member traversals, and relevant original Path/Profile digests. An
inapplicable optional recovery returns no signature; a valid earlier chart
cannot be rejected merely because this unused late scope is inadmissible.
Strict construction still rejects an invalid scope. PLAN, Surface and seam
epochs include the global policy token independently.

## Evidence and remaining gate

Development evidence is under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114/`.

The genuine Beijing 691 owner first passes fresh source preparation and normal
source conditions. The earlier normal Product attempt 176 then rejects the
available chart after about 218 s: 12,673 pavement quads still contain cells
with edge ratio 46.89 and corner sine 0.00372. Neither existing late atlas is
eligible, because the physical range pattern lacks two terminal prefix claims.

The scoped prototype in `genuine691_owned_native_chart185_report.json` passes
with 5,941 pavement and 142 grass quads, using the first real role. Chart
construction takes about 8 s. Stored pavement ratio/sine are 7.87278/0.134917;
grass values are 7.72024/0.09685. The 276-cell continuous source-width proof is
unchanged. `genuine691_owned_chart_admission186_report.json` also passes the
normal necessary entrance/body admission. These are chart experiments, not a
normal Product or native import qualification.

The maintained scope, optional-signature refusals, original role selection and
nearby source-family regressions pass 18 tests in
`owned_bundle_normal_tests188_report.json`. The later optional-cache test also
covers missing legacy interval metadata without rejecting an unused recovery.
Python 3.10 run 191 exposed synthetic fixtures using a literal 10 m bound where
that runtime's native Path is a few ulps shorter; the strict finite guard was
retained. The new scope fixture now uses the actual Path length.

`genuine691_normal_product189_report.json` regenerates preparation from genuine
original raw Tasks, then calls the normal Product/chart/source-height path.
No geometry, mesh, chart, cuts, conditions or field are supplied. It passes with
5,941 pavement and 142 grass quads, in about 190 s including fresh preparation
and the earlier chart families that must fail before this late recovery. Peak
resident memory is about 255 MiB. Original inputs and the immutable numerical
source snapshot are unchanged.

`genuine691_normal_product190_independent_report.json` separately checks stored
ratio/sine 7.87278/0.134917, stored physical slope 0.086350 below 0.1, entrance
longitudinal slope 0.014227 below 0.03, and zero drainage fill or interior
minima. All 119 original incoming XY/height samples are exact. Source-width,
finite UV, pavement/grass overlap, physical queries and geometry-release query
invariance pass. Independent source-body reconstruction reports a complete
supported scope with no unsupported fields. This owner has **zero hard body
nodes** after excluding its real junction support; no nonzero hard-target
accuracy is claimed for this case. Query residual is below 1.72e-15.

These results close this genuine owner's local unsupplied constructor and
independent mesh gates. The complete prepared-neighbor inventory, native
coordinator Task and whole-city import remain separate pending runtime gates.
