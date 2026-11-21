# Original source lead and core rules

`source_offset_rules.py` supplies finite, declared jet pieces for original
intersection boundaries. It changes no physical curve, chord, corner radius,
mouth plane or profile. `SourceCurbConditions.derived_rule()` uses the helper
where an original boundary lead crosses a real lead/core contact or needs its
actual parameter-map seams for a residual-qualified field fit.

## Original offset leads

`BoundaryCurve` defines each lead by its original approach, offset, source
station parameters and cumulative offset-curve distances. Its own evaluator
maps curve distance to source station with that table. Those actual curve-map
seams must be retained when approximating a gradient along the lead; they are
different from a native Path's adaptive integration lookup table.

`SourceOffsetRule` declares:

- Those existing curve parameter-map seams.
- Original native Path polynomial/Exact control boundaries, mapped to metric
  source stations with `path.station(path.knots)`.
- Real Route piece boundaries, with the original forward/reverse ancestry.

For an ordinary Smooth source, the original native unit normal gives its
transverse profile covector directly. Exact sources retain their interpolated
miter normals and use the analytic one-sided source Jacobian on the selected
original segment. Complete incoming `ProfileEdge` and `Approach` planes remain
unchanged. A derived approach with its own physical field continues to supply
its own physical gradient.

Profile lookup follows the actual native Path and its orientation. It introduces
no nearest-source assignment or unowned station extension. Prepared dimensional
`OffsetCurb` rules retain their existing `_OffsetRule` protocol and results.

### Exact zero seams and original auxiliary leads

The lead's declared construction-table and native metric knots retain their
exact source stations for one-sided jet evaluation. If interpolation of a true
zero knot returns a roundoff-sized negative number, that result cannot choose a
different construction branch. Published exact knot parameters identify the
station; no arbitrary near-zero clamp is used. The primal station mapping and
all physical curve evaluations remain unchanged.

On a wholly negative original auxiliary interval, `Approach.frame` defines the
actual affine construction `p = p0 + s T0 + d N0`. Its constant physical
covector g is obtained from

```text
dot(g, T0) = 0
dot(g, (width / 2) N0) = profile.derivative(q)
```

This uses the actual original frame basis, including source reversal and rigid
transforms. It replaces a noisy finite-difference stencil on that exact affine
extension. It does not assign negative stations to a native body. A genuinely
undeclared interval crossing station zero remains an error. Positive native
pieces keep their analytic native Jacobian and exact one-sided endpoint jets.

Portable tests are in `test_python_street_source_offset_auxiliary.py`. The
actual Bergamo 133-Section component's original site 3107 previously failed two
lead pieces around zero with polynomial residuals above `1e-8`; the same saved
actual group now rebuilds normal conditions and the complete boundary field
without those failures. No physical curve, boundary sample, fit threshold,
polynomial degree or prediction budget is changed.

## Lead/core and composite spans

A physical chord can span a declared curve contact without containing a mesh
vertex at that contact. `SourceCurvePieces` retains the original physical chord
and a complete finite parameter ledger for its source pieces. Each lead keeps
its original source rule; each core keeps the existing bounded quintic source
rule. The class exposes internal `polynomial_breaks()` and one-sided
`sample_piece()` evaluations to `BoundaryJetField`.

Only polynomial prediction is partitioned. Physical boundary coordinates,
closest-edge projections, weights and topology remain unchanged. The existing
20,000 prediction-branch, degree-20 and `1e-8` fit limits remain mandatory.

## Paired incoming field chord domains

A paired curved entrance has four actual Cartesian curb polylines. Their
ray-intersection derivatives change at each declared chord endpoint. Fitting a
single polynomial across such a derivative change can fail even on a short
physical boundary edge. Increasing the polynomial degree or accepting its
residual would hide that change rather than describe the supplied field.

`paired_source_jet_domains.py` admits the exact `PairedAngularEntrance` and
`AngularCrown` producer. It verifies both canonical native Paths, real directed
run endpoints and adjacent visits, original half-profile orientations and the
shared rise/datum, coordinate-frame identity and finite chord lattices.
Ordinary native Paths retain their existing branch. Unsupported derived fields
do not gain this producer certificate.

The entrance's outward station `s` and chart station `ξ` satisfy

$$
\xi=\xi_0+\sigma s,\qquad \sigma\in\{-1,+1\}.
$$

Every actual chord endpoint `ξₖ` therefore supplies the exact declared station

$$
s_k=\sigma(\xi_k-\xi_0).
$$

The existing monotone source-offset parameter map carries these stations into
prediction breaks on the physical edge. Existing source-map seams remain
declared as well. Only the positive finite incoming domain uses this branch;
negative auxiliary construction keeps its existing separate treatment.

`RayCurb.jet(interval=...)` selects the actual chord containing the declared
interval. At its endpoints it returns that chord's one-sided derivative.
Intervals crossing another undeclared chord knot or queries outside the
declared interval fail. Representability margins concern only comparisons at
the exact shared knot; they do not accept a geometric boundary displacement.
`AngularCrown.gradient(chord_interval=...)` evaluates its existing analytic
height gradient with those four consistent one-sided chord jets. Default
gradients and the existing linear-field subclass interface stay unchanged.

No physical boundary vertex, profile, height value, field projection, blending
weight or mesh topology is changed by this prediction partition. The existing
`1e-8` residual, degree-20 and 20,000-branch limits still apply. The real authored
curved-Bundle failure drops from a `5.98e-6` residual to `9.91e-15` with 183
prediction branches on its unchanged 100 physical boundary edges.

`SOURCE_PAIRED_JET_POLICY` identifies this typed producer;
`SOURCE_HEIGHT_POLICY = 5` invalidates the earlier surface implementation.
Portable tests are in `test_python_street_paired_source_jets.py`. They cover the
normal paired source conditions, one-sided endpoint derivatives, undeclared
crossings, borrowed Paths/frames, altered profiles/end identities, finite
station bounds and preservation of ordinary/linear interfaces. A scoped
numerical candidate passes these tests with the source-offset, auxiliary,
physical-binding and boundary-piece suites (27 methods). Complete current-source
and native qualification remains in the city acceptance report.

## Whole original curve authority

The generic clipped-corner reconstruction must not replace a complete original
junction's known lead/core construction. Such authority is admitted only when
`whole_original_curve()` independently verifies:

- One actual original junction, no consumed claim, no composed physical layout,
  one original physical footprint and one outer loop.
- The exact physical layout object supplied by the runtime's
  `physical_intersection` dependency.
- The original distinct directed End inventory, actual Section/Path identities,
  original oriented profiles, widths, offsets and authored radius inventory.
- Complete contiguous exposure of the same original curve over parameter
  `[0,1]`, retaining every source span.

Clipped curves, unions, wrong scopes and stale layouts retain the existing
reconciliation or refusal. No `original_junction` metadata is fabricated. The
existing core rule is retained only after this proof; it is not replaced by a
new taper or solver.

## Validation and remaining gates

Focused tests are in `test_python_street_source_offset_rules.py` and the existing
source-curb, returning-road and continuation suites. They cover actual runtime
construction, exact lead-map knots, source Jacobians, contact-crossing chords,
and rejection of clipped, union, stale-layout and altered-radius authority.

The near-straight width continuation, Exact returning road and side-junction
runtime construction pass with this policy. The separate fixed grass quad-ear
failure has a bounded [direct-quad recovery](STREETS_FIXED_GRASS_QUAD_EARS.md).
Its physical boundary and ten-percent slope limit remain unchanged. Full city
and combined current-source qualification is tracked separately.
