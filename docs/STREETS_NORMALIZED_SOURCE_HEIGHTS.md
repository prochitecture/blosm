# Height jets on a normalized source rim

## Physical and computational references

`source_rim_normalization.ResampledSourceCurb` supplies a derived polyline
carrier within its existing one-micrometre additional projection allowance.
The original analytical curves, original polygons, source-span ledger and full
incoming profiles remain the physical references. The carrier is not an
assertion that an original circle or offset curve equals a straight chord.

`source_curb_conditions._MappedRule` transfers each original source jet through
the carrier's explicit parameter map. Its internal polynomial pieces follow
actual construction seams. These seams are distinct from adaptive arc-length
integration tables; no unrelated quadrature knots are introduced.

## Maps and Jacobian

Let s be carrier polyline distance, and let S(s) be original sampled-rim
distance. `curve.stations` and `old_at` define the piecewise affine map S.
`old_stations` identifies the original source-span rules. Each source-span rule
may also declare exact native Path metric knots or Route seams.

Within one carrier and original-source piece, let T be the carrier unit
tangent, t the original sampled-rim unit tangent, and R(v) the clockwise
right-normal operation. With original physical gradient g, the transferred
gradient is

```text
lambda = dS / ds
g_carrier = dot(g, t) lambda T − dot(g, −R(t)) R(T)
```

Height is sampled from the actual original rule. This preserves the along-rim
derivative under the known distance map and the original inward-normal
derivative in the physical carrier frame. A carrier interval spanning a table
seam must use each piece's own tangent and lambda.

## Declared prediction pieces

The rule declares only:

- actual carrier polyline knots;
- original source-span knots mapped through the explicit distance ledger;
- any nested native rule's exact declared piece knots.

`sample_piece(t, first, last)` evaluates a consecutive declared map interval.
When a clipped carrier interval covers only part of a native rule branch, it
resamples inside that native branch while preserving the native branch's own
declared bounds. Thus an Exact native control receives its true one-sided jet;
two sides are not averaged.

Published normalized table parameters identify their exact authored station
knots. The implementation recognizes these exact parameter values directly.
Multiplying such a fraction back by total length must not select the adjacent
source rule because of roundoff. At a physical edge endpoint, the source jet
comes from the branch inside that edge. Existing interior arithmetic and
already-correct endpoint arithmetic remain unchanged.

`BoundaryJetField` keeps its original physical edge projection, origins,
closest-edge weights, degrees 5/10/20, `1e-8` fit gate and total 20,000 prediction
branch cap. Pieces change polynomial evaluation inside a physical edge; they
do not create physical vertices, alter sampled geometry or increase a budget.
Slope, entrance, drainage, body-target and stored-float publication guards
still run after field construction.

## Regressions and evidence

Maintained tests are in `tests/test_python_street_mapped_source_rules.py`.
They cover real source-map seams, a clipped one-sided native branch,
published-knot endpoint selection, the map Jacobian and byte-identical sample
arithmetic for an unaffected single old piece. Existing boundary polynomial
and explicit-piece tests run with them.

The actual Bergamo handoff investigation found that all 668 original sampled
source conditions compiled before normalization, while the accepted normal
chart's mapped carrier conditions produced 476 polynomial-fit failures. The
problem was both undeclared source-map seams and incorrect endpoint-side
selection; it was not evidence that the tiny original corner needed a looser
fit tolerance.

The qualified algorithm experiment compiled 656 physical edges into 678
prediction branches, with 16 split edges and a maximum sampled fit error of
`6.241048997046583e-11`. Projection and physical weights were unchanged.
The normal Product subsequently passed independent body, profile, slope,
drainage, mesh-quality and stored-float checks. See
[Exact original-circle bank handoffs](STREETS_ORIGINAL_BANK_HANDOFFS.md) for the
scope and limits of that experiment.

Regenerate the synthetic tests from maintained source. Actual-network
qualification must use fresh original input and a coherent runtime candidate;
development pickles under `C:/dev/workspace/artifacts` are optional diagnostic
evidence and are not a recovery dependency.
