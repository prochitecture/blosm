# Smoothing dense and short source controls

The Python path builder retains every authored Section ID and control point.
The changes below affect derived evaluation curves only, and are attempted only
after the previous construction rejects a Smooth path. Exact polylines and
previously accepted Smooth paths retain their established construction.

The relevant implementation is in
[paths.py](../way/street_network/python_backend/paths.py) and
[local_smoothing.py](../way/street_network/python_backend/local_smoothing.py).
Regression cases are in
[test_python_street_dense_controls.py](../tests/test_python_street_dense_controls.py).

## Dense controls

Several real Bergamo roads encode a rounded corner with very short consecutive
segments. Treating each intermediate point as an independent circular-style
bend can demand overlapping setbacks even though a smooth physical road fits.

After ordinary construction fails, the builder tries coalescing evaluation
controls within successively larger portions of the existing deviation budget.
The source data remain unchanged. Every replacement chord must satisfy all of
the following:

- Its finite capsule contains the complete corresponding source subpolyline.
- Source projection onto the replacement chord progresses monotonically.
- A sharp authored turn or reversal is retained.
- The neighbors required by fixed shared endpoint bends are retained.

The curve retry receives the original deviation budget minus the coalescing
budget. A second, independent check proves that every Bezier piece remains
inside the original finite-polyline corridor: a piece's entire convex hull must
fit one source-segment capsule, or both recursively subdivided halves must pass.
Uncertainty at the finite recursion limit rejects the candidate. This check is
continuous, rather than a test of a few sampled points.

The original source identity and all controls remain available as `path.section`
and `path.source_points`. Pavement metadata records the retained evaluation
indices, source control count, budgets and original failure under
`centerline_control_coalescing`.

## Finite bend reserves

Local rounds normally reserve a radius of at least 1.5 times the pavement's
lateral reach. When this layout preference cannot fit, the failed path tries
1.45, 1.4, 1.35, 1.3 and 1.25 times that reach, in that order. The largest fitting
reserve is retained. These attempts do not relax the continuous offset-fold
certificate, the source-deviation limit or published mesh-quality checks.

Shared bends previously reserved at most 40% of the shorter adjacent source
leg. On the rejected branch only, a bend may use the leg's available 90% budget.
The complete Section subsequently budgets that fixed bend together with its
other bends; overlapping reservations still fail. A small positional bow is
also no longer taken as sufficient evidence of safe curvature on a very short
endpoint span. When the curvature estimate indicates risk, both canonical sides
are checked before selecting one shared round.

These policies do not replace explicit intersection corner-radius overrides.
They govern the derived centerline bends between original Sections and controls.

## Sampling tight offsets

A recovered tight round can be geometrically valid while an overly conservative
sampling estimate creates very short quads on its inner curb. The original
estimate combined an inner curb's maximum curvature with an outer curb's maximum
speed. Newly recovered tight rounds use a bound for each actual offset instead.
Previously accepted paths retain their original row layout.

Let a unit-speed centerline be \(P(s)\), with signed curvature \(\kappa\), and
let \(N_R\) be its rightward normal. At signed lateral offset \(r\):

\[
Q_r(s)=P(s)+rN_R(s),\qquad
\frac{ds_r}{ds}=1+r\kappa(s),\qquad
\kappa_r(s)=\frac{\kappa(s)}{1+r\kappa(s)}.
\]

The unchanged fold certificate ensures the denominator remains positive.
For one row interval, the offset's arc length satisfies:

\[
L_r=\Delta s+r\Delta\theta.
\]

The endpoint tangent turn is unwrapped only when the continuous curvature bound
proves a turn smaller than a right angle over that interval. Otherwise the
interval is shortened. Bernstein bounds on the cubic's signed curvature
numerator identify positive and negative curvature possibilities. If their
absolute upper bounds are \(k_+\) and \(k_-\), then:

\[
K_r=\max\left(\frac{k_+}{1+rk_+},
              \frac{k_-}{1-rk_-}\right),\qquad
E_r\leq\frac{K_rL_r^2}{8}.
\]

Each signed branch of this error bound is convex in \(r\). Testing the two
outermost lateral offsets therefore also bounds the intervening profile rows.
The implementation includes floating-point margins and retains the existing
5 mm default tolerance. It advances by the largest certified interval up to the
nominal mesh spacing, then redistributes a short terminal remainder over a
bounded tail only if every replacement interval is certified again.

## Verification and limits

Tests cover the original Bergamo dense controls, reversal and rigid transforms,
unchanged Exact and previously accepted Smooth curves, preserved endpoint frames,
sharp-turn rejection, continuous-corridor rejection, and stored-float32 quad
quality on a tight returning loop. Dense independent samples also check the
offset chord bound on both curbs and several intermediate offsets.

Passing a path certificate is not qualification of its entire junction or
Bundle. All final profile, ownership, overlap, self-intersection and published
quad checks remain active. A short Section wholly inside a junction needs an
explicit ownership/absorption plan; removing its curvature guard is not an
acceptable substitute.
