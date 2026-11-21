# Finite dead-end pads inside a junction

This is a bounded fallback for a short, true dead-end carriageway that cannot
have an independent rounded entrance because its complete footprint is already
within the through-road junction. It keeps the authored Section and its width,
points, profile settings and source-relative anchors. It publishes one physical
surface instead of an overlapping pad road and junction.

## Eligibility and ownership

`free_pads.py` tries this construction only after an ordinary circular-contact
or finite-approach failure. There must be exactly three original approaches,
with one unique qualifying dead end:

- One straight, finite original Section, ending at a genuine free end.
- Width no greater than either through approach, and length at most half the
  smaller through width, with no Bundle membership or half profile.
- No selection cut, grade-scope end, absorbed seam or external junction.
- Two other approaches within 60 degrees of opposite directions, with a common
  original node and minimum-height datum.
- Style radii throughout; an explicit strict radius does not enter this fallback.

These eligibility bounds admit a narrow stub inside a wider bent road; they do
not establish containment. The complete rectangular footprint must fit strictly inside the
new pool, including all four boundary segments. Its full original station
interval is then owned by the junction. The ordinary road provider publishes no
second pad. Furniture still refers to the original Section and station; its
height comes from the physical junction surface.

The consumed pad has no independent incoming terminal profile. The two exposed
through-road entrances retain their complete original profiles. This distinction
is explicit in the generated `finite_free_pads` metadata.

## Smooth finite curb envelope

The two through curbs start as certified quintic Hermite curves. Mouth distances
are tried within 40 percent of each actual finite approach, with four bounded
distance candidates. The physical curb radius must remain at least the larger
of the two requested style radii.

When necessary, one curb receives a small outward bell:

$$
\mathbf C(t)=\mathbf C_0(t)+64h\,t^3(1-t)^3\mathbf n,
\qquad 0\leq t\leq1.
$$

Here the unit vector $\mathbf n$ points outward and $h$ is the maximum added
distance. The six distance candidates range from zero to one quarter of the
smaller through width. The opposite curb remains its original polynomial.
The bell has zero value, first derivative and second derivative at both ends,
so it retains the incoming positions and tangent planes.

`CurbEnvelope` represents this curve directly as a degree-six Bernstein
polynomial. Degree elevation preserves the original quintic, and only the
middle control coefficient receives the bell. Positive speed projection proves
forward progress. Restricted Bernstein bounds certify

$$
\kappa(t)=
\frac{\left|\mathbf C'(t)\times\mathbf C''(t)\right|}
{\|\mathbf C'(t)\|^3}
\leq\frac1{r_{\min}}
$$

over the complete interval, using at most 128 subintervals. The ruled map
between the two curbs also requires a complete positive-Jacobian certificate.
No triangle mesh is used.

Continuous segment containment is checked against the certified curb chords,
and the minimum distance between the complete pad boundary and pool boundary
must exceed the larger of 20 mm and four curve tolerances. The curve sampler
has a 2,048-interval limit per candidate. These checks prevent unsampled source
exits from becoming accepted geometry.

## Heights, mesh and update behavior

The existing compound boundary-jet and direct-quad pipeline generates the pool.
Boundary profile heights stay fixed. Height reconciliation, drainage, entrance
longitudinal slopes, publication precision, UVs and overlap checks remain
mandatory. Finite pools run the drainage and entrance checks even when the
initial field is already below the general slope limit.

The construction does not move an original source endpoint, extend a road past
a scope boundary or create an inferred OSM connection. A changed source is
reconstructed through the normal dependency graph and Update transaction.

### Source-labelled envelope heights

`free_pad_curb_conditions.py` validates the producer's exact `EnvelopeBoundary`
and `OwnedFreePad` record before the normal source-height rules can use it.
This is a generated polynomial through-pool curb, not an original circular
corner or a guessed native bank. Its two actual incoming profiles determine
the boundary values and inward crossfalls. In particular, ordinary full crowns
keep their real outer curb heights at the original minimum datum.

The validator checks the actual raw source-layout identity, the real three-ended
node, STYLE radius policy, canonical original Paths and profiles, complete
dead-end interval, two through-profile planes, certified degree-six envelope,
native endpoint positions/tangent directions and whole-width pad containment.
Scope/selection cuts, Bundle participants, copied envelopes, altered claims and
explicit strict corners do not inherit that certificate. The original 40-percent
finite lead bounds and one-quarter-width envelope cap remain unchanged.

The existing `_SourceRule` supplies the generated curb rule, with quintic
blending of the two real inward crossfall magnitudes. It remains a soft
generated-chain condition: real incoming profile data stay authoritative. No
original End or circular parent is invented. Existing height, drainage, slope,
quad, UV and overlap checks still apply to the resulting normal Product.

Canonical Path identity is strict, while an equal reread of an immutable Section
record is allowed. See `source_identity.matches_current_section`; a new Python
object address from scene input reconstruction is not a changed street source.

The actual Bergamo 3130 footprint/layout/surface constructors pass with 299
quads. An independent local audit checks all 34 incoming profile samples,
complete original source widths, stored-coordinate quality and slopes, drainage,
UVs, overlaps and surface queries. Stored slope is 0.06586 under 0.1 and entrance
longitudinal slope is 0.00183 under 0.03; no interior minima or drainage fill
remain. Four dedicated source-condition tests plus the existing free-pad and
dependency tests pass together (14 methods). These local results do not replace
a fresh whole-component/native import.

## Evidence and regression tests

The original Bergamo case in full source context produces a 299-quad pool with
maximum edge ratio 5.974, maximum physical slope 0.06585, entrance longitudinal
slope 0.00183 and zero drainage fill. These are focused numerical results;
whole-city qualification is recorded separately in
[the import report](STREETS_BEIJING_BERGAMO_FIXES.md).

`tests/test_python_street_owned_free_pads.py` covers the original translated
controls, complete source ownership, rotated kilometre coordinates, stored
float32 quad quality, exact incoming profile samples, scope-boundary rejection,
the envelope's endpoint jets/curvature and incremental mesh/furniture updates.
`tests/test_python_street_curb_leads.py` separately covers finite through-width
lead extension. Run them with the existing junction and ownership suites.

Generated diagnostics and source-context replays belong under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114/`.
