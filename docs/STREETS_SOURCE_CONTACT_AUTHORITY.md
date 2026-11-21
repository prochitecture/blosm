# Original corner authority and displaced physical contacts

An intersection's corner-radius setting belongs to an ordered pair of its
original Section ends. A physical circular corner can nevertheless reach a
different Section through roads wholly owned by the same compound. The later
Section keeps its actual width, profile, native path and endpoint identity.

These are two distinct facts. Requiring both physical contact ends to be at the
original node rejects legitimate displaced contacts. Relabelling the later end
as an incidence of that node would create false topology. The implementation
preserves the original logical corner and validates the physical source walk.

## Code and dispatch

`way/street_network/python_backend/displaced_corner_authority.py` supplies:

- `validate_displaced_corner(context, group, node, corner)`, returning immutable
  `DisplacedCornerAuthority` with the real logical ends, physical contact ends,
  source walks and attempted original radii. Invalid records raise
  `DisplacedCornerAuthorityError`.
- `original_parent_corner(resolver, node, corner)`, retaining the existing direct
  `junction_bank_transitions.original_corner` check and additionally admitting
  the explicitly validated displaced producer. A displaced check requires
  `resolver.group` as well as `resolver.context`.

The shared parent dispatcher also accepts the separately validated
`OwnedPrefixRoute` producer. Its real logical radius end and actual later
physical carriageway remain distinct; the original native order and the
producer's physical proxy order are recorded separately. See
[owned-prefix corner authority](STREETS_OWNED_PREFIX_CORNER_AUTHORITY.md).

The general `original_corner` function remains strict. A consumer requiring
actual original-node incidences must continue to use that function. Retained
circle geometry and its height conditions can use `original_parent_corner`
because they need original radius authority and actual physical source profiles,
without requiring the physical contact end to be incident at the logical node.

## Required provenance

The additional branch accepts only a corner produced with
`DISPLACED_ORIGINAL_DIMENSIONAL_SOURCE_CURBS`. It requires the producer's two
typed `SourceContactWalk` records and a genuine width-changing owned prefix.
It checks all of the following before accepting that record:

1. The original node is the actual network object and belongs to the compound.
   Both logical ends belong to it, in their original counterclockwise order.
2. Logical and physical arms are ordinary original `Approach`/`Route` objects.
   Every route piece retains the context's identical native `Path`. Widths,
   offsets and oriented profiles equal the original Section's values.

   `source_identity.matches_current_section` requires that canonical Path's
   exact Section object and exact authored-value equality to the current input.
   This permits equal scene rereads to reuse cached native resources while
   continuing to reject copied Paths, stale geometry, dimensions, profiles or
   source membership. Pointer equality to a newly allocated input Section alone
   would incorrectly invalidate such a reread.
3. The source-incidence walk is recomputed using the producer's existing
   bounded `_walks` routine. Every crossed road is consumed by the compound;
   real junction incidence and the original bank side choose each continuation.
   The existing scope, selection, profile and strict-sector exclusions apply.
   The recorded piece directions and native path identities must match.
4. The requested radius and mode come from the logical node. This producer
   supports only `STYLE`. Its attempted radii must be an exact prefix of the
   existing `CornerRadiusPolicy.contact_candidates` sequence, ending at the
   published radius. No radius reduction or alternative contact search occurs
   in this validator.
5. Both contact stations are finite, inside the actual arm and away from
   internal route seams. Their native directed intervals must lie within the
   compound's actual disjoint owned ranges. The test does not take a hull over
   a gap in ownership.
6. The original circle endpoints coincide with the actual native offset curbs,
   preserve the radius and remain tangent to those curbs in the original
   circular sector. Existing numerical tolerances remain in force.

## Contact equations

For the physical first and second arms, let their actual offset curves be
`B₀(s) = P₀(s) + (offset₀ − width₀/2) N₀(s)` and
`B₁(s) = P₁(s) + (offset₁ + width₁/2) N₁(s)`.
The accepted circle has centre `C`, radius `R`, start angle `α` and sweep `θ`:

$$
A(t)=C+R\begin{pmatrix}\cos(\alpha+t\theta)\\\sin(\alpha+t\theta)\end{pmatrix},
\qquad 0\leq t\leq1.
$$

At the recorded actual stations `s₀` and `s₁`, validation checks:

$$
A(0)=B_0(s_0),\quad A(1)=B_1(s_1),\qquad
\lVert B_i(s_i)-C\rVert=R,\qquad
(B_i(s_i)-C)\cdot B_i'(s_i)=0.
$$

The endpoint, radius and normalized tangency residual limits are respectively
`10⁻⁸ m`, `10⁻⁵ m` and `10⁻⁷ m`, consistent with the existing circular contact
and prepared-circle checks. The earlier complete source-width, topology,
sampling, full-mouth and final prepared-neighbor gates still apply afterwards.
This authority check does not replace any of them.

## Verification and remaining qualification

`tests/test_python_street_displaced_corner_authority.py` verifies logical versus
physical incidence, reversed storage, rigid motion, unchanged source controls,
and refusal of false prefixes, copied native paths, unowned intervals, changed
radius policies, wrong source nodes and changed circular geometry.

The saved real Bergamo failure has logical ends `164/0 → 163/0` at node 333,
but physical contacts `164/0 → 165/0` through owned Section 163. The validator
accepts that exact record, preserving its original STYLE sequence from 10 m to
0.524288 m. The direct original-node check correctly remains false for that
physical pair. This local authority result alone does **not** qualify the
component's surface, native import or progressive user experience; those gates
are tracked in [the import qualification report](STREETS_BEIJING_BERGAMO_FIXES.md).
