# Sampling compound boundaries without phantom internal profiles

This document describes the Python backend's bounded correction of microscopic
clipping fragments. It is part of the ongoing
[Beijing/Bergamo qualification](STREETS_BEIJING_BERGAMO_FIXES.md), not a claim
that either full import has passed.

## The defect

Each junction footprint initially contains transverse mouths on its incoming
roads. When two junctions own an overlapping interval of an original Section,
that Section becomes internal to their compound. Its mouths must not remain as
external height conditions.

Clipping sampled curbs can leave a microscopic part of such a mouth on the
combined boundary. In the captured six-junction example, the two remaining caps
are only 58.6 and 64.5 micrometres long. Nevertheless, treating them as real
profiles makes analytical boundary sampling reject the compound. Rebuilding an
otherwise correct exterior can unnecessarily change its geometry and ownership.

Other failures are substantial: a real external profile, a larger fragment or
missing complete road width requires construction work. They cannot be removed
by this correction.

## Bounded correction

`compound_port_fragments.normalize()` runs after ordinary footprint grouping and
before complete-exterior planning. Successful existing groups are unchanged.
A candidate fragment must satisfy all of these conditions:

- Every original Section in its approach route already has an ownership claim
  in this same compound. No claim, Section or junction is added.
- The retained mouth parameter range is partial, rather than a complete profile.
- Its entire contiguous physical run is at most
  `min(1 mm, 0.2 × curve tolerance)` long.
- It lies between two existing curb intervals. Absorbed dimensional seams keep
  their dedicated construction path.

The normal analytical sampler then refines those curb contacts on the original
curves. Their radii, source paths and profile parameters stay unchanged. Every
removed cap point must remain within the existing curve-error allowance of the
refined boundary. Hole count and orientation remain unchanged, and all surviving
incoming profiles must remain complete.

Finally, continuous convex source cells prove containment of every claimed road
across its complete width, with explicit exclusion of every real hole. The
existing maximum 2 mm combined source-curve/boundary allowance is retained. A
failure returns the original group for ordinary reconstruction; it never becomes
an accepted partial surface.

The retained physical layout aggregates the original exposed approaches and
their actual outward mouth distances. Original source/end aliases stay intact.
It also requires the established drainage, height, join-slope, stored-coordinate
quad-quality and overlap checks. This operation qualifies preparation only;
surface generation and publication remain separate gates.

## Coincident original curb leads

Two footprint leads can follow precisely the same physical curb from opposite
ends of an original Section. Newton's intersection Jacobian is singular in that
case: there is a finite common interval, rather than one isolated root.

After ordinary analytical intersection fails,
`coincident_source_curbs.contact()` checks original path identity and oriented
offset. Equal coordinates of unrelated paths do not qualify.

For a route piece with original length $L$ and outward route start $b$, an
outward station $d$ maps to original station $s$ as

$$
s =
\begin{cases}
d-b, & \text{forward traversal},\\
L-(d-b), & \text{reverse traversal}.
\end{cases}
$$

The physical offset changes sign under reverse traversal. Both leads must use
the identical original path and matching oriented offset. Their finite station
domains must overlap:

$$
I = [\max(s_{0,\min},s_{1,\min}),\;\min(s_{0,\max},s_{1,\max})].
$$

A station near the two initial contact estimates is clamped to $I$ and mapped
back through each existing lead's arc-length table. The resulting points must
agree within the unchanged intersection tolerance, stay inside both finite
parameter windows and pass the existing source-approximation check. Negative
auxiliary extensions do not become original source intervals.

## Redundant ownership islands on the same curb

A short boundary run can alternate `A → B → A` even though both owners follow
the identical original Section curb. Preserving those incidental owner changes
as sampling constraints can create millimetre-long edges beside normal cells.

`canonical_source_curb_runs.canonicalize()` maps the entire B interval into A's
original lead parameterization. It requires identical original path identity and
oriented offset, continuous coverage of every station by A's finite domain,
consistent traversal direction, and matching evaluated endpoints within the
unchanged contact tolerance. The island is bounded by the current physical
sampling spacing, with at most 32 rewrites. Original path, width, radii and
incoming profiles stay unchanged.

Only ownership labels become canonical. The ordinary analytical sampler can
then distribute samples across the continuous original curve without retaining
the two incidental owner-change vertices. This is not geometric simplification
across unrelated owners. Metadata records `source_curb_owner_canonicalization`,
including the original Section, previous and canonical owner, finite run length
and endpoint error. Mesh quality and complete surface checks still apply.

## Evidence and regressions

The pavement metadata records `internal_profile_fragments`, including original
owners and Sections, the clipped parameter ranges, measured fragment lengths,
boundary displacement and the complete-width certificate. Development captures
remain under `artifacts/streets/beijing_bergamo_fixes_20261114`.

Focused tests are `test_python_street_ordinary_ports.py` and
`test_python_street_coincident_curbs.py`. They cover full external profile
preservation, source-width and hole checks, original-coordinate preservation,
reverse traversal, unrelated coincident geometry, different offsets and disjoint
finite station domains. Existing compound, source-support, Bundle and publication
tests remain required before release.
