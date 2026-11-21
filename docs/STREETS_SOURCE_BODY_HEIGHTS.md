# Original-source body height guides

The source-body helper is an **opt-in preparation API**. It supplies height
seeds and exact body constraints for an already accepted physical quad chart.
The normal surface and Product constructors expose an explicit opt-in path.
Their default remains `False`; experimental runtime revision `0.1.0.dev85`
explicitly requests source heights from ownership Tasks. This activation is not
a complete city-import qualification. Chart construction, source ownership and
actual-city progressive gates remain separate requirements. See
[the current city-import gates](STREETS_BEIJING_BERGAMO_FIXES.md).

## Source contract

[`source_body_references.py`](../way/street_network/python_backend/source_body_references.py)
retains the actual canonical Path and Profile objects. A reference covers one
positive, finite, owned original Section interval. Disjoint owned intervals
stay disjoint; the helper does not fill an unowned middle. It requires the
explicit `source_cycle_full_width_proof` interval ledger and independently
reruns the existing continuous complete-width proof before and after preparing
the guides. A 2 mm source-containment allowance is unchanged.

Full-crown references use the original Profile directly. A paired reference
requires the actual declared two-member Bundle, unique original visits and
their authored traversal directions. Both oriented half profiles must face
inward and retain the same original datum and rise. No lateral pair is inferred
from positions. The existing shared angular coordinates are preferred; a
source-directed linear frame is available only when all original offset banks
have certified positive progress in that frame.

All four original curb rays define the paired crown. Its domain is the finite
intersection of their longitudinal ranges. The helper does not extrapolate a
missing ray or supply a mixed-datum/unequal-rise policy. Dependencies include
the Bundle and every original Path and Profile used to derive its frame,
including remote member visits outside the locally guided range.

## Where a body profile can be preserved

The helper reads original corner contacts and compares **whole transverse
source cells** against actual foreign source and corner area operands.
An individual vertex lying beyond a scalar contact is insufficient: its entire
source row must clear the junction's support.

Only explicitly declared opposite-role Bundle body strips are excluded from
this influence test. Another same-role Section, a derived corner or a nearby
road remains an influence operand. Source cells require a clearance of
`2 * source_tolerance + 1e-7` metres, accounting for both sampled boundaries;
the default is 0.5001 mm. Interior influence blocks produce separate core
intervals and remain visible in the ledger.

The node-selection mask is also important. A guide node must belong to an
explicit represented source operand and must be outside every foreign operand.
This exclusion applies to **both** hard targets and soft seeds. No corner's
current surface is overwritten merely because a source strip crosses it.

The supplied `boundary_refs` must equal the actual one-face mesh boundary
vertices, including holes. Interior indices, missing boundary nodes, floating
indices and nonmanifold edge incidence are rejected. Body targets never become
drainage outlets.

## Heights and existing constraints

[`source_body_heights.prepare`](../way/street_network/python_backend/source_body_heights.py)
returns a `SourceBodyGuides` containing:

- A seed mesh with unchanged XY, quads, physical boundary XYZ, seams, materials
  and UV maps.
- Existing `HeightSlabs` rows for exact original full crowns, complete declared
  paired cores and remaining unique original half-profile cores.
- The exact target-node ledger, immutable physical fixed mask, original source
  references and dependencies.
- Continuous source-width proofs, complete-row influence witnesses, finite
  field domains, work charges and explicit unsupported scopes.

Over a measured junction-support interval, soft seeds use the quintic

$$
S(t)=6t^5-15t^4+10t^3,\qquad 0\le t\le 1.
$$

The weight reaches one only at a certified contact-free body. Within an interior
influence gap it falls to zero at the computational midpoint and rises toward
the next body. These knots are height-seed coordinates; they are not source
Ends, logical junctions, physical mouths or ownership cuts.

For original target height $h_i$, choose distinct actual fixed physical profile
anchors with heights $z_L$ and $z_H$, and define

$$
\alpha_i=\frac{h_i-z_L}{z_H-z_L}.
$$

The existing four-vertex affine row enforces

$$
z_i-(1-\alpha_i)z_L-\alpha_i z_H+0\,z_A=0.
$$

Targets must lie in the true fixed-profile envelope, so
$0\le\alpha_i\le1$. The additional distinct fixed vertex $A$ satisfies the
existing four-vertex row contract. Conflicting hard rows are rejected. An exact
unique-source body may supersede an opposite member's soft transition seed;
conflicting soft domains otherwise require an explicit junction blend.

The helper does not solve these rows. A caller must concatenate the existing
entrance rows and retain all body rows through the established bounded
`SlopeProjection`, physical-boundary-only `BoundaryDrainage` and
`PublicationSlopeGuard` passes. If later work changes XY or topology, the guides
must be rebuilt from the resulting chart. `physical_field()` constructs a
bilinear `QuadHeightField` from the final mesh only after confirming unchanged
mesh contracts and exact hard rows; it does not independently certify drainage
or slopes.

## Explicit Product construction

The maintained constructor accepts an explicit Boolean policy:

```python
product = CompoundProduct.build(context, owner, group, source_heights=True)
```

`CompoundSurfaceBuilder(..., source_heights=True)` supplies the same lower-level
path. `False` is the default; other truthy values are rejected. These options
are not a GUI setting. The runtime Task must explicitly request the policy and
declare its complete source dependencies and policy signature. A requested source
condition or body reference that cannot be supported raises its normal error;
the constructor does not silently retry using the old height conditions.

The builder resolves actual source-curb conditions first and retains that
exact `BoundaryData` tuple in the returned surface. `CompoundProduct` reuses
the tuple for physical rims and source-profile aliases, avoiding a second
boundary construction with a different policy.

The source-body stage runs after all planar refinement, shape correction and
source-cycle reconciliation, immediately before publication. It is requested
only when the group has its actual `source_cycle_full_width_proof`. The helper's
body rows are combined with existing physical-entrance rows and passed through
the existing eight-round slope/drainage reconciliation. Publication receives
both sets of rows on its initial pass and on its entrance-recovery fallback.
Its physical outlet mask stays equal to the real mesh boundary. No new solver,
chart family or resource allowance is introduced.

After publication, the builder verifies the unchanged body contract and creates
the actual physical quad query. Existing hole construction then reads the final
physical pavement rim and creates the corresponding physical grass queries.
Display exaggeration is applied afterward. Source body targets and the physical
boundary remain authoritative even at zero exaggeration. At nonzero
exaggeration, float32 publication certification can adjust free physical
heights slightly; this is recorded in `publication_roundoff`. Whole-mesh
physical byte parity between two exaggerations is therefore not assumed when
that existing correction is needed. Exact XY, physical boundary XYZ, body
targets and the final-mesh physical query are checked independently.

## Bounds and failure behavior

`BodyGuideLimits` retains finite maxima: 20,000 source samples, 8,192 source
cells, 2,000,000 operand pairs, 256 paired-field attempts, 100,000,000 charged
classification operations, 20,000 source-polygon vertices, and 50,000 physical
mesh cells. Callers can lower these bounds. Source chord tolerance is at most
0.25 mm, source containment at most 2 mm, and affine rows retain the existing
`8 * quad_count` allowance. Scratch work is charged before the corresponding
operation, including all points entering a source-polygon bounding-box test.
This also bounds cached node masks when many finite source ranges are present.
Exhaustion is a refusal, not a request for a higher global allowance.

`require_complete=True` is the default. An unsupported paired scope raises a
typed `SourceBodyUnavailable`. A review caller may explicitly use
`require_complete=False` to inspect supported rows and the unsupported ledger;
such a partial result cannot qualify the complete source-body policy.
Disjoint finite paired ranges and absence of a complete contact-free paired
core are recorded separately: no paired field is extrapolated into those gaps.

## Maintained verification and current evidence

Run the self-contained numerical tests with Blender's bundled Python in the
chat's isolated, resource-guarded environment:

```text
python -m unittest discover -s tests -p test_python_street_source_body_heights.py -v
```

The six tests cover transformed and reversed full crowns, disjoint original
ownership, whole-row junction support, seed exclusion inside foreign corners,
opposite-role-only exemptions, finite paired domains and reversed traversal,
wrong facing halves/datums/rises, remote dependencies, false physical outlets,
changed original profiles and exhausted work budgets. No generated fixture file
is needed.

[Product-constructor tests](../tests/test_python_street_source_height_products.py)
exercise the explicit source policy without monkeypatching boundary conditions:
real body targets and physical queries, a zero-height road block at zero and
tenfold exaggeration, unchanged default behavior, explicit unsupported errors,
and body-row retention through publication's entrance fallback. They also
check that publication preserves authored UVs and keeps body nodes out of the
physical outlet mask.

The five-site Beijing diagnostic was independently replayed through this helper
against its preserved physical chart. All 8,584 hard nodes and target heights
match the previous proof exactly; the final physical boundary and existing
solved hard rows are unchanged. The current helper charges 1,190 source cells,
1,203 source samples and 2,532 operand pairs, with the same 0.5001 mm clearance.
The older prototype charged 2,518 pairs; the narrower same-role exclusion adds
14 conservative checks without changing the accepted targets. This is an
opt-in helper comparison, **not** a newly constructed normal Product, a native
city import, or global source-curb policy admission.

A second replay starts from that group's original pre-body surface, whose
maximum body-target error is 11.768 mm. It uses the helper's seeds and rows,
concatenates the existing entrance rows, and runs the existing slope, drainage
and publication passes. The resulting 19,442-quad physical mesh exactly matches
the earlier independently solved mesh, including every height. All 8,584 body
targets have zero error, the combined rows have zero violation, there are no
isolated interior minima, and the maximum measured entrance longitudinal slope
is 0.8117 percent. Maximum physical and stored slopes are 9.1473 and 9.1461
percent, within their existing limits. This fixed-chart solve still does not
admit a normal city Product or resolve the remaining boundary-construction
gates elsewhere in the network.

The explicit maintained Product constructor also passes on that same five-site
group when its already accepted planar chart is supplied as the diagnostic
input. Boundary conditions, source-body preparation, physical-height solving,
publication, hole fields and Product assembly all use the maintained code.
The output retains 19,442 pavement quads, 2,566 grass quads and all 8,584 exact
body targets; physical queries at those targets have zero measured error. Its
physical XYZ exactly match the earlier solved surface. The ordinary block has
zero height and the genuine separator retains its original 0.0799999982 m
plateau, including the Product's physical rims. Only the planar chart input is
supplied; this does not qualify normal chart dispatch, an ownership Task or a
complete native city import.

Generated reports and diagnostic captures belong under
`artifacts/streets/beijing_bergamo_fixes_20261114`; the source implementation and
these tests are sufficient to reconstruct the helper and explicit constructor
path. Runtime admission still requires normal chart dispatch, source-height
Task signatures/dependencies and actual-city progressive coverage/cancellation
qualification. Earlier immutable diagnostic snapshots through source02 retained
the disabled runtime policy; their explicit constructor checks must not be
reported as genuine Task qualification. The development runtime activation in
revision 85 requires its own fresh Task and complete-city checks.
