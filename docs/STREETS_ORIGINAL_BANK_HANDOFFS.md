# Exact original-circle bank handoffs

## Purpose and scope

`way/street_network/python_backend/original_bank_handoffs.py` handles one narrow
failed analytical contact. A complete original circular corner can already be
present, while the sampled polygon union assigns its exit to the wrong, nearly
coincident original bank. Adding another fillet would duplicate the existing
circle and can leave the new patch entirely unexposed.

This is an ownership correction. It preserves the actual original circle,
radius policy, centerline Paths, complete cross-section profiles, native bank
operands and incoming mouths. It does not accept a geometric gap, infer a nearby
bank, shrink a strict radius, increase a numerical allowance or add a solver.

Successful ordinary sampling follows its existing path. The additional work
starts only after `AnalyticalContactError` in source-bank preparation. Other
failures remain failures.

## Certificate

The candidate must establish all of the following independently:

1. The failed predecessor is a complete original `cycle_corner` circle. Its
   original Corner belongs to the group's actual physical corner inventory.
2. The two distinct original Ends belong to the actual original junction. Each
   physical Crown equals the exact native `Crown.oriented(end)` definition.
   End 1 intentionally creates a new Crown object; identity of that temporary
   object is not the cross-section contract.
3. Authored radius and mode still agree. A STRICT radius remains exact. When a
   producer `PreparedCircularInterval` exists, its original Corner/vertex ledger
   must validate and describe the whole circle without auxiliary leads.
4. The complete arc is exposed continuously from parameter 0 to 1. Its entry
   already has the exact declared first-bank handoff. Its exit equals the actual
   second Approach bank contact bit for bit.
5. A finite downstream source edge begins at that same exact stored exit
   vertex, belongs to the declared second native bank, and rejoins the unchanged
   exact other source vertex. The intervening run belongs only to the declared
   first bank. At most eight source spans are inspected.
6. Both the false owner and the intended owner have validated finite native
   `CurbCoordinates`. There is no coordinate-only or nearest-bank assignment.
7. The entire piecewise-linear run lies in the intended finite chord's
   parameter range and the existing `RegionUnion.epsilon` welding band.
8. No unrelated source operand's bounding box intersects the finite changed
   run. This conservative condition refuses cases requiring another ownership
   analysis.

For a source chord from A to B, put D = B − A. For every actual polyline knot P
of the false run and the existing union run, the certificate evaluates

```text
u(P) = dot(P − A, D) / dot(D, D)
d(P) = norm(P − A − u(P) D)
```

It requires `0 <= u(P) <= 1` and `d(P) <= RegionUnion.epsilon`. The range and
distance bounds hold over each complete linear interval: u is affine and
distance to the fixed chord line is convex. This is a whole-run bound, rather
than a midpoint or sparse visual check. The band only qualifies the already
proved ownership ambiguity; it never establishes source identity or permits
the final analytical sampler to leave a gap.

## Construction and final guards

The normalizer replaces the false-owner run and the already retained partial
exit edge with the complete intended original edge. Both endpoints are existing
source vertices. All original footprint objects, their arrays and their curve
ledgers remain unchanged; obsolete union clipping vertices can leave the outer
boundary while remaining present in their original operands.

The ordinary analytical sampler runs again. Only after it succeeds does the
pipeline accept a fresh continuous 2 mm native-width proof, unchanged hole
topology, complete incoming profiles and unchanged mouths. The existing future
source-strip overlap check runs with `max_candidates=0`; any required mouth
advance is rejected. No more than sixteen failed handoffs can be normalized in
one sampling operation. Original 20,000-sample and existing proof budgets remain.

`SourceBankReconciler` retains the corrected sampled source base even when no
additional bank fillet is needed. Its final sampling after other bank infills
uses the same guarded path.

## Metadata and invalidation

`physical_layout.source_bank_handoff_canonicalizations` contains local
construction certificates and validation reports. `source_bank_reconciliation.signature()`
includes each record's `construction` token even when no new fillet records
exist. Diagnostic counters and proof-scan counts do not define that local token.

The policy is `SOURCE_BANK_POLICY =
'2-typed-original-banks-with-exact-exit-handoffs'`. Unchanged cases keep the old
signature tuple shape; a handoff-bearing group adds its local construction
records. Existing PLAN/layout/surface dependency closure still supplies actual
Path, Profile and ownership-range dependencies. Global prepared-neighbor and
normal Task validation remain required.

## Verification and resuming

Maintained tests are in
`tests/test_python_street_original_bank_handoffs.py`. They construct a small
native three-arm network from source, including source reversal and STRICT
radius preservation. Negative cases cover changed profiles/Paths/radius,
partial arcs, unrelated third operands and runs outside the numerical band.
The handoff-only signature test prevents a missing dependency token when no
new bank fillet is generated. No development pickle or checkpoint is needed to
regenerate these tests.

The actual Bergamo 1,255-Section failure was captured from a fresh normal PLAN
attempt. The original radius-2 m corner's 0.33648 mm arc had a 7.076 nm sagitta.
The union retained a different bank's tail at its exit; the maximum certified
whole-run chord distance was 62.885 nm inside the unchanged 100 nm welding band.
The proposed extra bank fillet was bit-identical to the existing original
circle and had no exposed arc.

An explicitly labelled immutable-runtime algorithm experiment then passed the
normal local constructor and normal Product with 2,865 quads, all 26 hard body
targets and 51 full incoming-profile samples exact, zero isolated minima and
unchanged physical/publication quality guards. It supplied no chart, mesh,
cuts or boundary conditions. That experiment included the mapped-rule fix
documented in [Normalized source heights](STREETS_NORMALIZED_SOURCE_HEIGHTS.md).
It was not a genuine whole-component or whole-city mesh qualification.

Fresh normal PLAN and native import gates must be checked against the final
coherent runtime source after integration. Development reports live under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114`; they are
evidence, not required recovery inputs. Do not restore numerical products or
restamp old Task digests to bypass those gates.
