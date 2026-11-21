# Corner authority through an original owned prefix

[`owned_prefix_corner_authority.py`](../way/street_network/python_backend/owned_prefix_corner_authority.py)
validates the existing `OwnedPrefixRoute` producer. It changes no source path,
profile, radius, physical curve, mouth or mesh. It does not broaden the direct
`original_corner()` predicate and does not accept arbitrary Route subclasses.

## The two identities

An original short Section can be wholly inside its real junction while a later
Section supplies the actual incoming carriageway frame. The route retains the
short Section's logical End for its original junction and radius policy; its
physical width, offset and profile belong to that later Section.

In the fresh native Beijing failure captured by
`runtime30_beijing691_closure`, the logical sector is `402:end0 → 401:end0` at
junction 802. Its 4 m STYLE circle physically contacts original Section 402 and
original `998:end1`. Source 998 is not an End of junction 802. The genuine typed
claim is:

```text
junction: 802
logical prefix: 401:end0
whole prefix length: 0.760006947 m
real two-ended seam: connection 162 = {401:end1, 998:end1}
actual external contact: 998:end1
```

The logical prefix is 10.8 m wide with a full crown; the external carriageway is
5.4 m wide with its own oriented half-crown. Substituting the logical profile
for that actual contact, or pretending that source 998 belongs to junction 802,
would both be incorrect.

The producer orders its actual Approach frames. A prefix replacement can
therefore change the geometric cyclic order while maintaining a one-to-one
mapping to the real original Ends. The validator records **native order** and
**proxy-layout order** separately. It verifies the selected sector against the
actual source layout used by the producer and obtains its requested radius/mode
from `node.corner_radius(first_logical_end)`. It does not claim the physical
proxy order is the order of the native prefix tangents.

## Mandatory checks

`validate_owned_prefix_corner(context, group, node, corner)` returns a typed
`OwnedPrefixCornerAuthority` only after checking:

- The real current junction and both logical Ends.
- The explicit `EXPLICIT_WHOLE_OWNED_PREFIX_CONTACTS` producer record and two
  real `SourceContactWalk` objects, including the actual prefix centreline.
- The exact canonical prefix Path, actual stored claim object and original
  two-ended seam. A copied claim or a different same-geometry Path is rejected.
- The unchanged external directed topology trace, canonical native pieces,
  physical widths, offsets, actual profiles and finite Route metric.
- Complete original prefix and external mouth aliases, original support,
  actual contact stations, the producer's proxy-frame sector orientation and
  finite circular tangencies.
- The original STYLE/STRICT candidate sequence. A STRICT radius remains exact.
- Fresh continuous whole-width containment, including real holes, under the
  existing 2 mm source-error bound.

Cached native resources may retain their original Section objects when an
equal Network is read again. The validator uses
`source_identity.matches_current_section()` to retain canonical Path identity
and exact immutable Section values, rather than requiring the newly allocated
Network Section to have the same Python address. A borrowed Path remains invalid.

## Native ranges without recursion

`verified_owned_prefix_ranges(context, group, arm, station)` is the narrow
admission API for source-range callers. It returns `OwnedPrefixRangeProof`:

```text
ranges                    exact native intervals for the requested station
requested_station         finite actual external Route parameter
mouth_station             unchanged final incoming profile plane
final_mouth_ranges        whole prefix plus actual external native intervals
final_mouth_width_proof    fresh continuous source support certificate
```

The helper never calls `source_bank_reconciliation.owned_ranges()`. It checks
the actual final Approach and mouth aliases directly, proves that the prefix
is wholly owned, and repeats centreline and complete-width support tests. The
requested contact ranges must be subsets of the certified final-mouth ranges.
Only a successful typed proof may contribute these intervals to the common
owned-range ledger; generic Route-subclass acceptance is not permitted.

For the actual Beijing source, the true final intervals are:

```text
401: [0, 0.760006947]
998: [34.255326579, 47.511052398]
```

At the actual corner contact, the 998 interval is the smaller
`[38.363606445, 47.511052398]`. Both independently pass complete-width proof.
The policy token is `SOURCE_OWNED_PREFIX_CORNER_POLICY`.

## Physical port bank coordinates

`owned_prefix_source_curbs.py` supplies the narrow metadata/source-condition
contract for a port whose exact type is `OwnedPrefixRoute`. The producer already
knows every outward bank station, side and full evaluation batch; coordinate
policy 6 retains that declaration. It does not reconstruct stations from XY or
accept arbitrary Route subclasses. The existing bank/cap vertices, source
labels and sampling-error budget remain byte exact.

The native pieces refer to the actual external Route, not the logical prefix.
For the reversed source 998 in this example, an outward station `s` and bank
side `q` map to:

```text
native_station = path_998.length - s
native_side    = -q
```

At real native Route seams, the existing explicit piece map splits a chord's
ancestry without adding or moving a vertex. Source 401 remains a separate whole
ownership witness. Its wider full-crown profile never supplies source 998's
physical half-profile height or gradient.

Before the existing `_NativeCoordinateRule` is used,
`validate_curve(resolver, source, index, curve, first, last)` checks the exact
port, Approach, original evaluation batch, actual canonical native pieces and
finite bank range. It repeats the current-group final-mouth/whole-prefix proof.
That expensive proof is shared only within the current immutable condition
resolver, group, exact Approach and complete profile plane. Every physical
curve still validates its own exact original points and source identities.
No proof is carried into another resolver or Update.

Legacy port-side inference remains unchanged and refuses missing/ambiguous
metadata. The native-coordinate rule supplies true external profile heights,
native gradient covectors and compulsory native-piece breaks.

## Evidence and remaining qualification

Reports are under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114/`:

- `runtime30_owned_prefix201_report.json` establishes the actual claim,
  topology, missing old range entries and independent contact/final-mouth
  support. The initial generic range code skipped this explicit subclass.
- `runtime30_owned_prefix_authority202_report.json` passes the new typed
  authority on that genuine fresh failed group, including an equal reconstructed
  Network with canonical cached Paths.
- `owned_prefix_authority_tests205_report.json` passes five portable tests for
  true half-profile contact identity, STYLE/STRICT policies, rigid motion,
  reversed storage, reconstructed inputs, and refusal of copied claims, borrowed
  Paths, missing aliases, exceeded mouths, changed radii and changed geometry.
- `owned_prefix_authority_tests206_report.json` repeats all five tests with
  explicit rejection of a changed proxy-frame gap or reversed circular sweep.
- `runtime30_owned_prefix_authority207_report.json` repeats the genuine source
  proof after that orientation check: all 36 contact-support cells and all 30
  final-mouth-support cells pass the unchanged 2 mm bound; an equal reconstructed
  Network again retains authority through its canonical cached Paths.
- `owned_prefix_authority_tests209_report.json` passes six portable tests after
  the public parent-authority and owned-range dispatch are integrated. The
  strict direct-original predicate remains false; the typed parent wrapper and
  exact original prefix/external ranges succeed without recursive range lookup.
- `genuine802_normal_product210_report.json` regenerates normal source-cycle
  geometry from the genuine runtime30 raw Task resources, with no supplied
  group, mesh, chart or conditions. Fresh preparation passes in 14.25 seconds,
  including one retained-circle transition at 0.65536 m. The subsequent normal
  source-height Product refuses a missing native port ledger, after 220.54
  seconds; this is not a Product qualification.
- `genuine802_legacy_port212_report.json` isolates that remaining refusal to
  `cycle_port:802`, curve 0. The exact producer station map is available, but
  the source-coordinate constructor rejects its explicit `OwnedPrefixRoute`
  type and omits the ledger. The physical profile belongs to source 998. A
  source-height remedy must retain that true native identity and declared
  station map rather than relax the legacy profile test.
- `owned_prefix_authority_tests214_report.json` passes 24 focused tests with
  coordinate policy 6: exact original evaluation batches/subsets and forged
  reference refusals, ordinary bank geometry parity, typed prefix metadata
  parity, reversed/rotated native half-profile values and gradient covectors,
  reread identity and portable authority/range integration.
- `genuine802_conditions215_report.json` regenerates preparation and normal
  source conditions from the genuine runtime30 raw resources. It passes in
  51.28 seconds (15.26 seconds preparation, 28.34 seconds conditions), with
  unchanged incoming profile values/gradients and zero height jump. Field fit
  error is 3.47e-11, the maximum boundary slope is 0.064 under the 0.1 limit,
  and the true height range is bounded by datum 0 and ceiling 0.08 m.
- `owned_prefix_authority_tests218_report.json` repeats all 24 focused tests
  successfully on Python 3.10; test 214 uses Blender's bundled Python 3.13.
- `genuine802_normal_product216_report.json` passes a fresh normal preparation
  and unsupplied source-height Product/chart constructor: 7,225 pavement quads,
  no grass, 14.97 seconds preparation and 521.61 seconds Product time (543.35
  seconds total including source capture). Peak resident memory is 261,627,904
  bytes, below the unchanged 3 GiB cap. The bounded normal chart recovery is a
  substantial runtime limitation; this is not evidence of fast city import.
- `genuine802_normal_product217_independent_report.json` independently passes
  in 12.77 seconds: stored ratio 7.67894, minimum corner sine 0.200452, maximum
  physical slope 0.066233 under the 0.1 limit, entrance longitudinal slope
  0.004397 under 0.03, zero drainage fill and no isolated interior minima.
  All 102 true native profile samples retain exact XY/Z, all 414 continuous
  source-width cells pass 2 mm, UVs are finite and overlap checks pass. Query
  error is at most 4.67e-15 and released-copy queries remain exact. The
  independent original-body scope is complete with no unsupported records;
  its **nine eligible hard body nodes** match true source targets exactly.

The authority/support gates, fresh normal source preparation, formal source
conditions, normal unsupplied Product and independent physical/body/profile
audits pass. The complete fresh native coordinator/Surface Task pipeline still
needs qualification; these local results must not be relabelled as a completed
city or GUI import.
