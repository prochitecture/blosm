# Original Bundle chains and ground extents

The numerical backend receives two ordered lateral paths, not an unordered
inventory of every Street heuristically associated with a native Bundle. The
mapping keeps every original Section, Street and junction identity. Coordinate
proximity never establishes incidence.

`bundle_source_chains.py` handles two complete opposed, named one-way source
paths fragmented across original Streets. Both source paths must have explicit
original endpoint incidence, and their original control outline must preserve
lateral order. A shared original terminal is evaluated using the actual local
junction trims before numerical centreline separation is certified. Its zero
authored endpoint separation does not by itself establish a crossing.

`bundle_main_chains.py` selects a unique highest road-class pair from the native
head and tail inventories. Each continuation must use an original incidence and
either its source-way continuation or a unique aligned source direction.
Ambiguity may be resolved by bounded reachability to an unused declared tail.
An original `_link` can belong to the same longitudinal path. Other owned roads
keep explicit interior/exterior roles from the complete source outline; none is
lost when the main path is normalized. A named legacy fan without a same-class
tail can use complete source endpoints only under the additional bounded
`FULL_NAMED_SOURCE_EXTENT` proof.

An unmapped Bundle still owns its selected ground roads. The Python adapter
uses the original owned/head/tail inventory and Street ownership to report
`UNMAPPED_BUNDLE_MEMBERSHIP` and block those roads. Empty mapped-member lists
cannot silently turn a genuine paired corridor into unrelated full crowns.

## Disconnected ground extents

Grade removal constructs the induced original ground graph. A bridge or other
excluded Section stays outside that graph, and each retained handoff uses an
existing original `scope_ends` incidence. A member path is never concatenated
across the removed span.

`ground_bundle_scopes` can return several disconnected paired extents. Original
ground runs must have positive finite axial overlap, strict lateral order, a
one-to-one pairing, and matching run order. A run with two possible partners is
rejected. Original controls establish selection only: actual curved paths,
finite caps, whole profiles, surface quality, UVs and final all-surface overlap
checks still qualify each generated extent. Unpaired ground runs remain source
roads and are listed explicitly in the scope report.

One finite extent keeps the native Bundle ID. Multiple extents use reserved
runtime task IDs:

```
(1 << 62) | (source_bundle_id << 16) | scope_ordinal
```

The source ID is a nonnegative INT32 value, the ordinal is in `1..65535`, and
the resulting task ID is positive and below `2^63`. Adapter collision checks
cover all native Bundle IDs. Numerical `Bundle.source_bundle_id` and
`Bundle.scope_ordinal` preserve the canonical lineage in the Network JSON and
published Part metadata. These runtime task IDs travel as Keys or strings;
they must not be written into Blender INT mesh attributes. Canonical source
IDs and ordinals remain separate INT32-compatible fields. Older Network rows
without lineage fields still load as ordinary Bundles.

## Interior roads after grade removal

A native interior Street can contain ground approaches separated by an excluded
bridge. The adapter inventories each retained original run separately. A ground
crossing belongs to an extent only when original physical incidence supplies
both opposite-member ports. A run with one member port or competing extent
ownership remains unsupported.

A retained run with no member incidence can be classified as exterior only
after its complete original control path is separated from every finite paired
control corridor, including the full original width allowance. Such runs keep
their complete ground road geometry, original scope handoffs and ordinary
profiles. Their scope records distinguish them from excluded bridge pieces and
state that numerical geometry qualification remains required. A disconnected
run inside the finite corridor is rejected instead of being omitted from the
separator inventory. Normal final pavement/grass overlap validation applies to
all retained exterior roads as well as the paired extent.

The regression suite is `tests/test_python_street_bundle_scope.py`; it covers
multiple disjoint extents, opposite traversal direction, complete half
profiles, retained ground approaches, intrusion/ambiguity rejection, reserved
identity ranges and collisions, and old/new Network JSON round trips.

## Physical source-cycle guards

An ordinary side-junction compound can retain an internal passage through one
complete Bundle member while consuming only an unbundled connecting road.
`retained_bundle_through.py` requires two consecutive original member
traversals, both unconsumed ports at the same owned physical junction, and
unchanged complete oriented half profiles, widths and offsets. It grants no
Bundle body interval and creates no new connection. The usual physical
compound, cap, complete-profile and separator checks still apply.

A merged entrance may sit inside a multi-Section Bundle run. Its finite local
field uses the two **crossing-adjacent** original Sections: the first visit at
a run head, the last visit at its tail. `PairedAngularEntrance.run` retains the
complete original `BundleRun`, and `source_run_sections` retains every member
Section for body ownership and dependency closure. The two local Paths do not
replace that whole run. Source mouths keep their actual Section/end/station
aliases and original half-profile parameter maps. Full Bundle body dispatch
must retain every other original Section and internal junction attachment;
local merged-entrance eligibility alone does not qualify that body.

`merged_run_inventory.py` partitions each complete run at every original
multi-arm junction. Only explicitly connected compatible two-ended seams may
remain in one source span. Complete incoming samples come from accepted
physical products; standalone inner curbs are walked on the actual published
mesh boundary, and compound inner curbs use the accepted physical curb chain.
The original whole-Path sampling phase remains intact, including authoritative
finite-prefix samples. A finite terminal's local chart is kept separately from
its independently certified whole-run frame. No local finite field is
extrapolated over a complete run.

Paired side-junction owners restrict **evaluation intervals**, preserving the
original Section and global metric station. `evaluation_spans` excludes only
the explicit owned opposite interval and retains all remaining original
visits. Its derived boundaries have no original `End`; their separate
`DerivedSourcePortKey(owner, section_id, direction)` identifies a physical cut.
The ordinary endpoint resolver is never used for these interior cuts. Both
the visible intervals and the owner's intervals must cover the complete
original metric support. A one-sided ordinary full profile which still
overlaps the opposite member is rejected until paired physical preparation
exists.

`hydrate_evaluation_spans` resolves each leaf against accepted physical
products. Original endpoints use the original registry; interior cuts use
only their owner-scoped derived registry. Each merged end records its complete
`PreparedMergedMouth`, physical approach index and exact source station.
Active `source_intervals` stop at these real mouth aliases. The two curb
chains retain the complete finite Paths of the leaf's original visits as
`curve_support`, so a full physical profile plane and its inward connector
have finite support beyond the centre station. This reserve grants no extra
pavement ownership: body output still stops at the accepted physical mouths.
All overlapping prepared prefixes must share byte-identical original station
partitions, and sampling adds no independent knots inside those prefixes.
Dimensional transitions retain their actual source curves and owners. The
helper does not apply a whole-run separation assumption to an already paired
side-junction owner; normal body, stored-surface and overlap checks qualify
each resulting piece.

Source-cycle pools keep the actual original carriageway regions and local
fillets, including physical median holes. A missing internal fillet can be
unexposed in one filled local corner while a different part of the pool still
has a median hole. `compound_cycle_hidden.py` checks each complete facing
original curb continuously against the pavement exterior and every hole. Any
positive-length contribution to a hole rim, including coincident ownership
ties, remains unsupported. The proof is recorded in
`source_cycle_covered_corners`; a global overlap observation does not replace it.

Small positive contacts on curved, nearly opposite approaches can leave a real
sliver of an original root cross-section outside an otherwise valid circle.
`corner_source_support.py` verifies these candidate contacts using the local
fillet patch and finite original road strips at that same source junction.
Only the candidate pair's root caps are acceptance targets; a third original
road can supply real coverage without changing the preferred radius. A failed finite-root check
uses the existing radius policy acceptance callback: STYLE may try its next
already permitted radius, while an explicit STRICT radius remains strict.
This does not expand epsilon, change widths, or accept a made-up closing curb.
The final union still rejects every exposed internal preparation seam, and
complete source-strip ownership and final mesh qualification remain required.

Secondary curved-contact candidates can request `audit_all_contacts=True`.
Their original road strips still start at station zero; negative stations
belong only to the bounded auxiliary fillet leads. An additional explicit
`audit_facing_halves=True` audits the candidate's own transverse sectors:
first-left and second-right. It intersects the full retained cap intervals
with those sectors, rather than testing a midpoint. Other sectors still need
their own prescribed local corners and the final whole-source union proof.
This sector check changes no incoming carriageway profile. Audit sampling
coalesces near-identical stations by the existing canonical `1e-9` station
rule, retaining exact zero, finite end and contact stations where distinct.

`compound_dimensional_corners.py` handles a centred near-through width step
whose accepted corner requires a certified dimensional continuation. The
quintic anchors at the **wider original root curb**, spends its bounded reach
on the narrower original source, and joins both actual finite curb leads. A
symmetric taper would cut the wider original cross-section. The helper keeps
the supplied radius as the full-curve minimum radius, uses at most sixteen
finite lead increments, and rejects scope/selection or real-junction
extrapolation through the ordinary mouth checks. Circular STYLE replacement
requires the caller's explicit failed root-support proof and existing
near-through policy; the helper does not replace a successful circular turn.
Its internal closure follows the two original centrelines, while physical
curb labels retain the actual certified curve. All complete source strips,
incoming profiles and final union guards remain required.

An original Bundle half profile can opt into that same dimensional helper only
after its actual facing root is exposed in the complete source union. The
nearly straight STYLE pair must meet at two zero-height outer curb endpoints;
the other plateau edge is ineligible. Original full widths, source-profile
halves and endpoint fields stay unchanged. The explicit
`original_outer_half_profiles` flag leaves ordinary helper calls and STRICT
half-profile rejection unchanged. The full source-width, hole, complete
incoming-profile and surface checks still qualify the resulting patch.

When all permitted circular STYLE contacts erode an actual exposed original
root, `compound_root_envelopes.py` can retain the preferred minimum radius with
a finite certified curve. It uses the actual source endpoint jets and adds a
bounded outward Bernstein bell which preserves the base curve's endpoint
jets. Both complete original root cross-sections must become internal. The
record reports the outward distance, its bound, finite stations and certified
curvature. STRICT radii remain circular. Original paths, widths, source
profiles and finite approach bounds remain authoritative, followed by the
complete source union and all surface/publication checks.

An already certified near-through width transition may need to retain an
additional original root at the same source junction. `prepare_adjacent`
requires that actual exposed-root diagnosis, the accepted dimensional record,
and unchanged local circular corners. It fixes the wider root endpoint and
spends finite reach only on the narrower source. A bounded outward envelope
must satisfy the original minimum radius and make every participating full
root cross-section internal to the local original strips and corner patches.
The neighboring circular corners, source half profiles and source controls
remain authoritative. This policy does not replace an acute corner by a
width continuation; it adjusts the adjacent existing dimensional curb. The
complete source union, holes, incoming profiles, heights and publication
checks still qualify the result.

`cyclic_root_support.py` extends that explicit exposed-root policy to a complete
authored cyclic incidence set. After the established circular and dimensional
families fail, an adjacent full-profile STYLE curb can support a wider original
zero-height outer half. The actual wider source sets the finite reach scale;
the sixteen-increment, original minimum-radius and outward-displacement bounds
remain required. The helper verifies original approach object identities,
preserves the other circular corners, and checks every complete root cap by
geometry regardless of boundary-owner tags. Its whole new analytical curb must
remain exposed in both the local proof and final source union, preventing a
clipped crossing from becoming an unrounded outer join. Plateau-side and
authored STRICT donors remain ineligible.

The source pool first retains its existing thirty-two complete-profile lead
attempts. `finite_profile_search.py` can then inspect at most 4,096 already
certified original-Path stations, with at most thirty-two further union rebuilds.
Only a profile actually clipped by the preceding full union can advance. A
continuous whole-segment test excludes every other current source operand.
If a rebuilt prefix clips another mouth, complete finite original strip cells
serve only as additional exclusion proofs, with two original chord-error
allowances plus numerical clearance. Those distant strips never enter the
physical union or acquire ownership. Original finite scope limits and known
opposite-owner closure remain mandatory. All original widths, complete mouths,
root seams, holes and final surface checks are repeated after selection.

An exhausted ordinary STYLE nose can be smaller than the usual 0.5 m fitting
heuristic when two complete original smooth roads enclose a genuine narrow
median. `thin_median_corners.py` requires both authored physical junctions,
full original profiles, continuously monotone facing banks and a certified
bank width too small for the old minimum circle. It then tries at most two
geometry candidates, 0.25 m and 0.125 m, in descending order, capped by `A/P`
of the original-strip hole. Here `2A/P` is a characteristic separator width;
the nose diameter cannot exceed it. This explicit style rule selects the
radius during source preparation, before any surface exists. Ordinary
successful contacts and authored STRICT radii retain their earlier behavior.
The final normal union must retain the true hole and expose each complete
tangent-to-tangent inner circle. Full-width, incoming-profile, mesh, height,
drainage, stored-coordinate and overlap checks remain mandatory.

The captured narrow Bergamo lens has about 5.361 square metres of original
area and 30.154 m perimeter, giving a radius cap of 0.1778 m. Its 0.125 m noses
retain 4.899 square metres, about 91.4 percent. The area reduction is intentional
rounded-end pavement geometry. It is recorded along with the original area,
perimeter, selected radius, bank certificate and retained-hole area; it is
neither numerical hole loss nor a filled separator.

If a retained original curb cannot represent overlapping endpoint fillets,
one compatible lead, or a failed original circular contact, the typed
`RetainedCurbOverlapError`, `RetainedSourceRouteError`,
`RetainedCornerContactError`, or `RetainedSourceWidthError` permits the same
fully guarded finite source union.
It does not catch arbitrary construction errors. Whole-Bundle graph admission,
prefix ownership, complete source widths, all external profiles and partial
Bundle claims still qualify that result. A failed contact retains its original
junction, incidence pair, radius, mode and underlying contact error. The union
must certify the complete original source widths; an omitted internal corner
must additionally pass its original facing-curb coverage and root checks.
Before accepting a rebuilt exterior, continuous finite source cells prove
the whole consumed roads and every claimed incoming prefix within the existing
2 mm bound. A raw ordinary footprint union with complete visible mouths still
enters original-incidence preparation if it clips a consumed road's width.
This is a failed-proof planning trigger; it cannot qualify the replacement.
The finite source union must repeat the complete width and prefix proof along
with its original profile, root, hole and neighboring-owner checks.
Deferred ownership closure keeps eight sites as its ordinary resource heuristic.
Its restarts are bounded by actual joins in the finite original source forest,
rather than eight passes shared by the entire network. An expansion must name
an available original physical neighbor and grow its component. During this
one immutable planning context, an unchanged member tuple reuses its already
certified exterior; failed attempts are not cached, and obsolete smaller
tuples are released. Final owner ordering and every geometric guard remain
unchanged.
`source_pool_budget.py` admits a larger connected pool only with a fresh
complete original-source certificate and measured work: at most sixteen
original physical owners, 256 operands, 20,000 sampled segments, 20,000
continuous width-proof cells, sixteen complete profiles and eight physical
holes. Every consumed road retains both original owners and its whole finite
interval. The certificate records the source identities, work counts and
geometry extent. The downstream 50,000-quad allocation and all source,
profile, height, stored-coordinate and overlap checks remain required.
Only this newly admitted input class can select the original-source corridor
atlas before the older eight-site chart families. It needs the already
attached `MEASURED_COMPLETE_ORIGINAL_SOURCE_POOL` certificate, atlas eligibility,
and exact agreement with a fresh work certificate and complete-width proof.
The admission check uses a copied layout, preserving the accepted PLAN data.
Ordinary groups without that certificate and every group of at most eight
owners keep their established route. Candidate corridor ordering, all physical
rims and claims, the 50,000-cell cap, and every final surface guard stay the
same; this changes construction priority for previously unsupported large
source pools only.
Source pools also try their existing
coupled slope/drainage reconciliation when an initial isolated slope
projection creates a basin, before exhausting spatial refinements of the
unchanged analytical field. The later drainage and join checks remain required.

If a compound consumes a complete short Bundle member, a partial claim cannot
deplete that member while leaving the opposite original terminal outside the
owner. `bundle_terminal_expansion.py` identifies only the missing original
opposite terminal junctions, provided both consumed-member terminals already
belong to the compound. These must be available physical sites, without a
scope/selection cut or a skipped unowned attachment. The ordinary eight-site
budget and source-cycle construction then decide whether the expanded whole
pair qualifies. Contiguity, complete-width, separator and overlap guards are
unchanged; no synthetic connection is created.

After the established physical-hole and original paired-interval families
fail, `source_cycle_atlas.py` can choose normals from at most three of the
longest actually consumed, unbundled original paths. All physical source
claims, curbs, profiles and holes remain unchanged. These are computational
axes; they do not assign lateral membership to roads. Optional visibility
endpoints include existing internal cut corners. Their diagonals extend on
the same straight line only to the first physical curb, with continuous
pavement and hole checks; incoming profiles remain exact. The ordinary
visibility angle, kernel, candidate and cut budgets are retained.

This last atlas retains the established energy attempts first. A best
rejected candidate with at most 128 ratio failures can receive at most four
further 500-step energy seeds and the existing bounded correction, always
with the physical rim and its vertex IDs fixed. If separator-rim attachment
exhausts its local guide and refinement family, a separate 300-step global
correction can release interior guides while keeping the whole physical
rim exact. Original and stored quad quality qualify both stages before the
normal height, drainage, profile, UV and overlap checks run.

A typed `VisibilityCutFamilyError` can also exhaust the existing kernel and
filled-source charts for a source-certified exterior without holes. The late
`source_reflex_atlas.py` family combines the existing original-road normal
cuts with the existing bisectors of actual outer reflex curbs. All endpoints
remain on the same physical rim; these extra computational cuts neither
connect roads nor move a profile. The original 128-cut and 50,000-quad budgets,
fresh complete source/prefix proofs, exact sampled boundary, and physical and
stored 8:1/0.03 quality checks still qualify the result.

After the first recognized typed primary chart budget or grading rejection,
a separate derived-exterior
sampling retry can remove an exposed generated contact phase of at most 5 mm at a
genuine near-through circular contact. The original circle, radius, controls
and source approaches remain unchanged. Facing endpoints must be original
zero outer curbs or equal original half-crown plateau peaks; plateau height
and zero derivatives are preserved. Negative auxiliary leads require the exact
stored original root and contact knots. Positive contact phases additionally
require the identical original Route, Section, profile, directed end and curb
side on the exposed finite source-prefix interval. Original control knots in
that interval are excluded using metric station, including reversed Smooth
paths; a merely short curb edge cannot trigger the retry.
The projected source-union polyline adds
at most one micrometre to the recorded original chord error. This is measured
polyline approximation, not a claim of analytical equality. The retry runs
fresh continuous 2 mm proofs for every consumed source and incoming prefix
before normalization, after normalization, and after the normal chart build.
It preserves complete physical profile XY/Z and every genuine hole loop
byte-exact, and cannot spend its projection allowance twice.
Successful primary charts return their existing geometry. A temporary attempt
budget spans nested and piecewise charts of that source pool, allowing at most
one normalized candidate per construction. The original complete-width proof
is outside candidate-decline handling and remains fatal. After that proof
passes, only a typed candidate chart rejection or candidate-only strip clipping
can decline back to the original secondary families. Unrelated errors still
propagate; no persistent geometry flag marks an original source as valid.

The consolidated exterior footprint is a boundary/curve carrier, not an area
operand. Its `Region.boundary_only` flag makes `RegionUnion` reject accidental
use as pavement area; `CompoundGroup.area_footprints` retains the unchanged
original physical operands. The actual domain remains the exterior minus its
original holes. Thus neither source inventories nor future union consumers
can silently fill a separator merely by reading the derived exterior disk.
All normal stored-coordinate, height, drainage, profile, UV and overlap guards
still qualify the resulting Product.

Publication slope projection can disturb an entrance which passed the earlier
physical join check. Only the typed `PublicationEntranceSlopeError` enables
one directional reconciliation on that fixed chart. Its affine height slabs
use the exact profile-edge midpoint Jacobians of the authoritative entrance
check, with a 0.95 search target and the unchanged 3 percent acceptance limit.
The established corner-slab helper keeps its original default. The publication
projection accepts optional slabs only for this explicit retry; its ordinary
arithmetic is unchanged. Fixed source XYZ, XY, topology, crown bounds, global
physical and stored slopes, drainage, and stored entrance joins must all pass
before that retry returns. Disabling height projection also disables it.
The publication search also retains its established first feasible target
above the immutable boundary-gradient floor. If stored coordinates still fail,
successive targets use less of that same finite headroom instead of repeating
an identical failed target. This remains inside the original eight rounds;
targets stay above the fixed-gradient floor, and every actual published slope,
drainage, entrance and exact-boundary check retains its original limit.

`InteriorQuadRefinement` is a bounded alternative for a positive planar atlas
whose fixed-rim correction has not qualified. It splits long **interior** edges
with the existing direct 1:3 quad templates, never subdivides a physical rim
edge, and retains every original vertex with an explicit old-to-new mapping.
The caller must remap boundary references through that mapping. At most three
passes, 512 candidate faces per pass, and the ordinary 50,000-cell allocation
are allowed. The ordinary edge-ratio and corner-sine limits still qualify the
result; immovable bad boundary cells remain unsupported. This planar operation
does not replace the grass, height, drainage, UV, stored-float32, or global
overlap checks. The actual expanded Beijing case contains seven original
junctions; its 24,526-quad planar atlas is a compound cost, not a single-junction
mesh budget.
