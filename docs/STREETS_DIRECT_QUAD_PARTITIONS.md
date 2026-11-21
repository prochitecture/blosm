# Direct quad partitions for difficult physical junction regions

The Python street backend has two additional **planar chart policies** for
physical junction footprints that cannot be graded satisfactorily by one large
radial chart. They change only the computational partition inside the already
prepared pavement region. Source Sections, their incidence, incoming profile
rows, analytical curb curves, height limits and publication checks are retained.

Implementation:

- `way/street_network/python_backend/longitudinal_cuts.py`:
  elongated corridor/fork cuts and compatible four-sided boundary parity.
- `way/street_network/python_backend/reflex_cuts.py`:
  inward cuts from substantial concave curb runs, including rounded return caps.
- `way/street_network/python_backend/physical_hole_cuts.py`:
  bounded cuts through a certified physical exterior with preserved holes.
- `way/street_network/python_backend/partition.py`:
  conforming graph, direct chart assembly and bounded fallback selection.
- `tests/test_python_street_longitudinal_charts.py`:
  curved/rotated forks, compact return caps, fixed profile rows, cut parity and
  explicit node/cut budgets.

## Selection and unchanged successful cases

An elongated region or a substantial rounded return cap receives the existing
radial construction with its local repair first, but without final edge
rotations or global repair. If it already meets the existing edge-ratio limit,
that result is retained. Otherwise, an eligible direct partition is attempted
before those more expensive repair stages. The ordinary construction and its
correction remain available if the new partition cannot qualify. Other regions
keep their normal path.

A previously accepted mesh that needed final/global repair can consequently use
a different interior quad layout. The physical boundary, incoming profile rows
and surface constraints remain authoritative; this branch does not promise
identical interior vertex arrays for every previously accepted region.

A direct partition is also attempted after an ordinary leaf chart fails its
quality contract. Failed candidates do not become publishable merely because
they have a different topology. Numerical uncertainty, cut-budget exhaustion,
invalid cuts and insufficient quad quality remain explicit failures.

## Longitudinal chart policy

For one exterior physical loop, let the centered boundary coordinates be rows
of the matrix $X$. The dominant eigenvector of $X^T X$ provides a computational
longitudinal direction $a$; the perpendicular direction is $n=(-a_y,a_x)$.
This is a chart orientation, not a replacement street centerline.

Projected extents are

$$
L=\max_i(X_i\cdot a)-\min_i(X_i\cdot a),\qquad
W=\max_i(X_i\cdot n)-\min_i(X_i\cdot n).
$$

The branch requires $L\ge\max(6W,8s)$, where $s$ is the requested mesh spacing.
It divides the longitudinal extent into approximately $L/\max(W,4s)$ bands,
with an explicit maximum of 64 cuts. Each transverse scan pairs intersections
inside pavement. Contacts on incoming profiles and nearly tangent curb
contacts are skipped. Candidate endpoints are snapped to existing curb
samples, avoiding tiny residual physical edges caused by arbitrary cut
fractions.

A substantial concave nose at a fork is split locally by two inward diagonals.
Their endpoints remain on the actual sampled physical boundary. They are
subject to exactly the same cut-intersection and area checks as every other
computational cut. They do not change the corner radius or shorten any source
Section.

## Reflex-curb chart policy

A compact physical region may have several substantial concave curb runs and
no useful single radial center. With CCW boundary unit tangents $t_i$, the
signed discrete turn is

$$
\theta_i=\operatorname{atan2}
\left(\operatorname{cross}(t_{i-1},t_i),t_{i-1}\cdot t_i\right).
$$

Contiguous negative-turn samples form a concave run. Runs with accumulated
turn at least 1.3 radians provide one candidate cut. The sample at half of the
accumulated turn supplies a tangent, and its left normal supplies an inward
ray. The first physical boundary hit determines the opposite endpoint.

The endpoint is snapped to an existing vertex. When it lies on an incoming
profile, the selected existing profile sample also preserves compatible
interval parity; the profile is never resampled or moved. The policy permits
at most 16 cuts and requires a single certified physical exterior. Logical
members remain attached to that exterior; the policy does not infer new
source connections or merge unrelated roads.

## Quad parity without triangles

Certified replacement layouts with physical holes use an explicit separate cut
policy. Two scan stations through each hole open its surrounding pavement into
simple blocks. Intersections are paired across the complete oriented boundary,
so a cut crosses pavement only; every physical hole remains excluded. Four
bounded scan orientations allow the cuts to reach owned curbs while preserving
all prescribed incoming-profile samples. Each resulting block receives the
existing visibility decomposition and direct four-sided chart qualification.
The default limits are eight holes, 128 cuts and 3,072 visibility candidates per
block. Failure to open a hole or qualify the resulting charts remains an error.

Cut endpoints retain exact references to physical boundary chords. Additional
curb samples subdivide those chords and do not move their endpoints or infer a
new hole. Ordinary source-ownership cuts around separator holes and successful
single-loop layouts retain their established policies. The preparation stage
must first supply and certify the outer and inner physical curves; this chart
policy does not convert a source-graph cycle into a hole by itself.

Circular contact fitting keeps its established $10^{-8}$-metre acceptance but
continues the same bounded Newton solve toward $10^{-11}$ metres, with a capped
world-coordinate precision floor. A stall retains the best already accepted
root. This prevents a few nanometres of contact residual from leaving a longer
false segment where a sampled arc meets its nearly tangent source curb; it
changes neither the requested radius nor the boundary welding distance.

A certified complete paired corridor can supply transverse cuts from one
ordered original bank. Using both nearly parallel banks can create tiny
computational slivers. All original source claims remain attached to the
region. If that atlas cannot qualify, the principal-frame atlas is one bounded
alternative over exactly the same physical polygon and source profiles.

After those alternatives fail, a source-cycle or source-region aggregate with
a complete original-width certificate can combine original corridor normals
with two connections through each real hole. Visibility decomposition follows
the combined cuts. The existing hole, cut, node and 50,000-cell budgets remain.

A cached four-sided quad pattern needs compatible opposite-side interval
parities. Cut endpoints and substantial physical turns identify structural
sides. An odd side on an owned curb may receive one midpoint on its longest
existing chord. This leaves the complete physical polygon exactly unchanged;
all original vertices remain. Incoming profile rows never receive a new
vertex.

The resulting physical side counts are even, allowing shared computational
seams to use even counts as well. The graph is rebuilt with exact directed
boundary-span references, so height and curb provenance survive the extra
sample. The declared graph-node budget is checked again after this refinement.

Each block first uses a direct four-sided atlas. Small leftover blocks can use
the existing bounded radial quad pattern instead. Neither branch creates a
triangular intermediate mesh or merges triangles into quads.

When assembled shape correction leaves a positive chart with poor interior
grading, at most three anisotropic refinement passes may split long interior
edges into thirds. Every original physical boundary vertex and edge stays
fixed. The alternative is bounded to 512 candidate faces per pass and the
existing 50,000-cell limit, then must satisfy the unchanged pavement quality
limits. Exhausting the alternative remains a construction error.

Up to 64 remaining positive cells can instead receive local direct-quad
collars before anisotropic refinement. Every earlier vertex stays byte-exact;
only new interior vertices and faces are added. Each repaired cell and the
whole assembly must pass the same physical and actual float32 limits. This
separates incompatible fixed boundary-edge lengths without changing their
source curves or incoming profile rows.

After every assembled topology repair, the output is reindexed so that every
earlier vertex keeps its original integer ID. Boundary references, shared rim
seams and incoming source-profile aliases therefore retain their exact lookup
meaning. The final product checks the original analytical profile samples,
profile peaks and physical boundary heights again.

A certified source union can also fill the interior of its incidence cycle,
leaving one exterior and no hole. If its ordinary atlas and bounded assembled
repair fail, an explicit filled-source policy can use normals of the original
consumed corridors. A fresh continuous full-width source certificate is
mandatory. The policy retains the existing physical rim and claims, permits
bounded fixed-rim direct collars in its leaves, and applies the same assembled
quality, height, drainage and entrance-join checks. It creates no hole or source
connection. Successful ordinary single-loop charts keep their existing path.

For a physical hole, a later bounded cut-family alternative can end
computational cuts at existing physical rim vertices. This changes only the
computational cut endpoints, avoiding artificial tiny separator chords while
retaining every earlier curb and incoming profile coordinate.

Thin physical grass holes can have unequal sampling on their two banks. A
bounded monotone pairing retains every original bank point and adds only
collinear samples on existing chords. Capped fibre strips produce the disk;
the pavement receives exactly those same samples through conforming direct
quad templates. Odd face parity propagates through internal edges, leaving
incoming profiles unchanged. Local graded collars separate adjacent immutable
boundary edges whose lengths cannot fit in one accepted quad. The final
shared seams must match byte for byte, including cyclic order.

## Qualification and height construction

These charts do not require every block to have a single visibility kernel.
That requirement belongs to the radial mapping, not to a general positive
four-sided map. Instead, each direct leaf must produce positive quads with
maximum edge ratio at most 8. Four-sided atlas qualification also requires a
minimum corner sine of 0.03. The assembly verifies shared seams, manifold edge
incidence, exterior boundary count and Euler characteristic. The partition
also verifies exact pavement area and refuses cuts through physical holes or
outside the pavement.

Grass retains its existing separate contract: strictly positive corner
Jacobians and edge ratio at most 8 in both physical and stored float32
coordinates. A sharply tapered fixed grass rim can force small boundary
angles; its minimum corner sine is recorded as metadata. This does not relax
the pavement atlas's 0.03 requirement. Disk qualification also checks its exact
rim, area, Euler characteristic, manifold incidence and transported UV data;
height construction must still preserve its plateau and zero grass slope.

The assembled chart still passes through `CompoundSurfaceBuilder`, including
its boundary conditions, incoming-profile joins, height ceiling, longitudinal
slope checks, drainage checks and stored-coordinate publication validation.
An acceptable XY chart alone is not an accepted street surface.

An otherwise qualified ordinary compound can have more than two exposed
ownership changes for a shared source Section. Only when the established cut
policy rejects that ambiguity, all incoming profiles remain complete, and a
continuous full-width source certificate passes, computational visibility cuts
may replace the failed ownership cuts. The physical rim, logical claims and
member junctions remain unchanged. The chart records the source-width proof
and still requires drainage and entrance-join reconciliation.

## Focused real-data evidence, 14–15 November 2026

These are focused numerical fixtures, not claims that the entire city import
has qualified. Development reports are under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114/`.

| Fixture | Result | Junction quads | Maximum edge ratio |
| --- | --- | ---: | ---: |
| Bergamo site 230, actual style-radius policy, Sections 114/115/116 | Complete focused preparation passes | 4,133 | 5.462 |
| Beijing site 15, Sections 6/7/8 | Complete focused preparation passes | 3,570 | 6.493 |
| Beijing owned-prefix compound, sites 794/802/810 | Compound planar and height surface passes | 2,314 | 7.622 |
| Beijing source cycle, sites 185/193/194/253/261/291/1211 | Compound planar, height and product checks pass, one physical hole | 25,017 including 623 grass | 7.748 |
| Bergamo source cycle, sites 903/1163/1166/1249 | Compound planar and height surface passes, two physical holes | 15,559 including 510 grass | 7.894 |
| Bergamo source cycle, sites 995/1182/1193/1198/1477/1478/1479 | Compound planar, height and product checks pass, three physical holes | 23,867 including 7,108 grass | 7.980 |
| Bergamo site 43, existing longitudinal cut family | Normal compound product passes after the bounded standard correction retry | 16,370 | 7.975 |
| Bergamo normalized source aggregate, sites 447/449/507/574/589/771 | Normal compound product passes, one physical hole | 4,041 including 38 grass | 7.649 |
| Bergamo physical compound, sites 510/1003 | Normal compound product passes after bounded interior topology repair | 8,371 | 7.749 |
| Beijing local merged mouth, sites 154/525 | Local compound product passes; full Bundle plan remains separate | 1,147 | 7.738 |
| Bergamo filled source cycle, sites 1238/1405/2457 | Normal compound product passes with original corridor cuts, no hole | 2,995 | 7.880 |

The owned-prefix fixture's complete reduced-network run subsequently reaches a Bundle
chain-scope failure; it is not evidence of full Bundle or full-network success.
The compound surface itself clears the earlier chart failure. Its focused
planar and height calculation took about 0.74 seconds in the development run;
timings are diagnostic and depend on concurrent system load.

The seven-site Beijing calculation retained its original 32.965-square-metre
physical hole and qualified the complete height/drainage/entrance checks in
about 106 seconds on the latest source preparation. This focused construction time is not evidence of an
interactive preview latency target or a complete city-import pass.

The Bergamo paired corridor retained its two 0.987- and 6.294-square-metre
grass holes. Its fixed-rim fibre disks share 192 added collinear samples with
the pavement, preserving every earlier rim coordinate. The complete focused
surface took about 181 seconds. Independent final checks measured stored
pavement ratio 7.967 and minimum corner sine 0.0388, stored grass ratios 7.722
and 7.658, zero interior minima and a maximum entrance longitudinal slope of
1.05 percent. All final shared XYZ seams and original boundary heights were
exact; projected pavement and grass interiors had no overlaps.

The seven-site Bergamo pool retained its three 116.810-, 40.258- and
10.524-square-metre grass holes. Its complete focused surface took about
77 seconds. Independent final checks measured stored pavement ratio 7.979 and
minimum corner sine 0.0412, maximum physical slope 0.0903, no interior minima,
and a maximum entrance longitudinal slope of 0.00690. Exact shared XYZ seams,
original boundary heights, UV maps and pavement/grass overlap checks passed.

The partition regressions include original compound source controls embedded in
the repository. They exercise physical preparation, the complete compound height
and publication guards, and stored-coordinate quality after a rigid transform;
they do not read a generated artifact or qualify the omitted remote network.

The ordinary partition, quad-atlas, radial-chart, quad-rotation and analytical
sampling suites run together with the source-rim and bounded-repair regressions.
Regenerate numerical evidence from the current source when qualifying a new
release; local section IDs are diagnostic labels, not hardcoded branch choices.

## Ordinary forks after a radial budget failure

A valid ordinary junction can exhaust the radial chart's cell budget before its
later grading fallback is reached. `QuadChartCellBudgetError` distinguishes that
computational failure from an invalid physical boundary. The chart then tries
the existing bounded reflex-curb partition, retaining every prescribed boundary
sample. This also permits ordinary forks below the cheap reflex preselection
threshold to try that construction. Successful existing charts stay first.

The translated acute half-profile regression uses that path to construct 5,480
quads with the original circular curbs and complete profile rows. Its complete
surface passes stored-coordinate ratio eight, the original crown envelope,
zero interior minima and the three-percent incoming slope limit. The fallback
always runs coupled height/drainage reconciliation before publication.

## Independent lead and arc sampling

Long acute forks can contain a short tight circular arc between much longer
source-curb leads. Applying the worst arc curvature to the complete curve can
force unnecessary rows on those leads. After established bounded charts reject
an ordinary single-junction exterior, `piecewise_sampling` retries the original
analytical curve using separate certificates for its first lead, arc and last
lead. It changes neither source controls nor the analytical fillet.

For each piece with curvature bound $k_i$, parameter speed bound $v_i$, spacing
limit $h$ and unchanged curve-error budget $\varepsilon$, the parameter step is
bounded by

$$
\Delta s_i \leq \min\left(
\frac{h}{v_i},\;
\frac{\sqrt{8\varepsilon/k_i}}{v_i},\;
\frac{\pi}{2k_i v_i}
\right).
$$

The curvature terms are omitted for a straight piece. Lead/arc contacts,
retained ownership contacts, Exact source knots and all incoming profile rows
remain compulsory. Allocation and parity checks remain bounded. Metadata records
the per-piece curvature, speed, step and chord-error certificates.

For a Smooth source lead, the retry uses the existing Bernstein curvature proof
restricted to the original metric interval that the lead actually traverses.
A tight bend elsewhere on the same source Path therefore does not force dense
samples onto the retained lead. Original Path controls, finite support, offset
fold checks and chord-error budgets stay unchanged. Exact source knots remain
compulsory, including reversed traversal. The legacy sampler still uses its
original global bound and remains the first construction.

The default sampler remains first, and the retry cannot recurse. It applies only
to one ordinary footprint without consumed claims or holes. It repeats the same
mesh, stored-coordinate, elevation, drainage, incoming slope and overlap checks.
The focused long-fork surface passes with 8,455 quads, stored ratio 7.629,
minimum corner sine 0.140, zero interior minima and maximum incoming
longitudinal slope 0.00424. Its physical slope 0.17459 is below the independently
derived incoming-profile allowance 0.176; that allowance is unchanged. Complete
component and native-import qualification remain separate.

## A failed radial center in a prepared compound

`QuadChartCoreError` identifies a radial construction whose final six vertices
cannot form valid positive quads. Reducing the positive-corner threshold would
accept very small or badly shaped cells. Instead, after the established radial
and reflex charts fail, a prepared multi-junction compound with one exterior may
try a bounded fixed-rim reflex atlas. Positive leaf charts may be retained as
intermediate construction data; the complete assembled chart must still pass
the original eight-to-one edge ratio, minimum corner sine 0.03, positive-corner
and stored-coordinate checks. Local collars are direct quad constructions.

No new physical region or source ownership is inferred. The prescribed outer
boundary, complete incoming profiles, source claims and finite cell budget stay
authoritative. The height builder then applies the same elevation envelope,
coupled drainage correction and incoming longitudinal slope limits. A focused
three-junction regression qualifies 4,609 quads and reproduces the original
radial-core failure from source controls embedded in the test repository.

## Deferred drainage reconciliation

A slope projection may lower a few interior vertices enough to create a small
basin. The builder rejects that candidate. Existing refinement attempts remain
first. If they also fail, the builder can reconcile drainage and slopes on the
saved original chart using the existing bounded spill-floor solver. Every
physical boundary XYZ remains fixed. The failed refined meshes are discarded;
they are not published or used as new source geometry.

`deferred_projection_reconciliation` records the rejected initial projection and
discarded refinement attempts. The focused three-arm regression qualifies 3,211
quads, with no interior depression, a maximum incoming longitudinal slope of
0.00414 and the unchanged incoming-profile physical slope allowance. Disabling
height projection does not enable this recovery.

## Several fork noses in a moderate-length junction

A footprint with several almost parallel entrances may contain more than one
substantial concave curb run. The original longitudinal family selects the
largest nose and normally requires a length-to-width ratio of six. Those
defaults remain unchanged. After its established charts, corner seed and
piecewise sampling have failed, a separate bounded recovery may try a ratio of
two and place the existing verified nose cuts at every substantial concave run.
All cuts remain computational: source controls, round fillets and complete
physical profile rows are retained.

Every block and the complete assembled mesh must pass the existing positive
quad, edge-ratio, angle, float32 publication and finite-budget guards. Heights
then pass coupled drainage and incoming-slope reconciliation. The source-derived
four-arm regression qualifies 3,743 quads, with stored edge ratio 7.862, minimum
corner sine 0.109, no interior minima and maximum incoming longitudinal slope
0.01085. It reproduces the difficult geometry without reading an OSM file or
saved scene. The production construction has no extract-name or source-ID test.

## Source-normal recovery after all physical-hole atlases fail

A compound's source-normal atlas can be useful even when its older preparation
did not record a full-width proof. `physical_hole_recovery` does not infer that
proof from a preparation label. After all existing hole atlas families fail, it
certifies every consumed finite source strip and explicit hole exclusion under
the unchanged two-millimetre combined budget. It then tries the existing source
normal cuts using existing boundary vertices, and repeats the complete source
proof on the resulting sampled exterior.

Both pavement and grass must meet positive-corner and physical/stored edge-ratio
eight checks. Pavement also retains corner sine 0.03. Grass retains its existing
positive-corner policy: a genuine narrow separator tip can have an angle below
that pavement threshold, and preserving its two physical banks makes such an
angle unavoidable. Actual grass angles are reported. The source proof is carried as
`source_width_visibility`, so ordinary coupled height, drainage and entrance
checks remain mandatory. This changes computational cuts and sampling only;
the original real holes, physical mouths, finite source paths and cell budget
remain authoritative.

If the final source-normal family itself exhausts planar grading, an explicit
last energy family may relax that rejected pavement candidate for at most 500
steps with powers one and two. It then reuses the general chart's finite
300-step correction relaxation, including its existing bounded direct-quad
topology repairs. The cheaper 100-step source child remains first and unchanged.
This extra effort is scoped to the final source-normal recovery; the global
chart correction is not changed. Every original physical rim coordinate and
source vertex ID stays exact. Complete-source certification runs before and
after, total pavement plus grass remains within the same cell budget, and the
normal physical/stored, height, drainage, profile and seam checks still apply.

The focused actual six-site Bergamo compound qualifies 10,541 pavement and
1,903 grass quads, retaining its 24.5078 square metre hole. Independent checks
pass stored pavement ratio 7.7498 and sine 0.05020, stored grass ratio 7.7322,
all 85 profile aliases and five peaks, exact heights and shared seams, finite
UVs, zero interior minima, maximum join 0.00379 and no overlaps. All eight
consumed source strips retain the unchanged 2 mm certificate over 1,488 cells.

## A source owner seam is not always a visibility partition

A valid ordinary compound can have a source-ownership cut whose remaining
block has no visibility kernel. `VisibilityKernelError` identifies this
construction failure without changing the kernel arithmetic, tolerance or
messages. Existing longitudinal, reflex, fixed-rim and all-nose alternatives
remain first. If those reject a complete ordinary multi-junction compound,
`visibility_recovery` can use the existing physical-rim visibility cuts.

This last policy requires complete original external profiles and fresh
continuous full-width certification for every consumed source Section before
and after sampling. The sampled physical coordinates and ownership spans must
match byte-for-byte; no replacement physical layout or source topology is
created. Final physical and stored pavement ratio/sine checks remain 8/0.03.
The source proof activates ordinary coupled height, drainage and join checks.

The two-junction source-control regression reproduces the missing kernel with
seven original Smooth Paths and two real compatible Connections. The actual
Bergamo compound qualifies 2,439 normal Product quads, stored ratio 7.7307 and
sine 0.07199, all 68 profile aliases and four peaks, exact boundary heights,
zero interior minima, maximum join 0.00694, finite UVs, no overlaps and the
unchanged 2 mm source certificate over 34 cells. Its rotated regression and
clipped-source rejection also pass.

## Refined height charts keep their original shape limits

Height refinement can leave a positive fixed-topology chart in a poor shape
solver basin. The first `QuadShapeCorrection` attempt remains unchanged.
Only after it fails, `refined_grading_recovery` tries the existing fixed-rim
energy relaxation with powers one and two, at most 500 steps each, followed by
the unchanged shape solver. No vertex, edge or cell is added. Original boundary
XYZ, seam indices and quad topology must remain exact; physical/stored pavement
ratio 8 and sine 0.03 still apply. The original analytical height field is then
sampled and the existing coupled slope, drainage, join and publication checks
run normally.

The five-arm Smooth source-control regression reproduces the stalled refined
chart. The actual normal Product retains 4,320 quads, all 85 source-profile
aliases and five peaks. Its stored ratio is 7.8563, sine 0.05012, maximum join
0.00439 and interior minimum count zero. Independent checks also pass exact
boundary heights, physical/stored slope limits, finite UVs and overlap. An
already qualified refined chart never calls the new recovery.

## Banded grass recovery after older physical atlases fail

A long separator with rounded ends can have millimetre longitudinal station
gaps beside much longer physical bank chords. Extruding each such station as
a complete fibre row can fold the tip transition. `hole_banded_chart` keeps
the original physical disk and divides it into bounded computational bands,
with cut endpoints at existing bank vertices. All original rim coordinates
remain byte-exact and in order. Any additional rim points lie monotonically
on existing chords; normal conforming pavement refinement shares those points.
No original curb, source path, mouth, real hole or source claim is replaced.

This recovery is deferred until every existing physical atlas family and the
source-normal recovery have failed. A typed grass grading failure identifies
one family to retry with `recover_grass=True`; ordinary successful families
retain precedence. At most 64 longitudinal bands and the existing compound
cell budget apply. The final grass disk must be a conforming direct-quad disk
with positive corners and physical/stored edge ratio at most 8. Grass tip
sines are reported under the existing grass contract. Pavement still requires
physical/stored corner sine at least 0.03 after the shared-rim refinement.
Source-width certification and coupled height, drainage and incoming-join
checks remain mandatory.

The actual Beijing four-site source compound retains all 610 prepared grass
rim vertices and its 402.6189 square metre hole. The focused normal Surface and
Product qualify 11,011 pavement and 2,363 grass quads after 28 shared collinear
rim insertions. Independent checks pass pavement stored ratio 7.8611 and sine
0.05305, grass stored ratio 7.8578, exact physical and stored seams, zero
interior minima, maximum entrance join 0.01365, finite UVs, no overlaps and
the unchanged 2 mm complete-width certificate over 2,648 original source cells.

## Paired terminal intervals after all older atlases fail

Normals from both opposed members of a partially owned Bundle can create poor
computational blocks. `paired_interval_atlas` keeps both original physical
intervals and every source claim, but tries the normals from one member at a
time, together with the existing physical-hole bridges and unpaired sources.
It rechecks the original terminal ownership and complete half-profile cuts
through `partial_bundle_claims.prepare`. Fresh 2 mm width certificates include
every consumed strip and finite external prefix before and after charting.

The family runs only after all earlier physical atlases, source-normal recovery
and any older grass recovery fail. Each of its two member choices uses the
unchanged 100-step child first, then the existing failed-only 500-step energy
seeds and 300-step general correction budget. Combined pavement/grass allocation,
physical and stored pavement ratio 8/sine 0.03, grass positivity/ratio 8, original
rim and profile identity, coupled height/drainage/joins and publication checks
remain unchanged. Successful older charts keep precedence.

A fresh seven-site original-source compound with three real holes qualifies
through the normal Product provider with 18,828 pavement and 6,342 grass quads.
Its current source plan advances a complete original road mouth past a consumed
crossing strip; the former cached surface belongs to the shorter source plan.
The new Product keeps all 68 original profile aliases and four peaks. Stored
pavement ratio is 7.9813 and sine 0.04285, physical slope is at most 0.09, entrance
join is 0.00669, and interior minima are zero. Independent source-width, exact
height/seam, UV and overlap audits pass. The focused normal Product takes about
419 seconds; this result does not establish interactive or whole-city readiness.

## Exact-rim standalone separator bands

`fixed_rim_bands.BandedSeparatorChart` is a final separator fallback after the
older layered, four-sided, partitioned and angular families fail. Up to 64
transverse computational bands end at existing physical bank vertices. The
partitioner certifies each complete cut and local polygon; all original rim
samples remain in their original order. This standalone family cannot add
collinear samples unilaterally to an adjacent carriageway's published edge.

Local charts first try the established builders. The failed-only grass builders
then use bounded fixed-rim relaxation and shape correction. If a four-sided
parity termination leaves three fixed collinear vertices in one cell,
`positive_quad_patches` replaces a two-cell computational patch with two
rotated quads or three direct quads around one new interior point. The twelve
new corner Jacobians are affine inequalities in that point; their half-plane
intersection supplies a certified positive seed. No triangles, physical rim
vertices, profile rows or source ownership are introduced or changed. The
same complete boundary edges and manifold incidence must survive, including
float32 storage. Shape, height, UV, allocation and overlap acceptance remain
separate mandatory checks after this seed repair.

The complete grass disk must retain positive physical/stored corners and edge
ratio at most 8, fit its original cell budget (at most 50,000 here), have Euler
characteristic one, and have exactly the original exterior edges. The ordinary
plateau extension supplies heights. The distinct pavement sine limit 0.03 is
unchanged; this family constructs grass only.

A focused normal separator construction over 6,756 prescribed boundary vertices
qualifies 21,432 grass quads without an added rim sample. Stored edge ratio is
7.85360, stored minimum corner sine is positive (0.00008166), and its original
plateau is flat. Exact boundary XYZ, conforming disk topology, UVs, material and
59,620 independent overlap candidate checks pass. Fifteen focused patch,
separator, band and existing seed-repair regressions also pass. This is a
separator qualification, not a fresh complete-city import or a performance
target for every Bundle.

The exact physical-boundary rules also apply to the separate
[gapped original-source body family](STREETS_GAPPED_SOURCE_BODIES.md), which
retains real one-sided junctions between independently meshed body cores.
