# Beijing and Bergamo import qualification

Work started 14 November 2026 in the isolated `streets_development` checkout.
Reference `buildings_gn` was already synchronized at
`386abdb15808ce397149b6c2c951b6bb65e6c87f`.
On 20 November the isolated development branch fast-forwarded to
`4520f85b8112be126a123912233e19aa3b341b41`.

## Objective

Correct the two real OSM imports, including relevant confirmed historical issue
classes, without weakening mesh validation or substituting hand-edited examples.
The starting evidence is the [GitHub audit](STREETS_GITHUB_ISSUE_AUDIT.md) and
[large-area import assessment](STREETS_PROGRESSIVE_LARGE_AREA_TEST.md).

## Work and acceptance

1. Preserve OSM topology at XY crossings: only original node incidence supplies
   road connectivity, with explicit treatment of grade transitions.
2. Support returning/closed Sections while retaining both endpoint incidences
   and avoiding duplicate numerical resource dependencies.
3. Define local, explicit boundaries at unsupported elevated transport features;
   a missing approach must not silently become a normal junction or suppress an
   entire supported ground network.
4. Reassess Bundle membership and rounded intersection construction after the
   upstream fixes, using captured source identities and focused reproductions.
5. Qualify original Beijing and Bergamo imports, existing street regressions,
   progressive/blocking parity and publication cleanup. Report supported coverage
   separately from deliberately excluded features. Inspect generated surfaces.

Rounded corners, carriageway profiles, source alignment, quad quality, UVs and
material behavior remain requirements. Production fixes must apply by geometry
and topology, never by extract name or OSM way ID. Existing validation failures
remain failures until their underlying cause is corrected.

### Verification on 20 November

The full 688-Section Bergamo run on revision 120 saved 19 completed
junction groups before two near-straight seams on the same Street failed
during construction of the group containing junctions 2998, 3001 and 3015.
The two style fillets had a 2 m radius. Their deviations from the endpoint
chords were 6.55 nm and 0.13 nm respectively. One auxiliary polygon fell
below the minimum area for a valid region and the other could not be
distinguished from point contact during the polygon union.

Revision 121 recognizes this case at the original shared Section endpoint.
It requires the same Street, matching oriented carriageway width, offset
and crown, and a successful fit at the specified style radius. Both bank
runs and the fitted arc must fit within the existing curve sampling
tolerance and a 1 mm diameter. A separate transverse check uses the stricter
source-bank sampling tolerance. Explicit radii still use the fillet builder.
On the original three-junction group, both boundary loops survived and
the pavement passed with 5,118 quads, a maximum edge ratio of 7.75, and
a 16.50% maximum physical slope under its 17.60% profile limit. UVs were
finite and 15,346 candidate overlap pairs passed. See
`nov20_688_exterior119_1/report.json` and `nov20_688_handoffs_mesh.json`.

That full run also recorded a 124.9 s numerical advance during visible-bank
lookup. The updated lookup visits only source edges overlapping the queried
interval, validates each original curve once per uninterrupted walk, and
reuses the visible interval across trial radii. All 612 endpoint queries
on 15 real Bergamo curves gave identical source witnesses, including at
one-ulp perturbations. In that isolated comparison, exhaustive lookup took
3.46 s and the indexed queries took 0.034 s. These are query timings.
The change in total import time still needs measurement.
See `nov20_real_bank_queries121.json`.

The next full 688-Section run passed the group at 2998, 3001 and 3015.
It subsequently stopped at a separate three-way junction, 1609. At the root
of the 7 m-wide Section 804, the proposed pavement outline exposes 1.737 m
of an internal transverse edge between two curbs. The source-containment
check rejects that outline. See `nov20_688_full121_1/report.json` and
`nov20_688_ordinary_seam.json`.

At 1609 the original 2 m circle between Sections 805 and 804 touched the
backward extension of the narrow bank, 7.224 m behind its endpoint.
Revision 122 first tries to connect the two finite banks when an ordinary
width check fails in this configuration. It anchors the transition at the
root of the wide bank and searches along the narrow bank, within its original
entrance. On the full component's original input the narrow-bank contact is
3.9375 m from the root. The maximum certified curvature is 0.498 per metre
against a limit of 0.5. The two other corners and all profile planes keep
their original positions. Both passes over the complete source widths
succeeded. Both incoming profiles are half crowns in the complete Bundle
network. The new curb joins the narrow carriageway's 8 cm crest to the wide
carriageway's outer edge. Explicit corner radii and corners between two
separator-facing crests are excluded from this construction.

The revised surface builder produced 1,434 pavement quads from the saved
688-Section context. The maximum edge ratio was 7.634 and the steepest
stored face was 5.259% under the 10% limit. UVs were finite, and the overlap
check examined 4,845 pairs. See `nov20_1609_full_context.json`.

The full Beijing 1,072-Section run on revision 121 generated the previously
failing surfaces at 553 and 691, then stopped at the pavement shared by
junctions 784, 925 and 1797. Its computed quad chart exhausted the shape
correction budget. This stops the component before publication. See
`nov20_1072_full121_1/report.json`.

A replay taken from 23 original Beijing Sections reproduces this junction's
entire planar outline byte for byte, both owned Sections and all nine recorded
mouth stations. The primary mesh had 2,817 out of 10,061 quads above the 8:1
edge-ratio limit after the ordinary correction. Reducing the exterior from
650 to 482 samples preserved full source width but the radial chart still
exhausted its cell budget.

For a filled junction with owned Sections and short generated corner contacts,
the Python backend can repair the dense, positive mesh on its original rim
before trying an alternative boundary sampling. It groups movable
vertices so that two vertices updated together never share a face. A move
must reduce the total shape error on that vertex's incident faces and keep
all of them positive. Each vertex can shorten its own trial move. After 58
sweeps on the actual rejected mesh, ten cells still needed the existing local
quad repair, which added 40 quads. The finished planar mesh retained all 650
boundary vertices exactly. The maximum edge ratio was 7.750 in double
precision and 7.910 in Blender storage precision (rounded up). The stored
minimum corner sine was 0.0494 against the 0.03 limit. The overlap check
examined 50,095 pairs. The repair yielded control 233 times, with a longest
advance of 0.214 s in this diagnostic. These timings cover the repair only.
See `nov20_784_compact_incidence.json` and
`nov20_784_independent_vertex_grading.json`.

On this mesh, float32 storage raised the steepest face from 13.037% to
13.320%. A height correction using the existing safety margin brought the
stored slope below the limit, but its last iteration made a small depression.
The publication guard now also tests the filled surface when a bounded
projection reaches its iteration limit. The test raised the depression by
0.0206 mm, kept every fixed vertex and met the original source-height
constraints. The stored maximum was 11.297% against the 13.037% limit,
and the drainage check found zero residual fill in both precisions.
See `nov20_784_spill_trial1/stdout.log` and `nov20_784_replay_surface.json`.

The production ownership provider also built the entire pavement directly
from the 23-Section reproduction. Junctions 784, 925 and 1797 shared one
surface. It retained the 650 boundary vertices and produced 10,101 quads
with 40,404 finite face-corner UVs. Its stored maximum slope was 11.297%.
The provider yielded 2,650 times in 68.7 s. Its longest uninterrupted step
took 2.73 s. These measurements cover that junction group on the bounded test
process. See `nov20_784_normal_product.json`.

The full 1,255-Section Bergamo run emitted the surfaces for junctions 199
and 261. At junction 43 the profiler measured a single 128.0 s numerical
advance. The chart that fits the street ends was calling the synchronous
mesh builder inside a progressive operation. The progressive path now
forwards the child chart's checkpoints during construction and grading,
and while checking heights against the incoming streets. Radial quad
construction also yields between layers and during shape correction.
The 128.0 s measurement is from the earlier code. See
`nov20_1255_full121_3/step_profile.json`.

The 688-Section Bergamo run also recorded a 50.86 s uninterrupted step
while fitting exposed carriageway edges. The step included multiple polygon
unions. The join builder now forwards each union's checkpoints and the
checks of the original street widths. During a union, long polygon rims use
the same bounded spatial index as the width checks to find candidate edges.
The final ray and distance tests retain the original edge order for ties.

A comparison on 12 groups taken from the full Bergamo ownership run covered
311 input polygons and 68,823 vertices. All 25 resulting loops matched the
previous implementation in vertex order, coordinates, source-edge identity
and edge parameters. Together the unions took 17.59 s, down from 25.15 s.
These timings cover the polygon unions only. See
`nov20_688_union_parity.json` and `nov20_688_full122_1/step_profile.json`.

For one 4,444-edge original Bergamo strip, 128 queries gave exactly the
same ray classifications, edge distances and source intersections after
batching. Peak temporary allocations measured separately for the three
operations fell from 32.0, 41.1 and 40.0 MB to 2.7, 2.4 and 2.9 MB
respectively. These are allocation measurements for those operations,
not whole-import memory figures. See `nov20_688_region_scratch.json`.
The final batched version also reproduced all 25 loops on the same 12
groups. Comparison time summed over the groups was 31.66 s for the
reference and 20.93 s for this version. The final group ran separately
after the system memory reserve interrupted the test. See
`nov20_688_union_parity_final.json`.

In the Beijing ownership pass, coupled Bundle junctions now forward
construction checkpoints whenever discovery expands to an adjacent junction.
Groups of more than eight original junctions also forward checkpoints while
checking carriageway widths. A cancellation test interrupts the width check
and verifies that the previously published proof remains intact.

The 1,072-Section Beijing trace also identified long advances while testing
whether a candidate curb remains exposed and while repairing poorly shaped
quads. Both operations now pass their internal checkpoints to the progressive
importer. On the four-arm reproduction of junction 1683, the previous curb
calculation took 10.81 s in one advance. The revised calculation took
10.21 s with a longest advance of 0.031 s. All 310 patch vertices, source
edge labels, contact stations, 1,025 curve samples and the construction
report matched. This timing measures the numerical curb calculation.
See `nov20_original_root_callbacks.json`.

Rebuilding Beijing junction 68 from the complete 1,072-Section network
reproduced the saved 2,205 vertices, all 2,058 pavement quads, UVs and
material exactly. The original full run recorded a 13.96 s advance at this
junction. The rebuilt surface yielded 2,335 times, with a longest advance
of 0.69 s. Its stored slope was 16.59% under the 17.60% source limit,
and the overlap check examined 6,058 candidate pairs. See
`nov20_full_owner_68_parity.json`.

For two other real Bergamo approaches, 164 and 163, the displaced-contact
search now yields during each radius trial and during the iterative fit.
The previous search occupied one 9.81 s advance. The same calculation
took 9.84 s across 2,332 checkpoints, with a longest advance of 0.158 s.
Both runs produced identical contact stations, circle, source evidence and
patch coordinates. These times measure one corner search. See
`nov20_contact_callback_measured.json`.

At Beijing junction 856 the previous trace spent over 40 s in one unfinished
advance of a nested junction chart. With the chart yielding to its caller,
80 successive internal checkpoints were visible. The longest of those
advances took 1.76 s in the saved run. Cancelling there closed the chart, released its retry
state, and left the junction unpublished. This checks cancellation at that
specific stage. See `nov20_856_callback120.json`.

The complete 1,072-Section Beijing input passed ownership planning from an
empty cache on revision 116. Progressive numerical preparation also emitted
the pavement surface at owner 91, its alias at 99, and the five-junction
pavement with two islands at owner 122. It subsequently built junctions 354,
531, 540 and 547 and stopped at a separator at owner 553. That separator
and the grass slope at owner 691 now pass individual physical mesh checks.
The component report is `nov20_beijing_solo116/report.json`. An emitted
surface precedes the final network-wide overlap check and Blender publication.
The newly qualified [side-street junctions](#side-street-on-one-carriageway-at-beijing-junction-354)
and the [Bergamo acute corner](#acute-exterior-corner-across-two-junctions)
have separate local surface checks.

A later run on the same 1,072 Sections generated surfaces past junctions
784 and 856. Its next failure was at junction 912, where the exposed outer
edges of Sections 457 and 455 meet at 89.53° with different crossfalls.
A nine-Section reproduction gives the same boundary-condition failure.
The revised corner test admits a 2 m style fillet after checking that both
finite outer edges meet at the original tangent intersection. The revised
outline passes the ownership checks on that reproduction. The full component
failure before this correction is recorded in `nov20_1072_dense122_3/report.json`.

At junction 912, the first eight grading projections reduced the residual
spill depth from 1.62 mm to 4.06 µm. Revision 123 permits up to eight further
projections when the slowest of the last three measured reductions would
bring the residual below 1 µm within the remaining budget. Each projection
still measures the actual slope and drainage. In the separate Blender
precision check, a projected surface meeting both slope bounds can supply
spill floors for the next round when filling a depression directly would make
a face too steep. That path checks the drainage of the stored vertices before
returning.

The nine-Section reproduction now generates the complete junction product.
Its 9,279 pavement quads have 37,116 mapped UV corners. The stored maximum
slope is 12.944% under the 13.037% limit, and the stored drainage test finds
zero trapped depth. The largest edge ratio is 7.864. The local overlap check
examines 50,668 candidate pairs and passes. Rebuilding junction 912 with the
physical neighbors saved from the original 1,072-Section component produces
the same counts, slopes, drainage and overlap results. The maximum entrance
longitudinal slope measured during grading is 0.638%. See
`nov20_912_local_surface.json` and `nov20_912_full_surface.json`.
The nine-Section fixture isolates this junction. Junctions 215 and 219 on the
opposite carriageway belong to the city input and supply additional incidence
for building the complete Bundle. City-wide overlap and publication are
separate checks.

As a further diagnostic, the current surface builder generated the next 35
independent Beijing junctions through ID 1223 using ownership layouts saved
from revision 122. All 37,469 quads passed the mesh checks and had finite
UVs. This pass covered individual junctions. It omitted street strips and
the final overlap check between published parts. See
`nov20_beijing_late_surfaces123_from_913.json`.

An independent replay of Beijing junction 122 from the 1,072-Section
network reproduced the saved pavement and grass exactly. The coordinates,
quad indices, material roles and UVs matched for 23,111 pavement quads
and 13,540 grass quads. The stored maximum slope was 12.864% against
the 13.037% limit. The combined overlap check passed for 244,667 candidate
pairs. Its longest numerical advance was still 10.75 s in this replay.
See `nov20_full_owner_122_parity.json`.

The 688-Section Bergamo ownership run reached a separate width failure at
junction 1700. Section 846 is 7 m wide and meets Sections 850 and 851,
which are 3.5 m and 2 m wide. Both outer corners touch extensions of a bank
behind its original endpoint. The fitted contact stations behind the endpoint
are −12.75 m for Section 850 and −7.62 m for Section 846. The ordinary
outline clips Section 846 and the source-union fallback rejects an internal
construction seam. The local reproduction preserves five original Sections,
including the connected pieces of the incoming routes. The geometry and
failure record are in `nov20_688_next_seam_network.json` and
`nov20_688_full122_6/report.json`.

The two affected banks form the same inner crest on consecutive Sections
of the same Street and Bundle member. Revision 124 checks the two directed
endpoint incidences, the authored profiles, and their common crest height
before fitting a finite width transition. The contact on Section 850 moves
to 4.8125 m along its real bank. The transition retains the 8 cm crest and
zero transverse gradient. The original entrance planes and both other corners
keep their positions. On the input saved from the full component, all source
width checks passed. The pavement has 3,221 quads, a maximum edge ratio of
7.492 and a maximum stored slope of 17.508% under the inherited 17.6% limit.
Its UV and local overlap checks passed. See `nov20_1700_full_context.json`.
Complete Bundle and city publication still require separate qualification.

The opposite member of Bundle 8 has a second width change at junction 2753.
Section 1377 protrudes beyond the entrance of Section 1375. Moving that
entrance along its original centreline clears the side street, but the initial
layout still clips the wide carriageway. The previous width-transition retry
started again from the original entrance, losing the clearance it had found.
Revision 125 fits the transition on the cleared layout and certifies all
three incoming widths again on the final analytical sampling. At this junction
the entrance moves from 0.587 m to 1.058 m along Section 1375. The transition
reaches Section 1372 at 4.375 m, with a certified maximum curvature of
0.4964 per metre against the 0.5 limit. The accepted surface has 352 quads
and 1,408 UV corners. The maximum stored slope is 16.998% under the 17.6%
limit. All 1,001 candidate overlap pairs passed. See `nov20_2753_provider_result.json`.

Both junctions now supply the separator with a continuous carriageway edge
at the original 8 cm crest. The edge at junction 1700 contains 331 vertices,
and the one at junction 2753 contains 39. Every consecutive pair is an actual
pavement boundary edge. Requesting each chain in the opposite direction gives
the identical vertices in reverse order. See `nov20_bundle8_crest_chains.json`.
This checks the interface with the separator. The complete grass mesh and
its overlap with all surrounding roads still need a whole-network test.

On Bergamo junction 1386, the earlier progressive producer spent 103.92 s
in one advance during refinement. Replaying that original junction with
the same source network and ownership after forwarding the refinement and
height-field checkpoints produced 4,475 pavement quads. The longest advance
in the replay was 2.72 s across 3,296 advances. It occurred while resolving
curb heights. The complete local calculation took 51.45 s. All UVs were
finite and the overlap check passed for 14,786 candidate pairs. This is
a numerical producer measurement; it does not measure Blender viewport
frame time. See `nov20_fixed_1386_result.json` and `nov20_fixed_1386.profile`.

Beijing junction 1145 had another synchronous call inside its longitudinal
chart recovery. Forwarding that generator's existing checkpoints reduced the
longest single advance from 92.19 s to 3.80 s on identical saved inputs from
the 1,072-Section component. Both replays returned the same 4,581 pavement
quads. Vertex positions, face indices and UV buffers matched byte for byte,
and the overlap results matched. The measurements describe one numerical
junction task. See `nov20_beijing_1145_baseline.json` and
`nov20_beijing_1145_current.json`.

The same component subsequently reached Beijing Bundle 16. One private
separator bank ended 3.75 mm after a regular sampling point. Its pavement
could not satisfy the 8:1 quad limit at the shared edge. The failed-only
smooth-bank retry redistributes the last two sampling intervals while
keeping the terminal profile and plateau height fixed. The skipped point
must project onto the replacement chord within one micrometre or one
thousandth of the configured curve tolerance, whichever is smaller. Both
grass and pavement consume the same resulting bank. The measured deviation
at the troublesome end was 0.285 micrometres. A separate failed-only tile
retry can place a computational cut on a nearby existing crown sample when
the original cut requires more than 64 transverse intervals.

The complete Bundle 16 provider, replayed with the original 1,072-Section
component inputs, produced 28,065 pavement and 6,180 grass quads. The stored
maximum edge ratios were 7.878 and 7.840. The UVs were finite and the local
overlap check passed for 96,361 candidate pairs. See
`nov20_b16_whole_provider126.json`. This test covers the complete Bundle
and its internal interfaces. The city import also checks neighboring outputs.

For Beijing Bundle 5, checking 13,980 quads against the same two original
strips and their enclosed hole took 126.19 s with the scalar check and
7.08 s with the new check. When an expanded quad bounding box contains no
boundary segment, the whole box lies on one side and a single classification
suffices. Ambiguous boxes still receive the continuous edge and enclosed-hole
tests. Both runs returned the same proof and error allocation. The comparison
times cover containment alone. See `nov20_full_strip_proof_b5.json`.

The next Beijing and Bergamo cases were reproduced on their original inputs.

- Beijing separator 1 has one end attached to a junction and one free end.
  The free end previously copied the opposite junction's 296 intervals,
  giving 16.7 mm edges beside road-edge intervals of about 0.70 m.
  The provider now samples a free end using its own length and the requested
  mesh spacing, with the parity required to fill the rim with quads.
  The 4.94 m free end uses eight intervals. On the original component inputs,
  all 296 attached-end intervals and all 178 road-edge intervals appear as
  individual mesh edges. The 1,520 grass quads pass the quality checks before
  and after conversion to stored coordinates, with a maximum edge ratio of
  7.50 or less. All 6,080 UV corners are finite. Five of 1,254 available
  neighboring meshes have overlapping bounding boxes. Testing the separator
  against those five found no overlapping faces. The component-level replay
  is recorded in `nov20_separator1_provider127.json`.
- In Bergamo, the ordinary-junction sweep reaches junction 2857. The corner
  between two 2 m streets (875 and 1429) cuts through the original 7 m
  width of Section 874 near its root. The union exposes fragments of the
  transverse end caps of both wide Sections, 874 and 1427. Part of an
  original circular curb lies between those fragments. A finite-difference
  endpoint tangent put the trial join across the last emitted bank chord
  near the end of the curve. For two consecutive Smooth Sections of the
  same Bundle member the builder can use their original analytical tangent
  directions and choose local contact vertices on the sampled banks. The
  original four street widths, profiles, mouth planes and area operands
  remain part of the qualification. The real captured group passes the
  continuous exposure and width checks with contacts at 1.23653 m and
  0.12564 m and a 1 m certified minimum radius. The full ordinary-junction
  footprint builder on the original 688-Section inventory then returned one
  closed outline. This check ends before surface meshing and city publication.
  See `nov20_2857_integrated128_run`.

At Bergamo junction 1485, the common mouth at the boundary of the supported
ground-road network had the correct bank elevations but assigned zero crossfall
where the two retained roads supplied 5.93%. The mouth now takes its endpoint
gradients from those roads and eases them to zero at interior bends and the
crest. The real surface passes with 593 pavement quads, a maximum edge ratio
of 7.53, finite UVs and no overlaps among 1,575 candidate face pairs. The four
finite-scope mouth regressions pass. See `nov20_1485_crown130b` and
`nov20_finite_scope130b`.

Continuing from that surface on the 1,255-Section Bergamo input exposed a
second failure at junction 1930. Its dead-end pavement has two different
curbs. The through corner has two circular tangent contacts. The cap on the
opposite side runs along retained road edges and a generated terminal curb.
The height solver previously searched that cap for the through corner's
contacts. They lie 2.47 m and 3.45 m away from its nearest vertices.

The union now records which original road bank supplied each surviving cap
edge and the stations covered by that edge. Those banks supply their physical
heights and crossfalls. The generated portions join the surviving banks and
mouths at the common curb elevation. On the real input, 202 retained edges
carry bank provenance. The footprint is identical vertex for vertex. The
complete surface has 2,465 pavement quads, a maximum edge ratio of 7.835,
finite UVs and no self-overlap in 9,424 candidate pairs.
See `nov20_1930_native132_run`.

The next downstream Beijing replay completed separators 1 through 9. At
separator 10, eight consecutive bank intervals cover just 20 mm. The
adjacent edges measure approximately 0.49 m and 0.47 m. A uniform interior
row produces long, thin quads at these two transitions. A local series of
quad bands now grows around the fine run, with each band sharing the same
two endpoints. The last band meets the ordinary separator mesh. On the real
rim this gives 6,739 quads, keeps all 1,666 original boundary vertices and
passes the overlap tests in both working and stored coordinates. The maximum
edge ratios are 7.963 and 7.986 respectively. UVs are finite. See
`nov20_separator10_verify129_run`. This is a separator result; the replay
does not establish whole-city coverage or progressive-import parity.

The subsequent component replay generated the separators through 17.
Separator 18 has a 1.4° grass tip. The two carriageway banks supply different
short intervals at that tip before returning to their ordinary spacing of
about 0.7 m. When ordinary charts fail, staggered rows can connect the two
banks while preserving each physical edge. On the original rim the complete
separator provider returned 2,153 grass quads, all 608 boundary vertices,
finite UVs, and maximum physical and stored edge ratios below 7.50. Both
self-overlap checks passed. See `nov20_separator18_verify131_run`.

The Beijing replay then built separators 18–20 and stopped at Bundle 21.
On one carriageway a short visible stretch of Section 156 enters and leaves
the same compound junction. Its four inner edges face an 80-edge uninterrupted
junction curb. Together they enclose a 48.39 m² grass island. This produces
three visible runs on that member and two on the opposite member. There is
no free terminal to extend.

The separator builder accepts this interior island when the two adjacent
regions on both members meet the same physical owner and its curb connects
the selected original bank vertices. The island and both adjacent regions
keep their original boundaries. The complete real Bundle has 9,414 grass
quads, finite UVs and a stored maximum edge ratio of 7.762. Checking it against
the 48 available neighboring road and junction meshes found no overlaps in
either working or stored coordinates. See `nov20_beijing21_provider133b_run`.
The 18 finite-cap and interior-Bundle regression tests also pass, including
the mirrored island, rejection when the neighboring junction identities differ,
reversed road directions and rigid transforms. See
`nov20_pocket_profile133_tests_run`.

The overlap validator collects city geometry directly into its final arrays
and quantizes bounded batches of quads. On a 131,072-quad test grid, peak
traced allocations for a whole-network check fell from 42.1 MiB to 28.2 MiB.
For a single changed road they fell from 37.1 MiB to 23.2 MiB. The old and
new checks returned identical results on both runs. Those measurements
exclude the input meshes and measure validation allocations only. The SAT
and publication regression suites also passed, including near-contact cases.
See `nov20_measure_overlap_memory129_run` and `nov20_overlap_budget129`.

The full native Beijing run in `nov20_native_beijing_progressive_guard1`
stopped after 38 min 42 s when free system memory fell below the 1 GiB
reserve. Blender used at most 681 MiB of sampled private memory in that run.
It stopped during the blocking comparison and yielded no city-wide parity
result. The next native check reserves the machine's street-test workers
for Blender and records the current numerical task and periodic stack
samples.

Separately, on the captured four-arm junction 1683, forwarding checkpoints
through the original corner-envelope calculation reduced its longest
uninterrupted advance from 10.81 s to 0.031 s. Total calculation time was
10.21 s. All 310 outline vertices, source labels, contact records and samples
matched the earlier run. See `nov20_original_root_measured1`.

In the 688-Section Bergamo run, one numerical advance took 17.04 s. Its
stack samples repeatedly showed the complete carriageway-width test. Revision
119 groups finite segment queries that use the same polygon. Each segment
keeps its original intersections, interval midpoints and tolerance-band
coverage test. Only the evaluation of the point queries is batched. On the
same 2,672 real source cells, the width check took 4.19 s, compared with
8.62 s in the previous measurement. The width report was identical, as were
the coordinates, heights and gradients on both boundary loops. These timings
measure this check alone. See `nov19_cached_source_proofs_result.json`,
`nov20_batch_width119_report.json` and `nov20_source_proof_regressions119`.

The progressive version of that same width check returned control 222
times. Its longest advance was 0.33 s. A separate qualification run on
Bergamo's complete 25-Section component saved all 91 completed tasks after
ownership planning, then resumed in a fresh process. The remaining 71 tasks
produced the same 31 parts as the uninterrupted run. Vertex and face buffers,
material names, UVs and the complete overlap report all matched exactly.
The restart harness checks the original network fingerprint and the hashes
of all frozen numerical sources before accepting saved work. This qualifies
the numerical test harness. Blender publication is checked independently.
See `nov20_688_width_steps119.json`, `nov20_25_cold119/report.json` and
`nov20_25_resume119/report.json`.

A second restart check covers the costly ownership-planning phase itself.
On the seven original Sections around Bergamo junctions 2831, 2871 and
4713, a fresh process restored 26 finished tasks and the exterior computed
for those three junctions. It ran the remaining 21 tasks. All six resulting
mesh parts matched the uninterrupted calculation in positions, faces,
materials, UVs and the full overlap report. Reuse requires the same original
layout objects and dependency digests after deserialization. The production
addon has no dependency on these test files. See
`nov20_orig7_exterior119_1/report.json` and `nov20_resumed7_exterior119_1/report.json`.

The numerical qualification of the 1,255-Section Bergamo component also
passed 20 more individual owners between 1,089 and 1,386 before the memory
reserve interrupted the next build. The Beijing and Bergamo full imports
still require end-to-end validation.

### Numerical checks, 19 November

| Location | What was verified | Scope |
| --- | --- | --- |
| Beijing junction 9 | The Python junction task in the development checkout produced the same 5,072 pavement quads and vertex coordinates as the independently audited mesh, including the coordinates stored for Blender. Its four approaches have genuine Bundle half profiles. | Used the ownership plan saved by the earlier empty-cache component run. |
| Bergamo junction 210 | The rounded join between a full crown and a half crown passed independent checks on 1,958 pavement quads. The maximum stored slope was 0.13113 against a limit of 0.14080. All 18 source identity and serialization checks passed. | Local junction construction. The maintained ownership tests exercise both orientations of the endpoint join and the join between two exposed carriageway ends. |
| Bergamo junction 236 | The default footprint clipped Section 108. Restricting the curvature bound to the actual owned interval made that defect detectable. Using per-piece curvature when the complete station grid exceeded its budget allowed the existing mouth search to correct it. The resulting footprint passed continuous width checks across all three Sections in 2,282 cells at the existing 2 mm tolerance. | Reconstructed from the failed component's original planner. The correction used the existing 20,000-point grid and 8,192-cell width-check limits. |
| Bergamo junction group in the 688-Section component | An independent audit of the experimental constructor output passed on 10,613 pavement quads and 3,639 grass quads. All 51 sampled profile positions and heights matched exactly. The remaining island covered 34.806 m². | The constructor received prepared junction geometry. This check covers the resulting mesh, including stored slopes, UVs, overlaps and height queries. Complete component preparation and Blender publication remain open. |

The reports for these checks are under
`artifacts/streets/beijing_bergamo_fixes_20261114`. The full-city acceptance
criteria include native publication of all supported ground components and
blocking/progressive parity.

#### Where the two component runs stopped

The Beijing numerical run built junction 9 and stopped at junction 68. At 68
the ordinary junction owns two divided-road half profiles and a full 2 m
approach. After ownership finalization the original single-site correction is
nested inside a summary record. The development checkout now reads that
record and verifies its site identity. A regression exercises the real
finalization function and rejects a mismatched site or duplicated record.
Replaying the junction still hits the quad layer limit at both supported cut
spacings along Section 35. Inspect the failing subregion and its boundary
before changing how this type of junction is divided into quad patches.

The Bergamo 1,255-Section component advanced past junction 236 and stopped
during footprint construction at junction 237. About 1.74 m of the transverse
end of Section 119 appears on the outside of the union. The other two
approaches have mouths about 133.5 m along their paths. Trace the adjacent
source banks and their intended join before deciding which boundary should
replace the exposed end. Increasing a meshing budget cannot resolve this
footprint question.

Reproduction data are in the artifact directory named above. Use
`nov19_beijing448_full01/ACTUAL_validated_failure_prefix.pkl` for Beijing and
`nov19_1255_full02/failure.pkl` for Bergamo. Both component runs stopped before
Blender publication.

#### Verified local constructions for the next implementation pass

Both saved failures now have numerical solutions on their original input data.
Each experiment supplied the candidate geometry to the development checkout's
`CompoundProduct.build(..., source_heights=True)`. Automatic selection during
import and a fresh run of each complete component are still pending.

**Beijing junction 68.** The layer limit was reached in a block with 78 rim
vertices. A 1 m-radius corner met a carriageway edge just beyond an existing
corner vertex. This left a 2.251 mm chord next to a row of 42.505 mm chords.
The corner producer had retained the analytical circle, its finite angular
interval and the mapping back to the original prepared curve. Using that map,
the sampler located the intersection with the adjoining curb again, then
spaced the vertices over the visible arc. The maximum angular step was bounded
by both the original corner sampling and the active spacing/error settings.
The existing mesher could then fit quads within its normal budget.

The contact moved 0.366 mm. The shortest arc chord became 39.120 mm and the
maximum arc-to-chord deviation decreased from 0.226 to 0.191 mm. An independent
check compared all 48 mouth edges, including coordinates and transverse
parameters, and verified a conforming disk. The resulting 2,058 pavement quads
had a maximum stored edge ratio of 7.755 and a minimum stored corner sine of
0.0401. The maximum measured physical slope was 0.16584 against 0.17600 allowed.

For integration, use the producer's declared circle only when clipping leaves
a short generated terminal chord. Retain genuine contacts and original input
vertices. Bound sampling on the exposed arc independently of its approach
segments. Recheck the complete source widths, mouth coordinates, topology and
stored mesh quality on the final geometry. The diagnostic touched one circle.

**Bergamo junction 237.** A fork angle of 8.79° puts the 10 m-radius turning
corner about 133 m along Sections 86 and 87. Those distances establish the
profile planes for the fork. They are a poor length for the width transition
on the opposite side. Starting at the 10.5 m-wide Section 119 bank at the
original node and ending 10 m along the 7 m-wide Section 86 gave a smooth
transition. Beyond these contacts the boundary follows each actual carriageway
edge to its existing profile plane. The other two corners keep their original
construction. The maximum certified curvature of the transition was
0.09644 per metre against 0.1 allowed. Complete carriageway widths passed
at the existing 2 mm tolerance over 538 cells.

The two ends of this curb have different elevations. The test preserved the
original heights and crossfalls on the approach segments, then blended the
0 m road edge into the 0.08 m median edge over the quintic transition alone.
For a curve parameter `t`, compute the physical along-curb grade as
`(dh/dt) / ||dC/dt||`, where `C(t)` is the actual planar continuation.
Combine that grade along its unit tangent with the interpolated inward
crossfall along its unit normal. On the adjacent road edges use
`SourceOffsetRule` so the evaluator sees the original path's parameter
intervals. The boundary-field fit passed the existing tolerance and the
surface builder accepted 5,390 pavement quads. The maximum physical edge
ratio was 7.950 and the minimum corner sine was 0.247.

For integration, select this local transition after an ordinary continuation
fails the complete-width check. Determine the reach from the actual source
widths, the permitted radius and available path length. Both directions of
the width change need tests. Reuse the original profile planes, source
offsets and certified continuation. Include the selected transition and
sampling policy in the dependency signatures before checking updates.

The reproducible local checks are in the same street artifact directory.
For Beijing use `nov19_geom68_product.py` and `nov19_verify_68_geometry.py`
with the saved component prefix. Their successful runs are
`nov19_max_68_product02` and `nov19_max_68_audit01`. For Bergamo use
`nov19_geom237_product.py` with the saved planner state. The successful run
is `nov19_max_237_product03`. These local results preceded automatic backend
integration, full Blender imports and progressive/blocking parity checks.

#### Automatic selection in the Python backend, 19 November

The regular junction task now recovers Beijing junction 68 from a quad layer
limit by resampling the visible part of the original circular fillet. On the
saved component plan, the task reproduced all 2,205 physical vertex positions,
the same 2,058 quads in the same order, and all 2,205 publication positions
exactly. The full-width and stored-mesh checks run on this new sampling.
The comparison is recorded in `nov19_high_exact68_02/stdout.log`.

For Bergamo junction 237, the ordinary footprint construction first tries its
existing width transition. If the polygon union still exposes the transverse
end of a Section, it tries a transition on the adjacent banks and certifies
the complete footprint again. This ordering preserves a previously successful
fit on the ordinary width regression. Automatic preparation and meshing of
the saved Bergamo junction produced 5,390 quads. The selected narrow-bank
contact was 10 m from the junction on the eighth candidate tested. See
`nov19_high_auto237_03/stdout.log`.

These are individual junction results. On a fresh run of the maintained
1,072-Section Beijing component, planning stopped before junction meshing.
At the join between Sections 428 and 431, only 44.77% of the proposed fillet
lies on the exposed boundary. The earlier experimental component used a
tangent connection between the incoming Section's bank and the visible part
of the original 2 m-radius corner. A separate replay of the fresh failure
reproduced the two tangent contacts. With the experimental fallback installed
in the diagnostic process, the four-junction group completed exterior
preparation and passed the width check on all nine covered Sections. That
construction still needs integration and qualification in the maintained
import path. The failure capture is in `nov19_high_beijing_full01` and the
successful diagnostic is in `nov19_high_861exterior01`. The subsequent
maintained implementation and its results are described below.

#### Beijing bank tangencies and native publication regressions

The Python backend now uses the visible portion of the original circle when
the ordinary bank-to-bank fillet fails its exposure check. At the Beijing
group containing junctions 856, 861, 935 and 1799, the accepted 1.6 m fillet
is tangent to Section 428 and to the original 2 m corner. Its Section-side
contact is at station 124.62590 m. Its circle-side contact is at 0.78352145
of the original arc. The closed footprint passed the original width checks
on all nine covered Sections.

The junction mesher lays out computational cuts along a fully owned,
undivided Section. It checks the resulting cuts against the original mouths
and ownership intervals. The maintained planner and product constructor
produced 13,685 pavement quads and 2,367 grass quads. Both meshes match the
independently audited construction exactly in face order and in every
physical and publication coordinate. The maximum pavement edge ratio after
storage was 7.891. The comparison and the slope, UV and overlap checks are
recorded in `nov19_high_native861_product02`.

The next empty-cache Beijing attempt reached junction 1014. The normal radius
search oriented a corner using a distant endpoint of Section 507, after that
Section had changed direction. The local fallback uses the directed
tangents of the two *exposed* banks. It searches within their visible station
intervals and checks their analytical tangent contacts and actual bank
crossing. The original 3 m style radius now fits on Sections 508 and 507.
The adjoining join with Section 71 still uses the ordinary construction.
Preparation and meshing produced 2,336 pavement quads at this junction.
See `nov19_high_1014_01`.

The isolated Blender 5.2 run in `nov19_progressive_current01/checks.json`
passed all 17 fixture tests, including exact blocking/progressive parity,
terrain, cancellation, scene changes, and preservation of user edits.
Those tests used the fixture inputs. Whole-city publication and real-extract
blocking/progressive comparisons are separate acceptance checks.

The next Beijing run reached junction 1229. An 18.4 mm segment of the
transverse end of Section 581 remained on the exterior. The numerical
replay `nov19_high_crest1229_01` used the two consecutive median-side
banks of one Bundle member, at their common crest elevation, to support
the adjoining corner. It produced 2,429 pavement quads with a maximum
edge ratio of 7.591. This replay used the experimental crest construction.

The maintained backend now admits this construction after the ordinary
outer-curb search exhausts its candidates. Both crest banks must belong to
consecutive Sections of the same Bundle member, meet at their declared
ends, and have the same physical height and zero gradient. The actual
incoming mouth positions bound the two contacts. The existing cyclic
search verifies the whole curve remains exposed and the transverse ends
are inside the pavement. The planner records the source paths, profiles
and approaches for invalidation after edits.

In `nov19_native_crest_replay05`, the normal planner selected Bundle 52,
member 1. Both contacts were at 1.75 m and the crest height was 0.080 m.
The 2,429-quadrilateral mesh passed the UV and overlap checks. After float
storage, the maximum edge ratio was 7.5830, the minimum corner sine was
0.08346 and the maximum slope was 0.09329 against a limit of 0.1. This
run replayed the accepted contact sampling on the final source ports.

The following Beijing run reached junction 2042 on Section 1021.
A trial curve touched another source polygon at one point. The search now
rejects that candidate and continues its finite search. Final unions still
reject point-only contact. The next construction exposed a separate issue
on the reversed route of Section 600. A generated sample and an original
circle contact differed by 5.33 × 10⁻¹⁵ m in outward station and mapped to
the identical native cut. The backend uses the original circle's contact
station when it can verify the source frame and complete cross-section.
It preserves source controls, other tangent contacts, mouth planes and all
existing polygon coordinates. The maintained product in
`nov19_2042_final_product` has 1,537 pavement quads, a stored maximum edge
ratio of 7.381 and a maximum physical slope of 0.17077 under its incoming
profile limit of 0.176. Its UV and overlap checks passed. The subsequent
cold full-component run reached junction 2089. A dimensional corner and
the incoming Section 1042 sampled the same curved bank at different
stations. They now share physical coordinates at 81 stations over
5.1891 m, including the original controls and tangent contacts. The
resulting 78.5945 m² exterior produces 1,430 pavement quads with a
maximum stored edge ratio of 7.712 and a maximum physical slope of
0.07378 under the 0.1 limit. UV and overlap checks passed
(`nov19_2089_maintained_mesh01`). The three-Section regression in
`tests/test_python_street_shared_bank_sampling.py` exercises the update
from its centreline inputs through the final mesh.

The next empty-cache Beijing run (`nov19_high_beijing_full07`) passed
2089 and stopped at junction 2383, between Sections 528 and 582. A 2 m
radius fillet there turns by 0.02575 degrees. The maximum distance
from its arc to its chord is 0.0000505 mm, below the polygon operation's
usual welding distance. A local retry now resolves its two contacts and
checks the result independently at two precisions. Every surviving edge
outside the two replaced bank spans has the same coordinates and source
as in the original footprint. The raw arc samples exhaust the quad
chart's layer budget. The backend now uses its existing exterior
resampler, with at most 0.001 mm additional error on the sampled
outline. It retains the analytical source operands, the complete
incoming profile samples and the original boundary height rules.
Junction 2383 produces 671 pavement quads. Its stored maximum edge
ratio is 7.934 and its physical maximum slope is 0.15849 under a limit
of 0.176. UV and overlap checks passed (`nov19_2383_maintained_mesh02`).
The 17 isolated Blender progressive
fixture checks passed after the 2089 fix (`nov19_shared_progressive01`).

The empty-cache Bergamo run `nov19_high_133_full01/report.json` completed
all 972 numerical tasks for the 133-Section component. It committed the
update after UV validation and overlap checks across 238,655 quads.
In the separate 688-Section component, `nov19_high_688_pair01` isolated
the failing island corner to the group containing junctions 624, 630 and
632. Each of the five allowed radii from 1.0 m down to 0.5 m failed the
same test: no circle of that radius could be internally tangent to both
original circles. In a separate numerical replay, the experimental search
identified the exposed circular and native bank intervals directly from the
failed import. It rebuilt both joins and left a 34.8056 m² island from the
original 35.5849 m² footprint. One 1 m arc joins the two native banks.
A 0.64 m arc joins the retained circle to the native bank at the island.
Across all 19 transitions the largest discrepancy between analytical
tangents was below 0.000001°. See `nov19_688_rederived_geometry.json`.
The maintained offset solver independently recovers all four tangent
points to within 6 × 10⁻¹³ m given those bank intervals and radii
(`nov19_688_native_solver01`). The remaining integration work includes
discovering the chains, binding their original height rules, and checking
the mesh in the full component. Passing the experimental footprint to
the current height provider stops at the first new arc, as recorded in
`nov19_688_current_boundary01`.

The 1,255-Section run `nov19_high_bergamo_full01` passed ownership planning
and stopped at junction 167. Its generated corner connects the inner bank
of Section 79 at 0.080 m to an outer bank at zero over 0.320 m. Its mean
slope alone would be 0.25, above the 0.10057 limit. The zero-height end of
that sampled span is 1.280 m from the incoming crest at the mouth.
That distance leaves room to investigate a transition inside the junction.
The opposite-bank tangent on Section 79 falls behind its origin, outside
the positive station interval required by the existing transition method.
The boundary measurements are in `nov19_bergamo167_bounds04`. For this
case the backend now fits a smooth curb from the Bundle's inner bank to
the zero-height outer bank. It uses the existing incoming mouth planes.
The original quintic can cut into either carriageway. One extra control
in the curve lets it move outward while preserving the endpoint
positions, tangents and second derivatives. The chosen offset is
17.5 mm. The complete-width checks pass and the curvature bound is
0.20890 per metre, below the limit of 1/2.62144 per metre. This fit is
available only for a style radius. Rebuilding the junction from its
original inputs produces 3,282 pavement quads with a stored maximum
edge ratio of 6.945 and a maximum slope of 0.09987 under the 0.10057
limit. UV and overlap checks pass (`nov19_167_maintained_mesh01`).

Beijing run `nov19_high_beijing_full08` completed ownership planning
for its 1,072-Section component. It built seven junction surfaces and
stopped while assigning boundary heights to the combined surface for
junctions 91, 99 and 579. This run used the first version of the 2383
footprint repair. It did not reach 2383 in surface generation.
Section 43 is fully absorbed by that group. Its preserved 2.7 m-wide
profile has an 8 cm crown with an 11.85% crossfall at the bank. The
outer entrances currently give the group a 10% slope limit. A longer
transition cannot resolve that difference at the bank itself
(`nov19_beijing91_source04`).

An experimental comparison on the same group tested the effect of including
the retained street profiles in the slope calculation. Using the existing
1.1 margin gave a 13.037% allowance. The complete surface builder produced
2,544 quads. The largest stored pavement slope was 11.722%, and the mesh, UV
and overlap checks passed. A second test allowed the same original boundary
gradients but kept a 10% limit on the pavement. It failed the coupled slope
and drainage checks.
See `nov19_beijing91_hypothetical_surface.json` and `nov19_beijing91_strict_mesh.json`.

Revision 113 implements the retained-profile rule in the Python backend.
For a fully owned Section, the maximum physical crossfall of
its parabolic crown is the larger absolute edge derivative in the normalized
profile coordinate, multiplied by 2 and divided by the width. The junction
applies the existing slope margin to the greater of that crossfall and the
measured entrance slopes. It checks the result against the existing minimum
allowance. The separate crown-height, drainage and longitudinal entrance
checks still run. The production builder reproduced all 2,544 quads at
owner 91 with a stored maximum of 11.722% and 6,823 candidate overlap pairs
checked. The regular product provider also resolves junctions 99 and 579 to
that same surface. The separate owner 315 passed with 3,432 quads.
See `nov19_113_beijing_production.json` and `nov19_related_slope_owners.json`.

### Real Blender publication on small components, 19 November

Blender 5.2 imported the complete original OSM file for each city. For this
check the harness selected connected components containing at most 24 Sections
after Blosm's usual network partition. Both import modes received the same
components and used the normal mesh builder, materials and Blender publication
code. Larger components still need separate qualification.

| Extract | Committed in progressive mode | First visible preview | Longest measured callback | Blocking comparison |
| --- | --- | --- | --- | --- |
| Beijing | 9 of 9 selected (16 Sections) | 6.44 s | 4.48 s | Exact match of source metadata, geometry, normals, materials and UVs in Blender |
| Bergamo before the width correction | 29 of 30 selected | 12.97 s | 26.47 s | Both modes reported the same preparation seam. The exact mesh comparison stopped on the partial result. |
| Bergamo after the width correction | 30 of 30 selected (45 Sections) | 12.08 s | 100.50 s | Exact match of source metadata, geometry, normals, materials and UVs in Blender |

These background Blender timings include source preparation. The long callbacks
can still interrupt interactive redraw. The test did not measure foreground
frame rate. The corrected Bergamo run took 233.78 s in blocking mode and
232.91 s in progressive mode. Its first committed component appeared at 13.92 s.
The maintained comparison harness disabled the import performance request.
These results establish publication parity and show that long uninterrupted
callbacks remain a responsiveness problem. The complete metrics and source
inventories are in `nov19_real_slice_blender01/checks.json` and
`nov19_real_slice_blender02/checks.json`.

On the same 30 Bergamo components, the ownership-planning checkpoints reduced
the longest observed callback from 100.50 s to 62.29 s. Adding checkpoints
inside junction surface generation reduced it further to 29.35 s. All three
runs had the same exact scene digest, and blocking and progressive publication
matched within each run. Progressive elapsed time was 232.91 s before these
changes and 237.87 s with both sets of checkpoints. The change improves the
intervals at which Blender can respond. Several uninterrupted operations still
take tens of seconds. See `nov19_real_sliced_plan_blender01/checks.json` and
`nov19_real_sliced_surface_blender01/checks.json`.

The Bergamo failure also reproduces outside Blender in the four-Section
component containing Sections 114, 115, 116 and 122. At junction 230 a 10.5 m
approach divides into 7 m and 3.5 m approaches. The ordinary layout clips the
start of Section 115. The fallback footprint exposes about 1.74 m of an internal
construction edge on the pavement outline. The single-component reproduction
is `nov19_bergamo_114_full01`.

The corrected width transition can advance the wider approach's mouth when
the opposite end of that Section has no junction. It preserves the original
cross-section at the new mouth, certifies the pavement against the full widths
of the incoming Sections, and checks the curvature of the new bank. The
candidate search is bounded to one carriageway width, stays on the original
Section and applies only to style-defined transitions. This component passed
all 26 numerical tasks from a fresh start, including the UV and overlap checks
(`nov19_bergamo_114_full02`).

The two other Bergamo components with 25 and 39 Sections completed all 162 and
296 numerical tasks respectively, including UV and overlap validation, in
`nov19_bergamo_25_full01` and `nov19_bergamo_39_full01`. A separate Blender
5.2 import selected those same two components from the complete OSM extract.
Both committed with zero errors and matched the blocking import exactly. The
first preview appeared at 34.48 s, the longest callback took 3.87 s and
progressive publication finished in 95.27 s
(`nov19_real_medium_blender01/checks.json`). Together the two Bergamo
publication checks cover 32 of 35 components, totaling 109 Sections. The
remaining three components contain 133, 688 and 1,255 Sections.

In the separate 688-Section Bergamo component, the maintained boundary tracer
reproduced the two proposed infill closures from their refined source operands.
It reproduced all 11 original spans at the exterior corner and all 132 at the
island corner, including each span's parent and parameters. Omitting the
intermediate curb rejected both constructions. This checks the extent of the
infill. Automatic construction and publication of that component are still open.

### Consumed carriageway at Bergamo junctions 199, 201 and 202

The empty-cache 1,255-Section run passed junction 167 and stopped later at
199, which shares a pavement with junctions 201 and 202. At junction 202,
the curb rises 8 cm from the ordinary bank to the plateau of the divided road.
The circular corner supplies only 0.428 m. With the incoming crossfall and the
10% physical slope limit, the transition needs at least 1.687 m.

Section 101 is entirely inside this combined junction. Its edge still carries
an exact map to the original carriageway. The corrected height calculation
uses the flat plateau edge from 0.597 m to 1.908 m along the Section. Both
ends of that interval are existing mesh knots. The calculation checks the
original contacts on both sides of the Section at both junctions and stops
before a competing corner can use the same bank. The body-height check still
rejects an overlap with a required interior carriageway profile.

The complete transition measures 1.775 m. Its continuous boundary-gradient
bound is 9.61%. The rebuilt group has 9,743 pavement quads with a maximum
slope of 9.999% after conversion to Blender precision. The worst edge ratio
is 7.84 and UV and overlap checks pass. This check uses the captured ownership
plan. The full city still needs publication qualification. The independent
surface check on the current code and the preceding cold run are recorded in
`nov19_199_maintained_mesh03` and `nov19_high_bergamo_full03`. A physical
height plot of the stored mesh and a close view of the plateau edge are in
`nov19_199_physical_inspection.png` in the qualification artifacts.

The next independent surface run passed groups 205, 210, 213, 224, 236 and
237. At junction 261 ordinary basin filling tried to raise a prescribed point
on the carriageway by 0.0774 mm. The resulting conflict was between the fill
floor and the original profile. All corner slopes could already satisfy the
group's measured limit. On the unchanged quad graph, lowering one neighboring
free vertex by 0.0774 mm gives the prescribed point a two-edge drainage route.
The final edge reaches a point 0.1008 mm below the prescribed point.

The new recovery searches for a path across the existing spill elevation when
ordinary filling conflicts with an exact carriageway height. It retains all
prescribed heights, recomputes boundary drainage, and checks the height
envelope, entrance slopes and quad gradients. Group 261 now produces 10,230
pavement quads through the normal surface builder and passes independent UV
and overlap checks. Eight subsequent groups passed on the same ownership
plan. At group 362 the two exposed source-corner fragments meet at the same
zero height, but their incoming gradient magnitudes are 9.14% and 4.57%.
This distinct corner-contact failure is recorded in
`nov19_362_contact_failure01.json`.
Evidence is in `nov19_261_probe02`, `nov19_261_spillway01` and
`nov19_261_current01`.

At group 362 the source union exposed a 48.38° turn between an 8 m circular
curb and the straight, finite extension of the Section 1440 corner. A 2 m
tangent arc fits on both source segments and leaves 0.744 m of the straight
extension. It adds 0.1514 m² of pavement. The generator constructs the arc
from the original circle and the producer's declared extension. Both contacts
use their original road profiles. The curb stays at zero height and the
inward crossfall changes smoothly from 9.14% to 4.57%. The construction
checks the source widths, the complete arc, and every original entrance.

On the actual three-junction group from the saved ownership plan, the
maintained surface builder produced 5,706 pavement quads. The largest
physical edge ratio was 7.831. The result passed the UV and overlap checks.
The independent plan-view fit and physical product are recorded in
`nov19_362_wedge01.json` and `nov19_362_maintained_surface01.json`.
The test module `test_python_street_bank_circle_contacts` also exercises
both traversal orders and a finite lead too short to accept a contact.
The next 40 physical junction surfaces passed the same UV and overlap
checks. The following group, owned by junction 958, still fails. At the
connection from Section 478 to the rounded curb, the original crossfall
points slightly outward relative to the exposed corner. At the other end
it points inward toward Section 745. The source curves and their parameters
are recorded in `nov19_958_contact01.json`.
A diagnostic on the complete source group measured a clockwise turn of
94.98°. Joining Section 478 to the remaining part of the 1.024 m original
arc failed at every allowed radius. A join between the two exposed Section
banks passed the source-width and geometry checks and added 0.740 m² of
pavement. The maintained generator now finds this join from the two Section
banks and the intervening original corner. It tries the repair only when the
ordinary profile calculation identifies that same curb contact as invalid.
It compares the repaired hole against the union of the original area operands
before and after analytical sampling. Comparing the original sampled outline
against the new source polygon previously rejected a valid hole.

The source-group replay accepted a 0.65536 m radius, retained the island and
measured zero curb elevation and 4.57143% crossfall at both contacts. The full
surface then failed the final check on the elevations that Blender stores.
The first stored buffer reached a 10.02694% slope against the 10% limit.
After correction one vertex still needed a 1.097 μm fill on the eighth
iteration. The guard now tests the spill elevations themselves on that
last iteration. It retains the original slope, source height, drainage,
entrance and fixed-geometry checks. On the captured 10,241-quad mesh the
stored maximum was 9.09472%, the largest longitudinal entrance slope was
0.37607%, and both the physical and stored buffers drained completely.
The complete junction builder on the captured source group then produced
10,241 pavement quads with a maximum edge ratio of 7.80614. The mesh, UV
and overlap checks passed. Evidence is in `nov19_958_repaired_precision01.json`,
`nov19_958_maintained_product01.json` and `nov19_958_complete_product102`.
Replaying junction 362 through the same revised code also passed with 5,706
pavement quads and a maximum edge ratio of 7.83105.

In a separate Blender test on the intact 133-Section Bergamo component, the
previous revision stopped during preparation after 301 seconds with a hole
containment error. The revised blocking import passed ownership and the
junction surfaces, including 995. It subsequently stopped at Bundle 10
while constructing its grass mesh. No component was committed. Capturing
the editable network directly from Blender confirmed that all Sections,
junctions, connections, Bundles and settings match the earlier 972-task
numerical input exactly. Only the generated network UUID differs.
The affected region is the second end cap of Bundle 10. It has 40 physical
rim vertices. On one bank, chords around 4–5 mm long adjoin chords hundreds
of millimetres long. An interior edge common to two neighbouring boundary
faces can only satisfy an 8:1 edge-length limit on both if their two
exterior lengths differ by at most 64:1. This cap exceeds that ratio.

After the ordinary fixed-rim and banded charts fail, the builder tries an
interior row with an extra parallelogram at each abrupt change in spacing.
This gives each neighbouring boundary face an inward edge of its own.
Two or four computational bands fill the rest of the cap. The complete
disk must pass the original limits in both floating-point representations,
preserve the entire physical rim and pass a global overlap check. The cap
has 174 quads and a maximum edge ratio of 7.682. The full Bundle built on
the current upstream geometry has 5,108 pavement and 1,080 grass quads.
Both UVs and the combined 18,109 candidate overlap pairs passed.

The focused build also exposed a missing dependency. The Bundle consumes
the layout of a two-arm seam between original member Sections 595 and 617
when proving their road coverage. Its plan now declares the layouts and
source profiles needed for these original spans. This allows the plan to
run from just its declared dependencies, independently of the full graph's
evaluation order. See `nov19_bundle_grass_resumed_guard`,
`nov19_graded_cap_full_plan` and the mesh detail in `nov19_bundle10_graded_cap.png`.
All 32 numerical tests for interior Bundles, fixed grass rims, gapped bodies,
terminal unions and grass recovery pass. The full blocking run of this
133-Section component in Blender 5.2 on frozen revision 103 stopped after
372 seconds when available system memory fell below the 1 GiB reserve.
Blender peaked at 499 MiB resident. The run produced no final mesh and the
paired progressive run did not start. See `nov19_133_guard103_blocking/resources.json`.

The same intact 133-Section component completed blocking publication in Blender
5.2 on frozen revision 104 in 664.6 seconds. Blosm reported one committed
component and zero errors. Peak resident memory was 634 MiB. Progressive
publication completed in 709.9 seconds. Both modes produced exactly the
same mesh buffers, normals, UVs and materials on 153 output objects
(239,613 polygons). Both published manifests match the frozen backend's
972 declared tasks and dependencies. Every one of the 133 source Section
IDs appears on a road, a junction, or Bundle pavement. The short dead end
on Section 1566 has its entire 2.452 m interval inside junction 3130.
These checks cover this component of the full extract.

The first progressive mesh appeared after 100.3 seconds. There were 40
preview batches. The 95th percentile among 7,909 callbacks was 25.3 ms,
but eight callbacks exceeded ten seconds. The longest took 370.45 seconds
in junction 995. A satisfactory interactive import still requires reducing
those stalls.

The strict scene comparison failed on two metadata fields. Bundle 10
reported 123 overlap-test batches in blocking mode and 126 in progressive
mode. All other published fields compared equal. Revision 106 visits the
contributing junctions in sorted order. Rebuilding that Bundle in separate
Python processes with hash seeds 1 and 91237 produced identical pavement,
grass, UV and overlap results, including 122 batches over 95,064 candidate
pairs. End-to-end comparison on the revised code remains open.
See `nov19_133_native_comparison104.json`, `nov19_133_attribution104_blocking.json`,
`nov19_133_attribution104_progressive.json` and `nov19_bundle10_order106_seed*.json`.

Sampling junction 995 on revision 106 found a single uninterrupted advance
lasting at least 193 seconds. Most samples were in fixed-rim quad relaxation
during edge rotations. Revisions 107 and 108 let the progressive importer
yield between chart attempts, per-block construction, relaxation iterations
and line-search trials. The blocking caller drains the same iterator.
The chart and grading suites passed 85 tests, including cancel and restart
with exact vertex and face equality. Full-scene timing and parity for these
changes still require Blender measurements.

The complete revision 107 progressive import produced the same 153 objects
and mesh buffers as revision 104. Its first mesh appeared at 99.7 seconds.
The longest callback fell from 370.45 to 151.52 seconds, still at junction
995. Thirteen of 8,388 callbacks exceeded ten seconds. Total elapsed time
increased from 709.9 to 807.2 seconds in these runs. The improvement concerns
the worst stall. The two runs do not establish a throughput gain.
See `nov19_133_progressive104_vs107.json`.

The matching revision 107 blocking run has also completed. Independent audits
of both manifests verified all 972 declared tasks and the attribution of all
133 source Sections. The outputs, mesh buffers and published metadata matched
exactly across blocking and progressive imports.
See `nov19_133_attribution107_progressive.json`.

In a separate Python replay of junction 995 on revision 109, all 16,868
pavement quads and 5,187 grass quads matched the previous vertex coordinates
and face indices exactly. UVs were finite and all 69,018 candidate overlap
pairs passed. The longest advance took 36.15 seconds and included source
coverage and height checks. This measures a geometry iterator outside
Blender's UI. See `nov19_junction995_step_result109.json`.

The same replay on revision 112 retained exact positions and face indices
for all 22,055 quads in 21,350 advances. Its longest numerical advance took
11.54 seconds. This still leaves long callbacks to address. See
`nov19_junction995_step_result112.json` and `nov19_junction995_step_profile112.json`.

### Rounded source boundary in the 688-Section Bergamo component

Junctions 624, 630 and 632 share a pavement. Along its island two exposed
original arcs change direction by 90.44°. Neither arc has enough exposed
length to fit the ordinary circle-to-circle join. The fallback follows both
arcs out to their original Section banks. It retains part of one arc and
joins it tangentially to the opposite bank at a radius allowed by the source
street styles. The other arc lies inside the resulting pavement. A candidate
must cover that entire intervening arc, keep the original road entrances and
widths, and leave the island identifiable in the union.

On the captured ownership group, the maintained constructor found a 0.64 m
radius and an island area of 34.806 m². A second, exterior curb met a
Section bank at a 21.66° angle. Its circle occupied only part of the
prepared curve, which also carries a straight extension. The constructor
evaluates the tangent on the declared circular interval and fits a 1 m
bank-to-bank arc at that exterior turn. Both loops passed the source-height
checks with the existing 17.60% slope ceiling. The full pavement builder
on this captured group passed publication checks, UV generation and
overlap validation with 10,613 pavement quads and a worst edge ratio of
7.86736. The full 688-Section Blender import still needs verification.
See `nov19_688_all_turns102` and `nov19_688_two_surface_102.json`.

A fresh import of the 688 Sections on revision 105 stopped during ownership
planning after 2,818 completed tasks. In the group containing junction 821,
the estimated parameter on the last chord of a rounded join was 0.955727.
The exact contact was at 1.0 on that same chord. The old neighborhood used
all 47 segments around the pavement wedge, including its internal closure,
while the rounded edge itself had only nine chords.

Revision 109 recognizes the clipped terminal chord from its source segment
index, original sampling parameters and orientation. The adjacent source
curve still has to satisfy its own neighborhood, ancestry and distance checks.
Rebuilding the captured group verified six contacts, both islands and 958
source-width cells at the existing 2 mm bound. Eighteen retained-arc and
hidden-tangency tests passed. This result covers the eight-junction group.
The complete component and city require separate verification.
See `nov19_whole688_revision105/report.json` and `nov19_retained688_reconcile109.json`.

The full revision 109 replay passed junction 821 and then reached a different
four-junction group (903, 1163, 1166 and 1249). Its exterior accepts the
2 m rounded corner. Four tips on two interior islands reject every permitted
style radius down to 0.5 m. The islands occupy 0.987 and 6.294 m² and fit
within strips at most 0.156 and 0.435 m wide. Every rim span comes from the
two opposing inner banks of declared Bundle 13 (Sections 449/452 and
453/582). Both banks carry the same 8 cm plateau with zero transverse
derivative. The group carries the Bundle's centerline-disjointness proof.
These are genuine separators. How their tips are represented is a distinct
decision from filling a small residual hole left by rounded road corners.
See `nov19_688_second_obstruction.json`, `nov19_688_local_bank_holes.json`
and `nov19_bundle13_rim_guard/stdout.log`.

Revision 113 keeps the original taper on a narrow, closed Bundle separator
when the style fillet cannot fit. It checks both original facing banks,
their Bundle roles and orientation, equal plateau heights, the centerline
separation proof and the complete loop width. An explicit radius still
requires a fit. Rebuilding the four-junction exterior retained both original
island areas and the rounded 2 m exterior corner. See
`nov19_113_bergamo_production.json`.

All 108 original boundary vertices on the first island and all 162 on the
second kept their exact coordinates. All 810 tested boundary samples retained
the 8 cm plateau and zero transverse slope. For inspection, the plan view
and enlarged sections are in the artifact `nov19_bergamo_bundle13_tapers.png`.

In revision 114, the conforming pavement builder can split a six-vertex face
directly along a diagonal or into three quads around an interior point. On
the affected pavement face the stored edge ratio fell from 11.58 to 5.45.
The complete planar junction chart passed with 11,896 pavement quads and
198 and 280 grass quads. Both grass meshes match the pavement seam exactly.
UVs were finite and all 40,858 candidate overlap pairs passed.
See `nov20_114_planar_attachment.json` and `nov19_113_separator_heights.json`.

Sections 445 and 452 meet on the exterior curb at zero height. Their
crossfalls differ and the curb turns through 12.06°. The source banks
become neighbours on the exposed pavement outline. One original side street
lies between them in the junction's angular order. The two adjacent style
sectors specify radii of 1 m and 2 m. The builder tries the smaller radius
on the exposed banks and checks both contacts against the original curves.
This gives each bank its own tangent and crossfall at the join. The 1 m
fillet passed the original-width and boundary-height checks. The ordinary
2 m corner on the same junction and both grass outlines also passed.
See `nov20_114_curb_profiles.json`.

The final rounded outline produces a planar chart of 11,892 pavement quads
and 478 grass quads. The grass surface builder reproduced the original 8 cm
plateau at every grass vertex. Both complete 3D seams match exactly, with
200 vertices on the smaller island and 282 on the larger. All 41,124
candidate overlap pairs passed. See `nov20_114_final_island_meshes.json`.

A side road can interrupt a Bundle member more than once inside one
compound junction. The body-height builder now intersects the clear
station intervals on the two original carriageways, accounting for the
direction of each Section. It constrains a whole transverse profile only
within the common intervals. Each gap retains the existing smooth weighting
for the transition. On this chart, 4,410 interior vertices have prescribed
body elevations. The preliminary entrance check passed at 1.592% against
the 3% limit. This check precedes final pavement grading.

At this junction, filling the lowest depression would raise a prescribed
carriageway vertex by 0.58 mm. If the lowest pass is sealed by prescribed
vertices, the drainage search can now try a path across editable pavement
to an existing physical outlet. The two searches share the same work limit.
The source heights, slope and drainage checks still run on the resulting mesh.
Here a 1.27 mm cut on a two-edge path drained the depression. The 3D
builder on the verified chart produced 11,892 pavement quads and both grass
meshes. The stored slope check measured 17.194% against a 17.600% limit.
The final spill check required zero fill. See
`nov20_protected_chart_final_guard/stdout.log` and
`nov20_paved_publication_guard/stdout.log`.

These checks use the real junction and the accepted production mesh chart.
A fresh run through chart selection and final grading in one process stopped
at the 1 GiB system-memory reserve. Full-city progressive publication remains
open. Separately, the 688-Section component reached an acute exterior corner
between Sections 1416 and 1417. Its 2 m style radius and the permitted smaller
radii down to 0.5 m failed on the visible original banks. The grass-tip rule
applies only to a closed separator. See `nov20_measured_acute_fits_guard/stdout.log`.

#### Acute exterior corner across two junctions

Section 1417 is 9.04 m long. The 2 m corner at its first junction could
shrink to 0.45 m and fit entirely on that Section, below the 0.5 m style
floor. Its outer curb continues through a convex corner at the next
junction onto Section 1437. The pavement already owns all of Section 1417.
Revision 115 follows that exposed curb through the two original junctions.
It fits a circle between the visible banks of Sections 1416 and 1437,
using the smaller of the two participating style radii as its starting
radius. Both tangencies must lie on the finite original banks. The curb
between the contacts becomes interior pavement. The complete circular arc
must remain visible, with the closure completely inside the merged surface.

The original three-junction constructor fitted a 0.64 m arc and verified
the full widths across 650 source cells. All 52 sampled arc intervals were
exposed. The resulting 4,544 pavement quads passed the UV and overlap checks.
The steepest stored slope was 17.305% against a 17.600% limit. A separate
empty-cache calculation using the seven original Sections completed all 47
backend tasks and checked 39,545 overlap candidates across 14,773 quads.
These runs also exposed a cache-signature serialization error for connectors.
The signature now records the original junction and Section-end identifiers
as scalar values. See `nov20_acute_full_product.json` and
`nov20_acute_seven_full2/stdout.log`.

A regression fixture translates the fork and assigns new Section identifiers.
It exercises a corner where both offset curves meet at their shared authored
centerline endpoint. Using that exact zero station avoids a rounding error on
a reversed curve that placed its contact just outside the Section. The fixture
reproduces the 0.64 m radius and checks the full arc through the progressive
planner. All 21 connector, bank and cycle tests passed in
`nov20_successor_regression04/stderr.log`.

The same rule resolved an exterior curb at Beijing junction 343. Its
original 3 m style radius fit. The resulting 3D product passed with 3,088
pavement quads, a maximum edge ratio of 7.25 and 9,082 overlap pairs
checked. See `nov20_beijing343_product.json`.

#### Side street on one carriageway at Beijing junction 354

Section 178 joins Sections 126 and 177. The opposing carriageway, Section
172, passes the junction without a shared node. Its curb contributes to the
outline of the merged pavement. Previously the planar outline was available,
but elevation construction stopped because the curb had lost the link to
the profile that defined it.

Revision 116 carries each original section, bank side and station interval
through both clipping planes. Each surviving part of a rounded junction
corner retains its original parameter interval. On a level exterior curb the
crossfall changes between the two incoming profiles over the exposed circular
arc. On the uninterrupted opposite edge it changes over the distance between
the two merged profiles. The transition uses the existing quintic smoothstep.
The clipping intervals distinguish the end of each constant-slope lead from
the start of the quintic when fitting the boundary gradient.

Rebuilding from the real entrance data reproduced all three planar curbs
vertex for vertex (26, 100 and 26 vertices). The normal surface builder
produced 1,088 pavement quads. Its maximum physical slope was 17.487% under
the 17.600% limit. The maximum edge ratio was 7.873, the UVs were finite and
all 3,465 candidate overlap pairs passed. Eleven other junction surfaces in
the same batch also passed. The 47 related profile, curb, boundary-gradient
and separator tests passed. See `nov20_beijing354_physical.json`,
`nov20_beijing_next_owners115` and `nov20_onesided_combined116`.

Three further junctions reproduced the same missing profile attribution. The
same construction qualified junctions 531 (1,090 quads), 540 (1,317 quads)
and 547 (1,184 quads). Each retained its original three curb polylines exactly
and passed mesh, UV and overlap validation. A plan view of the physical
heights and mesh topology at all four junctions is in the qualification
artifact `nov20_beijing_paired_surfaces.png`.

A separate four-Section extraction also emitted the junction mesh during
progressive preparation. Cutting the Bundle to these four Sections produced
a self-touching separator elsewhere in that reduced input. Whole-component
coverage and Blender publication must be checked on the intact network.

The intact 1,072-Section Beijing run subsequently built the four affected
junctions and stopped at junction 553. A side street enters one carriageway
here and leaves 22.96 m² of grass between the carriageways. The two halves of
the grass rim occur after the three entrance curves in the boundary list.
Looking up an entrance at those two indices caused an indexing error.

The surface builder now traces both halves to the facing carriageway banks.
It uses their common elevation when both banks and the surviving parts of
the original corner agree in height and have zero transverse slope. External
entrance curves continue to use their own profile contacts. On the real
junction all 870 original rim vertices had the same 8 cm height. Both the
706-vertex and 128-vertex sampled rims preserved that height and zero slope
and passed boundary-field construction. The complete surface has 5,463
pavement quads and 434 grass quads. Both physical meshes contain the same
436 seam vertices. The maximum pavement slope was 17.600% and the grass
was level. UV and joint overlap checks passed. The 54 related tests also pass,
including rejection of a substituted outer bank and cancellation while
building a chart. See `nov20_553_qualified.json` and `nov20_final_sources117`.

The same original junction also exercises the progressive sampling path.
The two rows below compare identical input geometry in separate processes.
Each resulting exterior and separator loop matched byte for byte in XY,
elevation and gradient. Times measure the two sampling operations alone.

| Operation at junction 553 | Previous time | Revision 119 time | Longest revision 119 advance |
| --- | ---: | ---: | ---: |
| Strict rim sampling | 5.74 s | 1.39 s | 0.226 s |
| Retry with the remaining curve tolerance | 2.26 s | 0.73 s | 0.240 s |

The 54 source, sampling, slope and cancellation tests passed on this change.
See `nov20_553_sampling_118.json`, `nov20_553_sampling_119.json` and
`nov20_serial119_one_sided_tests_1`.

An independent check of the next 24 Beijing pavement owners passed for 23.
Junction 691 initially failed the fixed-boundary grass slope check. Its
steepest grass rim edge measures 7.40%. Three consecutive rim vertices
belonged to one quad whose implied corner gradient was 21.54%. In revision
118, when immutable corner gradients fail, the builder inserts a ring of
quads inside each affected face. The original boundary, including its UV
coordinates, stays fixed. The builder checks the full revised surface against
the same cell, slope, shape and drainage limits. A rim edge that already
exceeds the slope limit still fails. The real grass mesh needed four more
quads. Its maximum physical slope measured 9.777%. Both materials passed
UV and overlap checks over 7,877 pavement and 514 grass quads. See
`nov20_beijing691_published118` and `nov20_fixed_hole_tests118`.

Twenty more
junctions in Bergamo's 1,255-Section component passed individual mesh, UV
and overlap checks. These runs examine individual products from saved plans.
The complete city imports and their progressive publication remain under
qualification.

Boundary assembly now reuses the same native road piece's endpoint gradient
when checking adjacent curbs and writing vertex heights. The cache belongs
to a single boundary assembly. The complete-width test computes curb and
end distances only for polygon vertices that occupy a source cell.
On the same junction, all 3,980 boundary positions, heights and gradients
matched byte for byte, and the full 2,672-cell coverage report was identical.
One isolated comparison measured 7.55 s and 3.39 s for boundary conditions
before and after the change. These are timings for this one operation.
See `nov19_reference_source_proofs_result.json`,
`nov19_cached_source_proofs_result.json` and `nov19_source_proof_regression_guard`.

### Independent Beijing surface check

At owner 343, the two original curb spans from Sections 172 and 128 meet at
zero height. Their transverse slopes at that point are 1.48% and 2.96%.
The physical boundary builder detects the disagreement on both the original
and resampled rim. Both slopes individually satisfy the 10% limit. This
case requires reconciling the two profile derivatives where the curbs meet.
See `nov19_343_native_banks.json`.

Eight additional physical junctions passed the pavement construction, UV
and overlap checks on the saved 1,072-Section ownership plan. Their owner
IDs are 145, 154, 185, 215, 219, 227, 239 and 250. Junction 185 has
22,537 quads and a maximum edge ratio of 7.943. The machine's free-memory
guard stopped this batch at junction 279 before it produced a result.
The per-junction evidence is in `nov19_other_owners102_Beijing`.

Revision 103 also passes owners 288, 301, 309, 317 and 332. Owner 279
exhausted its nearest-height iteration budget and three refinements before
reaching a satisfactory surface. Revision 104 retries the existing coupled
slope and drainage solver on the original chart after those attempts. The
2,125-quad pavement reaches 9.9321% physical slope under its 10% limit,
0.1806% maximum entrance longitudinal slope, zero drainage fill, and a
maximum edge ratio of 7.349. All 5,576 candidate overlap pairs and the UV
check passed. The 28 related numerical tests pass, including a fixture
that deliberately exhausts the initial solve. See `nov19_confirm104/279.json`
and `nov19_deferred_suite104_guard`.

At owner 315 the generated
arc belongs directly to original intersection 735. A diagnostic verifying
the original incidence, profiles and corner radii reached a retained bank
with an 11.85% crossfall. This surface remains unresolved. See `nov19_inspect_original`.

Revision 105 passes the actual boundary and pavement of owner 325. The
2.16 cm station on Section 169 is both a 4 m-radius corner contact and
the end of a clipped outer bank. The original road crossfall is 1.48%,
while the merged carriageway field prescribes 1.65% there. The height is
zero on both. The solver blends the two bank gradients over the original
arc after verifying the Section, bank side and native station. At chord
knots it uses the derivative from the side of the edge being evaluated.
The junction has 3,566 quads, a 3.252% maximum sampled physical slope
under the 10% limit, a maximum edge ratio of 7.955 and finite UVs.
All 9,515 candidate overlap pairs passed. This checks the physical owner
on a saved original network plan. See `nov19_325_complete105` and
`nov19_source_jets105_final`.

The saved ownership plan allowed independent tests of later junction surfaces
before the crossfall change at junction 91. Junction 122 shares one pavement
with junctions 125, 126, 137 and 821. The exterior covers about 4,999 m² before
subtracting two genuine islands of approximately 253 and 109 m². Its current
quad-chart construction exhausted its cell budget. This is a second Beijing
blocker. Both islands and all incoming profiles must survive a replacement
chart. The failed surface and its measured rim inventory are recorded in
`nov19_beijing_independent_surfaces01.json` and `nov19_beijing122_rims01`.

The source-bank partition exhausts its 50,000-cell allowance at block 16.
The principal-axis partition fits 23,067 pavement quads before grading, but
150 faces exceed the edge-ratio or corner-angle limit. A bounded NumPy
adjustment of the interior vertices, followed by 24 extra quads at six faces,
qualifies the planar pavement with 23,091 quads. Every outer and island rim
vertex retains its original coordinates. The stored maximum edge ratio is
7.984 and the minimum corner sine is 0.0423. Both grass disks qualify on
those exact same rims with 7,408 and 6,132 quads. The long tapered island
uses the existing grass rule of positive stored corners and edge ratios at
most 8. Its worst stored ratio is 7.990. The combined 36,631-quad planar
mesh passed the independent overlap and exact seam comparisons. Elevation and
final Blender publication still need qualification. Per-block evidence is in
`nov19_122_field_blocks01` and `nov19_122_principal_blocks01`. The verified
NumPy pavement mesh is in `nov19_122_numpy01`, and the complete planar
validation is in `nov19_122_qualified01`. The inspection artifact
`nov19_122_plan_inspection.png` shows the five-junction footprint and the
pavement quads at the rounded end of the tapered island.

The normal Beijing 122 surface build now reaches elevation. Section 54's
continuous profile has a maximum gradient of 5.93%. Two quads that interpolate
its prescribed heights would require a gradient of 21.60%. Subdividing them
within the original body lowers the bound to 8.60%. A subsequent long-edge
split to improve quad proportions raises it to 49.90% on two new faces.
The reproduction and measured per-pass bounds are recorded in
`nov19_beijing122_current02` and `nov19_122_body_stencil01.json`.
An alternative keeps the first subdivision and moves only its 34 new
vertices. It balances quad proportions and the slopes implied by the road
profile, then evaluates the exact nonlinear profile at the final positions.
On the measured 42-face patch, ten optimization iterations brought the
stored edge ratio below 7.984 while the required gradients stayed below
8.61%. All 24,221 original vertices kept their positions. The local evidence
is in `nov19_122_coupled_xy_guard01`. Passing the original body targets to
the drainage solver also allows it to preserve prescribed crown heights.
The maintained elevation recovery passed on all 23,125 pavement quads
in `nov19_122_body_replay04`. The normal junction surface builder then
completed on the saved Beijing ownership plan. Independent checks passed
for mesh validity, finite UVs and pavement overlap. This run took 246 s
for the single compound in the guarded test environment
(`nov19_beijing122_current04/report.json`). Full-city Blender publication
and interactive timing remain open.

### Acceptance status as of 17 November

| Required result | Current evidence | Status |
| --- | --- | --- |
| Original topology, returning Sections and explicit elevated exclusions | Original-OSM preparation audits and focused source/dependency regressions | Preparation checks pass; complete mesh gates below remain required |
| All supported Beijing ground components published | Earlier whole-city attempt: nine of ten components READY. The fresh candidate-35 composite-contact run accepts the complete 1,072-Section PLAN, all six requested repaired surface Tasks and all five final-owner inventories. The separate complete 7,698-Task attempt independently accepts PLAN but stops at Surface 9's quad-chart construction budget | Incomplete; repair that chart, finish every mesh/UV/overlap task, independent whole-component checks, maintained integration and fresh native publication |
| All supported Bergamo ground components published | Earlier whole-city attempt: thirty-two of thirty-five components READY. Candidate 35 completes all 972 tasks of the 133-Section component and independent checks of all 29 source-derived junction Products. The 1,255-Section mixed root passes local construction and independent checks; the fresh complete attempt stops later at a different exposed root. The 688-Section owned-boundary repair passes source/condition and planar-chart construction, but normal Product regeneration fails a physical slope constraint | Incomplete; resolve the newly exposed source and height failures, complete both large components, final owners and fresh native publication |
| Mesh quality, rounded analytical curbs, full profiles, UVs and material roles | Corrected source-labelled conditions pass independent audits of the affected local Products. The approved 0.07806 m² residual gap is filled with all other loops and profiles preserved | Local evidence passes; current whole-city stored-mesh, source-height and visual qualification remains required |
| Blocking/progressive equivalence, cancellation and publication cleanup | Candidate 34 passes nine native methods on actual OSM adaptation and the authored lifecycle matrix, including exact blocking/progressive parity, cancellation, staged source edits and publication-failure recovery | Original-city parity and callback timings remain pending |
| Existing street behavior retained | Candidate 34 passes 310 combined numerical methods and a fresh Karl-Marx-Allee native import. The next immutable candidate-35 native run accepts PLAN and actual Surface 175, then stops at Surface 191's finite height domain. All six accepted source-derived Products pass independent checks, as does the complete final-owner inventory | Not ready for maintained integration; resolve the remaining height-domain failure, complete native publication and repeat final combined checks |

These are acceptance gates rather than estimates of percentage complete. Final
qualification requires fresh imports of both original OSM files from the frozen
implementation. A failed-context replay, direct constructor diagnostic or local
surface audit must not be relabeled as a successful complete native import.

### Earlier complete-planning failures and their repair scope

The three preceding frozen attempts preserve their failure evidence. Those
attempts did not reach a completed PLAN or the requested surface Tasks; the
newer evidence below supersedes that result only within its stated scope.
Concurrent elapsed times measure diagnostic cost, not release performance:

- Beijing stops after 1,152.55 s at an original two-bank fit. The requested 2 m
  corner has a roughly 0.899 mm arc and a 50.5 nm sagittal displacement, below the
  unchanged 100 nm union welding precision. Its base operands are valid. A local
  diagnostic reconstructs the true arc between the two exact native contact
  vertices, retains every original and generated area operand, and passes
  continuous coverage and the unchanged complete-mouth/width/future guards.
  The newer typed native-contact/composite sampling construction preserves the
  same analytical curves and bounds every new chord by the original 0.25 mm
  approximation budget. Normal local meshing and the fresh combined PLAN now
  pass; complete-component and native publication remain required.
- Bergamo's 688-Section run stops after 2,243.79 s at another three-site group.
  Two retained circles would require an auxiliary radius strictly below
  0.415315 m, incompatible with the existing 0.5 m floor. This is not a sampling
  failure. The next construction must use a larger actual owned boundary chain;
  its genuine 35.5849 m² island is outside the approved tiny-gap fill policy.
- Bergamo's 1,255-Section run stops after 1,211.32 s at a mixed-profile singleton.
  Two original full-profile root-cap segments remain exposed. The unrelated
  Bundle half profile is genuine; it cannot be flattened to fit the existing
  full-only root contract. The previously qualified local Surface 43 is not
  reached by this attempt.

The fresh Beijing closure is recorded in `runtime35_composite_PLAN414/report.json`:
4,303 rebuilt Tasks, 76 groups, six actual owners (1014, 1229, 1441, 2042, 2089,
2383), 9,297 pavement quads and 323 exact incoming-profile samples. Actual final
inventories pass the original prepared-neighbor guard, including 447 owners
and 826 area operands for the composite construction. The PLAN takes 322.17 s
and the selected closure 422.81 s total during concurrent diagnostics. No old
Products, fingerprints, mesh, charts or Task builders are supplied. This is
complete planning plus selected surface validation, not a full-component mesh
or native-city import. The separate `runtime35_composite_FULL420` run starts
from empty state, independently accepts PLAN at 372.29 s, then stops at the
first surface, owner 9, after 406.17 s total. Its four genuine Bundle half-profile
inputs do not authorize an unbundled full-profile chart axis. The unchanged
quad-grading/cell limits remain enforced; no whole-component mesh is accepted.

Earlier candidate-34 PLAN/target acceptance retains its original scope and does
not establish acceptance of the stronger candidate-35 native-half closure.

### Additional local qualification, 17 November

The latest qualified full-source candidate is **34**, runtime `0.1.0.dev93`,
manifest SHA256
`3214e57f8cbebb7a09bc7ff3d2457ff8670d3d5c11167cd2d5a33808332a46cf`.
It includes the typed residual-gap rule, faithful retained-circle sampling,
equal-reread identities, finite free-pad and owned-prefix source conditions,
published-mouth identity binding, paired chord-jet domains, bounded fixed-rim
grass-cap recovery, and native-curvature root templates with staged recipient
validation. All 1,584 copied source files remain unchanged during qualification.
Historical results below retain their original scope; they are not evidence of
complete current-city acceptance.

### Integrated candidate 34

All **310** selected numerical methods pass, including the added root-template,
cache-epoch, paired-jet, physical-mouth, complete merged-Bundle, update/restore,
query, boundary and height checks. The Blender-only dimension-propagation module
is explicitly excluded from this numerical count. Nine native lifecycle methods
pass separately. A fresh native Karl-Marx-Allee import also passes: successful
finished operator, no errors, two READY components, all 254 supported Sections
covered and 353,443 total quads. Its reported 85.10 s is concurrent validation,
not a controlled performance comparison or proof of city-wide source heights.

The actual original 1,255-Section Bergamo network now passes a fresh Surface-1947
dependency closure on this exact source, including the complete PLAN: 5,211
actual Tasks, no restored products/digests, supplied geometry or replaced
builders. Its normal Product has 5,359 pavement quads, 42 exact eligible hard
body nodes and 51 exact original incoming profiles. Physical/stored slopes are
0.049691/0.049702 under 0.1; entrance slopes are 0.012222/0.012249 under 0.03,
with no interior minima. Field query error is below 1.4e-17.

The final generated-root owner guard is independently revalidated against the
**complete** eventual inventory: 665 original junctions, 97 groups covering 145
original sites, 617 physical owners and 904 area operands. Exactly one typed
generated root area is present, and its published guard matches this inventory.
This closes the actual target/PLAN gate, **not** the entire 1,255-Section mesh
backend or a native city import. The full 9,318-Task backend is checked separately.

The original 133-Section component completes all **972** numerical tasks from
empty state in a scoped candidate comprising full candidate 32 plus the exact
maintained field/binding/cap changes. Complete backend UV and overlap validation
pass; source and input hashes are unchanged. That source intentionally excludes
the later root-template integration, so it is not labelled a complete candidate-34
or native-city acceptance result. Its initial instrumented run had a native
access violation near a periodic stack dump; a separately labelled fatal-only
run completes in 697.56 s at about 244 MiB peak resident memory. The failed run
is retained and is not reported as a geometry refusal or memory-cap event.

Current Beijing contact and Bergamo-688 topology experiments remain artifact-only
until their generic derivation, complete eventual owners, normal Tasks and final
combined/native gates pass. In particular, successful local geometry or saved
grass repairs must not be substituted for a full composed Product result.

### Latest local gates: rounded Bergamo joins and Beijing root contact

The Bergamo-688 experiment now passes a complete composed constructor, rather
than a repaired saved grass mesh: **21,955 pavement and 4,480 grass quads**.
An independent audit passes on those actual buffers. Stored pavement ratio/sine
are 7.81304/0.05004, grass ratio/sine 7.94171/0.03411, maximum pavement slope
0.173269 under its original 0.176 limit, and both drainage fills are zero. All
119 incoming profile samples and 352 eligible original body nodes remain exact;
continuous full-width, logical-owner, rounded-tangent, UV, overlap and retained
query checks pass. Both genuine larger holes remain. This is a supplied-geometry
local constructor gate, not a completed PLAN or native publication.

The same three rounded-join families are also derived independently from ordinary
pre-repair source resources without selecting a junction, Section or extract ID.
The generic local preparation preserves the holes, claims, full profiles and
finite original contacts, and a second call returns the same prepared group.
Final neighbor ownership still requires the actual complete 688-Section PLAN;
the local reconciliation does not claim that inventory.

Beijing's later root refusal is a duplicate generated sampling station at an
exact retained-circle contact: the two outward numbers differ by 5.3e-15 m but
map to the same native cut and identical full-profile XY. A typed contact
rephase retains the true original circle station, every original unique point,
mouth plane and authored control. It does not claim that the two stations have
identical derivative callbacks, and the existing full-station-alias policy
remains unchanged. A separate failed-candidate dispatch requires both original
base unions and the complete typed native-root reference to validate before a
new root candidate is attempted; malformed bases continue to fail.

The unchanged ordinary source/admission builder and normal Product then pass
on full candidate 34 with these artifact helpers: **1,537 pavement quads**,
51 exact incoming profile samples, complete 2 mm source-width proof, stored
slope 0.17077 under its inherited 0.176 limit, entrance slope 0.00948 under 0.03,
and no minima or drainage fill. The body ledger explicitly has zero eligible
hard body nodes because all four owned native ranges lie within original corner
influence. No body-core preservation witness is claimed for that example.
The fresh current-source PLAN passes that root repair, then refuses a later
source union. Two operands sample the same original native bank on different
station lattices. A local common-lattice test restores a valid union while
retaining all authored controls and true contacts and tightening both analytical
error bounds. The complete sampler/coordinate implementation, normal Product
and fresh PLAN remain required before maintained integration.

### Independent complete-component source-surface audit

The completed scoped 133-Section run is independently inspected across all 29
source-derived junction Products, rather than only the difficult named targets.
They contain 51,037 pavement quads, 5,656 grass quads and 1,683 original incoming
profile samples; the maximum checked original curb-height error is zero.
The physical/stored quad, original curb-height, incoming profile, drainage,
entrance slope, UV, local overlap and retained-query checks pass throughout.
The audit distinguishes an unchanged analytical boundary-jet field from an
adjusted bilinear quad field: nodal values remain exact in both, while the
unadjusted field is checked against independently reconstructed analytical
conditions. Its interpolation difference from a quad-centre mean is reported,
not confused with a failed adjusted-mesh query.

The owner-specific full-width gate exposes one exception: an ordinary junction on a
Bundle member retains two native half profiles, which the existing ordinary
mouth-clearance eligibility excludes. A narrow incoming original source strip
extends up to 0.372 m outside its declared junction owner. Most of this strip
is covered by the adjacent Bundle pavement; the sixteen sampled points strictly
outside all published pavement are at most **0.637 mm** from it, within the
existing 2 mm approximation band. This establishes a local owner/source-range
contract failure, not a 0.372 m whole-network hole. Attachment-query and dependency
consequences are being checked separately.

The attachment-height consequence is subsequently reproduced with valid
hypothetical placements through the actual `FurnitureProvider.build` API.
Among 4,043 probes on the old component, 74 lie strictly on neighboring Bundle
pavement outside their declared junction owner. Every one receives datum height;
the largest mismatch against the actual published quad is **52.73 mm**. The
same probes on fresh candidate 35 have zero such foreign-pavement height
mismatches. These are authored test placements, not existing user furniture.
This establishes a height-query ownership defect; a stale-cache failure is not
claimed without an edit/rebuild test.

A local diagnostic of the existing clearance algorithm restores the complete
width by advancing one physical mouth from 0.396609 to 0.820610 m. Every original
corner, analytical curve, radius, Path and profile remains the same object.
The normal constructor and independent audit pass on the corrected group:
636 pavement quads, 51 exact profile samples, stored ratio/sine 7.45718/0.24376,
slope 0.161174 under the inherited 0.176 limit, entrance slope 0.004083,
and no minima or drainage fill. This is supplied local geometry, not a fresh
Task. Eligibility is being restricted to profiles authenticated through the
actual original Bundle membership; merely supplying a half-profile value must
not authorize the extension. Combined current-source and complete-component
qualification remain required after that change.

The extension subsequently passes nine portable eligibility checks, including
reversed/translated native Bundle roles and forged profile/width/ownership
refusal. Full-source candidate **35** is an immutable copy of all 1,584
candidate-34 files with only the self-contained clearance module changed; its
manifest SHA256 is
`287195ab4808e52e82adbc0552ec5b47ff539fff1f7b5519efab544cabcad5cc`.
The ordinary-mouth policy is explicitly versioned
`12-authenticated-native-bundle-half-profiles`.

From empty state, that candidate completes all **972** actual component Tasks,
normal UV mapping, whole-update overlap validation and numerical commit, without
restored products, supplied geometry or replaced builders. All 238,655 quads
pass the backend overlap stage. Source and input hashes remain unchanged.
The 900.81 s run is concurrent functional qualification, not a controlled
performance comparison. Independent checks then pass on **all 29** actual
source-derived junction Products: 51,063 pavement quads, 5,656 grass quads,
1,683 exact profile samples, zero original curb-height error, complete local
source widths, physical/stored mesh checks, drainage, UVs, overlap and retained
queries. This closes that numerical component gate, not native city publication
or final combined-source regression qualification.

Static local Blender inspection also covers the accepted rounded Bergamo-688
buffers at physical scale and tenfold height exaggeration. The analytical
silhouette is rounded, while exaggerated curvature lighting exposes residual
height rippling near fixed grass rims. This observation is retained separately
from the passing numerical limits; it does not constitute full native-city
visual acceptance.

The first fresh 688-Section PLAN with generic rounded-join discovery refuses a
different narrow retained-circle/auxiliary switch after 322.42 s. That unrelated
group passes its existing full-width proof but has a real sharp change of
boundary tangent. The local three-family success is therefore not a proof that
every discovered switch admits the same finite fillet. Neither a relaxed radius
floor nor silently skipping a real sharp switch is accepted. An additive,
original-bank chain construction is under review, retaining the original source
operands and circle records.

### Superseded full candidate-34 large-component run

The fresh 1,255-Section numerical run is stopped externally after 4,457.88 s.
Its complete PLAN and the actual first surface return normally, but Surface 43
has not finished after approximately 3,954 s in that builder. The worker
accumulates 4,386.84 CPU seconds across the whole run, with low observed memory
use. This is a substantial atomic-task latency, not a successful mesh timing.

Before stopping, an independent read-only check of the exact saved PLAN proves
that this ordinary junction's prepared owner boundary clips the complete original
Section-18 width. This is a local ownership result; no whole-network gap is
inferred from it. Two original approaches have real Bundle half profiles,
so candidate 34 skips the same clearance predicate identified in the 133-Section
component. The authenticated policy-12 extension enters that predicate and
correctly refuses the unchanged boundary after its finite existing correction
family is exhausted. That establishes changed admission, not a repaired Surface
43. Its original native source-union preparation is the next construction gate.

The stop record is explicitly `EXTERNALLY_STOPPED_SUPERSEDED`, with source/input
hashes unchanged and the validated 5,209-product PLAN prefix, task definitions,
resource observations and logs retained. It is not labelled a builder exception,
native crash, memory-cap event or GUI-cancellation test. The separate fresh
Surface-1947 qualification and complete final-owner proof remain valid for
their declared candidate-34 scope.

### Current local mesh and regression gates

The long Bergamo junction's default radial chart exceeds the unchanged 50,000
cell cap: its next layer would request 50,127 cells. The limiting 4.325 mm
boundary interval joins an original native source endpoint to a true prepared
circle contact. Both are compulsory; removing an incidental sample cannot solve
this case. A different computational partition uses only the genuinely
unbundled original Section's finite owned prefix. The two real Bundle-member
inputs remain present with their actual profiles and dependencies.

That partition regenerates a normal local Product with **19,461 pavement quads**
and no grass. Independent checks pass: 51 exact incoming profile samples, 32
exact original body targets, zero original boundary-height error, stored slope
0.098895 below its inherited 0.100572 limit, entrance slope 0.003688 below 0.03,
and no minima or drainage fill. The full normal constructor takes 570.53 s in
concurrent validation; this is still substantial runtime. Its independent audit
does not rebuild the mesh. A generic mixed-native input/signature contract
subsequently passes seventeen portable checks and actual signatures for both
the long junction and the new Beijing case. The original 1,255-Section backend
is now being evaluated from empty state with that contract and the existing
typed original-mouth seed policy. No local chart or mesh is supplied to it.

The later mixed-root singleton also regenerates normally with **7,559 pavement
quads** and no grass. Its exposed zero-height caps belong to two genuine full
profiles; the third input retains its actual Bundle half profile. The typed
owned-cap chain uses only original source bounds and a continuously certified
generated C2 curve. Independent checks pass: 51 exact incoming profiles, 22
exact hard body targets, zero original boundary-height error, stored ratio/sine
7.95731/0.031718, slope 0.067889 below 0.1, entrance slope 0.000922 below 0.03,
no minima or drainage fill, and UV/overlap/height-query checks. Eleven saved-scope
authority tests and eight separately authored portable construction/negative
methods pass. Evidence is in
`source35_mixed_prefix1255_full_backend01/normal_owned_cap_chain_product02/report.json`
and `actual190_owned_cap_chain02_independent_report.json` under the task artifact
root. A new empty-state complete run uses an immutable seven-file helper bundle;
neither this local Product nor its saved context is injected into that run.
That run subsequently stops during PLAN after 400.98 s at another singleton,
site 210. A 0.253 mm exposed full-profile root interval remains beside an actual
native Bundle half. Tighter observation converges to a positive analytical
gap, so the failure is not waived as a sampling artifact. The accepted local
7,559-quad mesh remains valid; complete 1,255-Section planning is not accepted.

For the separate 688-Section repair, normal regenerated Product 360 constructs
its planar chart after 755.52 s and rederives boundary conditions in 21.28 s.
It then fails physical entrance/body reconciliation: slope 0.226441 exceeds the
unchanged 0.175999996 limit. The 850.93 s attempt is preserved as a failure.
Neither its accepted planar chart nor the earlier source-jet checks establish
a qualified height field or Product. This diagnostic cost remains substantial.

The newly discovered two-owner Bergamo corner group also passes generic local
preparation, normal construction and independent checks. It uses the preferred
R2 native-chain infill and a mirrored R0.64 circle/auxiliary fillet. Two physical
profile planes move by 8.465 m and 9.605 m along their same original Paths; those
changes are explicit, and the old planes are not claimed unchanged. All 85
profiles at the new planes and 49 eligible body targets remain exact. Its normal
mesh has **2,851 quads**, stored slope 0.171648 under 0.176 and zero drainage
fill. The affected streets are ordinary full-profile roads, not Bundle members.
The fresh 688-Section run must still establish complete final owners, all affected
actual surface Tasks and the corresponding road trims.

Beijing's later mixed-profile root uses a distinct authenticated constant-crest
contract. The donor follows consecutive visits of the same actual Bundle role;
its real curb height is 0.08 m with zero gradient. No half profile is flattened,
and the old full-profile root reference continues to reject half profiles.
The normal local mesh has **2,429 quads**; independent checks pass on 68 exact
incoming profile samples, every generated crest edge, complete source widths,
stored slope 0.093289 under 0.1, entrance slope 0.002402 and zero drainage fill.
There are zero eligible body-core nodes in that case, stated explicitly. Fresh
candidate-35 preparation and the affected actual surface Tasks remain required.

Finally, the fresh native Karl-Marx-Allee test on candidate 35 fails at one
four-arm junction during cyclic-root preparation; only the other component is
published. The failure is reproduced from the genuine native failure inputs.
A diagnostic candidate-only rejection correctly propagates a subsequent invalid
fixed-base union, rather than accepting it. A separate actual-native-port
reference preserves the full zero-height donor and canonical half profile on
an unrelated approach; its reference checks pass, while contact refinement then
exposes an exact duplicate generated station at a retained original circle.
The separate mixed-native reference subsequently passes the unchanged
native-curvature construction after an authenticated retained-circle rephase.
The whole local singleton then completes normal construction, admission and
continuous source-width checks. A diagnostic signature stage initially lacks
upstream digests; the follow-up restores only the genuine native-input digests
after verifying all 253 Path and 253 profile resources are unchanged. Its normal
mesh and independent audit pass: **6,636 quads**, 68 exact incoming profiles,
31 exact body-core targets, zero source-curb height error, stored edge ratio
7.34896, minimum corner sine 0.244116, slope 0.064453 under 0.1, entrance slope
0.001465, no interior minima or drainage fill, and valid UVs, overlap and queries.
Five authored native-role/zero-donor authentication methods and all three
existing full-profile root methods also pass under the separate integration.

A shared structured incidence exception now carries the rejected union operands
and graph while preserving the existing union arithmetic, epsilon, budgets and
failure predicate. The mixed-root diagnostic uses it only inside the legacy
synthetic-strip proof scope; admission still independently validates the actual
fixed operands. No public error-message-based recovery policy is promoted.

The subsequent fresh native import gets past this repaired root but stops at
another exposed internal source seam. It remains partial, with one component
published, after 156.49 s of concurrent diagnostic execution. The captured
506-task update is Path/Profile source preparation, not a complete mesh PLAN.
The new seam is at a genuine consecutive Bundle width change: 2.5 m and 5 m
half profiles share the 0.08 m inner plateau, and only their positive-height
root halves are exposed. The third street's root is already internal. It cannot
be admitted by inventing a non-donor zero-height root for the preceding repair.
The existing inner-plateau width-continuation primitive is geometrically
applicable, but its first local exterior diagnostic clips the old narrow-street
mouth. Its proposed contact is at 4.625 m, beyond that old 4.35906 m mouth.
The separate exterior authority now retains the genuine dimensional primitive,
reconstructs the normal ports and passes the unchanged narrow-recipient proof.
The 160 mouth advances from 4.359056 to 4.75 m: its actual old exposed profile
measure is 0.499992 and the new measure is 1.0. This is recorded movement, not
a station alias or an unexplained exemption.

Normal local generation and independent checks pass with **1,527 quads**, 51
exact incoming profiles, stored ratio/sine 7.49227/0.299351, slope 0.065865 under
0.1, entrance slope 0.005778 and no interior minima or drainage fill. A separate
semantic audit checks 116 generated curb edges and 580 rule samples: the actual
0.08 m plateau and its zero gradient are exact. Both original outer circles
retain exact centers, contact stations, radii and sweeps. Eight authored
construction/authority methods, including reversal and rigid transformation,
and seven existing hole-only/full-root methods pass. All new exterior areas
participate in the ordinary final prepared-owner exclusion checks.

The next fresh native run applies this general width rule at both the original
case and another matching junction. Its complete main-component PLAN is accepted.
An independent read of the genuine validated native prefix revalidates all
eventual owners: **75 original junctions, 11 groups, 70 physical owners and 113
area operands**. Published generated-area guards agree with this full inventory.
The import remains partial after 200.90 s of concurrent diagnostic execution:
normal surface 175 refuses a generated bank arc in the old original-corner-only
height-transition authority. Its two real native banks have heights 0 and 0.08 m;
one fitted curve traversal runs opposite to its actual junction End. A native
bank-origin authority passes the actual native-incidence checks. A separate
finite crossing-road corridor also proves the complete selected inner-bank
interval lies within another original incoming road, with explicit analytical
approximation margins and no intersection with the unchanged hard body core.
It does not substitute the mouth for a missing opposite-bank contact. Normal
construction then exposes a second, direct native-bank height conflict elsewhere
on the same loop: a 0.08 m inner bank meets a zero-height full-profile bank.
At that checkpoint the conflict remains rejected and no successful surface is
claimed. Geometry, profile and slope limits are not relaxed.
The subsequent separate exposed-sector construction now addresses that second
switch. It authenticates the actual neighboring original banks and both real
intervening STYLE sector operands, then uses their 3 m preferred radius. The
intermediate incoming road remains an incoming road; it is not claimed or
consumed. The old fully-claimed-sector rule is unchanged. Normal source-bank
admission preserves all profile planes, original source vertices and circles,
and independently revalidates all eleven final prepared-neighbor groups.

Two distinct finite inner-bank transitions use complete continuous overlap
witnesses in original incoming-road polygons, with explicit native and foreign
chord-error margins. Their native intervals are 8.91497–9.97415 m and
0.340852–1.76534 m, both outside the unchanged hard body core. The second witness
belongs to a third original incoming road; it is not the generated fillet or
an accidental reuse of the first fillet's cached authority. Domain caching is
bound to the exact generated-fragment owners. Sector references additionally
bind the original graph UUID, input/settings fingerprint, source identities and
unchanged pre-construction group.

The local normal Product and independent audit pass with **9,310 pavement quads**,
68 exact incoming-profile samples, 34 exact original hard body targets, zero
source-boundary height error, stored ratio/sine 7.86554/0.093917, slope 0.0996613
below 0.1, entrance slope 0.006182 below 0.03, and zero minima or drainage fill.
UVs, overlaps, height queries and geometry-release behavior also pass. Normal
generation takes 64.40 s and the independent audit 23.27 s during concurrent
diagnostics. Twelve saved-scope affinity/authority tests pass, including foreign
UUIDs, changed settings, borrowed sectors and rejection of self/generated-area
overlap witnesses. Seventeen existing authored full-root, hole-width and finite-
bank methods also pass with the complete helper pipeline. Evidence is in `kma175_occluded_banks12/report.json`,
`kma175_occluded_banks12_independent_report.json` and
`kma175_authority13_report.json` under the task artifact root. This explicitly
reconstructs original source operands from the genuine completed PLAN; it is
local construction evidence, not a fresh PLAN, actual Task or native import.
The next fresh native run is `kma35_mixed_zero_donor_native03`, with sixteen
immutable helpers (manifest SHA256
`e33a30e6b8eb238e204dd321b4c38b1f51f096f37017ae59bbbdb8f15e409890`).
It supplies no Products, fingerprints, geometry or replacement Task builders;
the earlier completed-PLAN context is used only for the separate local evidence.
This run ends after 299.23 s with one component published and a main-component
failure at actual Surface 191: its existing finite inner-bank domain cannot fit
the required height transition. Surface 175 is now genuinely accepted by the
normal native Task, with 7,155 quads; its incidental computational sampling
differs from the separately reconstructed local diagnostic. Independent checks
of all six accepted source-derived native Products pass: owners 19, 61, 79, 141,
151 and 175, totaling 22,858 pavement quads and 493 exact profile samples. The
accepted native PLAN's complete final-owner inventory also passes: 75 original
junctions, 11 groups, 70 physical owners, 115 area operands and 864 validated
Tasks. Evidence is `kma03_native_prefix_source_products19_independent_report.json`
and `kma191_native_failure18/final_owner_inventory_report.json`. These are real
accepted native Tasks in a failed overall import, not a whole-import pass.

The new height-domain failure retains a real 0.244567 m opposite-bank contact
bound, but needs 1.63282 m of transition length. A separate read-only feasibility
proof finds an existing 1.69603 m source interval with continuous overlap in a
third original incoming road, before the real same-bank contact and fixed mouth,
with no hard-body overlap. The old opposite-bank domain remains unchanged.
Typed admission and normal Product construction for this stronger independent
overlap domain are still required; this feasibility result is not a repaired
Surface 191.
Complete native publication, all final owners and actual target-surface Tasks
still need completion for the final integrated candidate. The accepted PLAN
and independent local meshes do not establish a complete native mesh import.
The successful candidate-34 regression and candidate-35 component results must
not hide this native regression.

Candidate 30's fresh Karl-Marx-Allee native regression imports both ground
components successfully and records complete supported coverage with no errors.
Its 116.81 s import was concurrent with other numerical checks and is **not** a
controlled performance measurement. The independent whole-scene source-height
audit was not part of this native regression. Its 216 selected numerical tests
pass with all source hashes unchanged.

Candidate 31 subsequently passes all **248** selected numerical tests, including
the new equal-reread source identities, owned-prefix authority and coordinate
ledgers, free-pad source conditions, and faithful retained-circle sampler.
There are no failures or errors and every frozen source hash remains unchanged.
The Blender-only dimension-propagation module is excluded explicitly rather
than counted as passing in this numerical run. A fresh native regression against
this same candidate also passes: the operator finishes successfully, both ground
components are READY, supported coverage is COMPLETE, and the error list is
empty. The reported 154.20 s import is concurrent validation, not a controlled
timing comparison. Every frozen source hash remains unchanged. Whole-scene
independent source-height inspection and Beijing/Bergamo qualification are not
implied by this regression.

Nine native candidate-31 lifecycle methods also pass: exact blocking/progressive
OSM-adapter parity, cancellation before preparation, hidden/staged cleanup,
preservation of committed results and user edits, prepublication mesh failure,
staged source changes, failure after numerical commit, and the captured material
override. Source hashes remain unchanged. These authored lifecycle checks do not
qualify actual large-city callback latency or cancellation during a long chart.

The fresh candidate-30 133-Section Bergamo coordinator advances beyond the
earlier boundary-fit failure at surface 3107, then fails at ordinary surface
3130. The cause is an unrecognized, certified degree-six finite free-pad
envelope, rather than an invalid source radius or insufficient source length.
The new [free-pad source rule](STREETS_FINITE_FREE_PADS.md) validates its actual
producer, complete pad ownership, original through profiles, endpoints,
tangents and bounded curvature. A fresh local preparation and unsupplied normal
Product now pass with 299 quads. Its independent stored-mesh audit passes:
maximum slope 0.06586 under 0.1, entrance slope 0.00183 under 0.03, 34 exact
incoming profile samples, no detected minima or drainage fill, complete source
ranges, UVs, overlap checks and query/release behavior. Fourteen focused tests
pass. This is a local constructor result, not a complete component or native
import. Candidate 31's fresh complete-component attempt subsequently passes
surface 3130 and advances to **953 validated tasks** before a different failure
at Bundle plan 10. The visible-port registry's index addresses the earlier
projected layout; the accepted junction Product has a different physical
approach inventory. Reusing that positional index selects an unrelated street.
An explicit identity-based physical-mouth binding passes this exact Bundle plan
as a diagnostic algorithm variant, with no geometry, source or radius changes.
The maintained correction and its authored hydration regressions are qualified
separately before another fresh complete-component run.

The maintained [published-source binding](STREETS_PUBLISHED_SOURCE_PORTS.md)
passes the exact local normal Bundle-plan builder in candidate 32 without a
replacement builder or supplied geometry. All eight new binding methods and
fifteen related methods pass. One additional older merged-Bundle method fails
in two subcases; both identical height/field failures also reproduce on the
unchanged candidate-31 baseline. Declared paired chord jets and a bounded
direct-quad grass-cap patch subsequently fix those causes; a scoped numerical
candidate passes all 43 selected field/binding/complete-Bundle/height methods.
These changes still need combined current-source qualification. The local builder uses
genuine restored upstream products and therefore does not qualify a fresh
candidate-32 coordinator or complete native import.

The Beijing owned-prefix case likewise has a passing local normal Product and
independent audit: 7,225 pavement quads, 102 exact incoming profile samples,
nine exact genuine hard body nodes, stored slope 0.06624 under 0.1, entrance
slope 0.00440 under 0.03, and no detected minima or drainage fill. Its fresh
candidate-31 main-network PLAN passes the previous logical/physical contact
authority failure, then rejects a different bank-to-clipped-corner shortcut:
part of the proposed arc is hidden by an original retained R2 corner. The next
construction must meet that actual original circle and preserve its retained
tail; reducing the exposure requirement or accepting the hidden arc would
violate the boundary contract. No complete-city or progressive publication
claim follows from these local results.

The 1,255-Section Bergamo component now reaches a distinct acute source-root
case at an ordinary junction. Both original contacts lie within their finite
source paths and complete-width containment passes. The actual refusal is an
artificial root cross-section remaining on the outer pavement boundary, not
an approach-length shortage. Bounded alternative curb constructions remain
artifact-only until their complete exposure and original-source proofs pass.

The retained-parent rejection in the 688-Section Bergamo component exposed a
distinction between a corner's logical radius authority and its physical
contacts. The original node has ends `164/0 → 163/0`; the valid physical corner
contacts `164/0 → 165/0` through owned Section 163, preserving Section 165's
actual 5 m width. A new typed validator independently reconstructs that original
source walk, verifies finite ownership, the unchanged original STYLE radius
sequence and the actual circular tangencies. The saved real failure passes
this authority check without altering its geometry or source. Thirteen focused
methods also pass, including reversed storage, false-source/ownership refusals,
prepared-neighbor proof publication and the new parent-chain dependency epoch.
See [original contact authority](STREETS_SOURCE_CONTACT_AUTHORITY.md).

Separately, the 1255-Section component's failed exact corner-exit handoff now
has a passing normal Product **algorithm-variant** diagnostic: 2,865 quads,
26 exact hard body nodes, 51 exact incoming profile samples, physical/stored
slopes below 0.053 under its 0.1 limit, joins below 0.0044, and no detected
interior minima. This preserves the original circular geometry and uses the
existing complete publication checks. Promotion and a coherent current-source
task/component run remain required. Neither local result is a complete native
import of Bergamo or a progressive-experience acceptance result.

The combined immutable candidate 28 passes all 148 selected numerical tests
with zero failures/errors and every one of its 1,548 source hashes unchanged.
These include ordinary ownership/update projection, radii, continuations,
returning roads, source-height Products, source placement, Bundle topology,
the adapter, normalized source rules, displaced authority and bank
reconciliation. The Blender-only dimension propagation module is explicitly
excluded from this numerical run; it is not counted as passing here.

Normal regeneration of the 688-component ten-site area from its genuine
upstream products passes the corrected original-parent authority check. It
then fails the unchanged finite retained-circle fit constraint after 142.12 s.
The first problematic inner loop has area 0.07806 m² and an approximately
0.299 × 0.583 m bounding box. Its three original circular boundaries have
radii 0.524288, 3 and 0.98304 m. None of their three switch pairs admits an
existing allowed finite fillet at or above the 0.5 m derived-radius floor.
The larger neighboring pair does admit existing candidates. At that checkpoint
the radius floor and holes remained unchanged; the result was not a passed
Product.

The user subsequently chose to fill tiny residual generated gaps with pavement.
The new [certified residual-gap policy](STREETS_RESIDUAL_ROAD_GAPS.md) accepts
only small convex clockwise loops bounded wholly by validated original STYLE
corner arcs. Native banks, incoming profiles, half-profile plateaus, participating
Bundle Sections and explicit radii cannot authorize a fill. Conservative
analytical common-disk bounds establish the small size independently of chord
sampling. Eight portable tests pass; a ninth update-policy test also passes
after the owned-Bundle normal-chart epoch is wired.

The real ten-site group's direct preparation proof passes in 43.73 s: exactly
the approved 0.07806 m² gap is filled, the other three holes and every surviving
loop remain identical, incoming profiles and original inputs are exact, and
fresh continuous complete-width cells pass the unchanged 2 mm limit. Independent
physical union reproduces the resulting topology with zero area discrepancy.
Normal current-source preparation and Product generation are being checked
separately; this local proof is not a complete native import.

The separate 133-Section component now passes a fresh empty-coordinator PLAN
with retained-arc policy 3: all 547 genuine preparation tasks rebuild, producing
nine physical groups in 94.73 s. No products, digests, repaired groups or meshes
are restored for that run. The source is immutable and unchanged; the largest
group preserves its three original holes and four complete retained fillets.
Mesh tasks and full native/progressive acceptance remain separate gates.

The same immutable 133-component implementation subsequently runs 720 genuine
validated preparation/mesh tasks before failing at ordinary surface 3107's
boundary-jet polynomial fit. Its normal source-height surfaces 447, 498 and 995
all complete. Their independent source/body/stored audits are separate from
the constructor results. The full component still fails overall until the
ordinary surface and final native/progressive gates pass.

The ordinary 3107 failure is localized to two explicit SourceOffsetRule pieces
at native station zero, including the original negative auxiliary lead. Exact
declared zero anchors and the analytical auxiliary-frame covector close the
actual BoundaryData/Field gate with fit error 1.06e-10, without changing the
source curves or station map. The next full current-source coordinator must
rebuild this task; the earlier 720-product prefix remains diagnostic evidence.

Independent audits now pass the actual 447, 498 and 995 Products' original
curb heights, incoming profiles, stored mesh/slopes, drainage, UVs, grass seams
and overlaps. The 995 core has 1,156 exact hard crown rows and 2,338 source-width
cells. Read-only reconstruction of the other two groups' 21 source intervals
finds no hard core remaining outside genuine original-junction influence under
the existing whole-area policy. Within that influence, original peak values
may differ by 13.060/12.353 mm; this is recorded rather than described as exact
body preservation. No new runtime body marker is attached to those groups.

The fresh candidate-31 coordinator's actually completed 447 and 498 Products
also pass independent checks without modifying their buffers: 4,566 and 4,182
pavement quads, respectively, 68 and 85 exact original incoming profile samples,
and zero source-boundary height error. Continuous source-width certificates
cover 1,070 and 992 cells at the unchanged 2 mm limit. Stored slopes are 0.09314
and 0.08156 under 0.1; entrance slopes stay below 0.00457 under 0.03, with no
interior minima or drainage fill. Quad quality, grass seams, UVs, overlap and
queries pass. These are validated products from the still-unfinished component
run, not evidence that its later tasks or complete native import have passed.

Its actual 995 Product also passes independent physical/stored checks with
15,895 pavement and 5,323 grass quads, 68 exact incoming profile samples, zero
curb-height error, 1,152 exact eligible hard body nodes and 2,338 continuous
source-width cells. Stored slope is 0.09044 under 0.1. The original audit driver
then encounters an unsupported ordinary-layout assumption at the separate 3130
Product; that is a test-driver error, not a failed 995 check. After correcting
the ordinary audit dispatch, the actual coordinator-built 299-quad 3130 Product
passes independently as well. The overall coordinator still fails later at the
Bundle-mouth binding gate described above.

The separate Beijing 691 normal preparation and unsupplied Product now pass
with 5,941 pavement and 142 grass quads using certified finite native Bundle-role
normals after the earlier chart families fail. Its independent audit checks
119 exact incoming XY/height samples, stored slope 0.08635 under 0.1, entrance
slope 0.01423 under 0.03, no drainage fill/minima, finite UVs, overlaps, original
276-cell width coverage and query/release behavior. Its reconstructed support
has zero hard body nodes; no nonzero body-target claim is made. Full native
runtime and prepared-neighbor qualification are still required.

### Latest bounded diagnostic checkpoint

The 133-Section Bergamo component completes its historical numerical backend,
including UV generation. Independent all-part overlap checks pass in physical
and stored coordinates, and the 48 junctions stay inside their elevation
envelopes. However, an additional source-height audit finds 218 unambiguous
ordinary zero-curb edges raised to crown height in owners 447, 498 and 995.
Its source-curb heights therefore remain unqualified. The numerical completion
is not a fresh whole-city native import or proof of the full height contract.

The strict height-feasibility correction at Bergamo surface 295 is integrated.
Its matching immutable-source coordinator replay passes in 31.15 seconds and the
independent product audit passes. Twenty-four height/projection tests pass.
The corrected continuation next failed mesh grading at the separate ten-site
surface 318. A selective computational cut-endpoint policy now passes its full
normal Product constructor and independent audit: 19,677 pavement and 1,588
grass quads, all three holes, exact original canonical rim points and profiles,
complete source-width proof, stored slopes, drainage, seams, UVs and overlap.
The production coordinator replay passes the actual surface task in 914.87
seconds from a matching immutable source variant. Its independent task-Product
audit passes, and XYZ/quads are byte-exact to the earlier constructor result.
The additional source-height audit also finds wrongly raised ordinary hole
curbs on this product and on surface 483. The cut policy is planar/source-width
qualified, while the shared boundary-condition correction and complete-component
and fresh native acceptance remain pending.

Beijing's five-site coupled owner has a validated source plan. Owned-prefix
normal cuts, existing physical vertices and a bounded direct separator-tip cap
produce a local combined diagnostic with 19,442 pavement and 2,566 grass quads.
Its normal Product constructor and independent planar/source-width, stored-mesh,
slope, drainage, seam, UV, query and overlap checks pass. Source-labelled boundary
conditions retain the ordinary hole at zero and the genuine median at its
original plateau. Physical coordinates and queries remain identical at display
exaggerations zero, one and ten.

This is an opt-in diagnostic, not an admitted production chart or a completed
normal coordinator task. Contact-free body checks identified approximately
6–11 mm of crown-height drift on some consumed source intervals. A subsequent
artifact-only reconciliation restores 8,584 existing body vertices to their
actual original or certified paired profile targets, retaining the same XY,
19,442 quads and physical boundary. Only real boundary nodes remain drainage
outlets. The existing projection, height slabs and drainage machinery passes in
2.26 seconds; no solver or validation limit was changed. Complete-row overlap
checks establish the body domains instead of assuming that a scalar contact
station ends all junction influence.

Independent source queries confirm that the earlier large height errors are
removed. The coarsest ordinary-road cell still has a measured 2.426 mm bilinear
interpolation deficit, with all four source-height vertices exact; this is a
resolution-dependent mesh approximation, not an exact continuous crown.

The maintained generic source-body helper and explicit
`CompoundProduct.build(..., source_heights=True)` path now reproduce that
corrected physical mesh exactly on the accepted supplied planar family. The
normal constructor, source conditions, body helper, slope/drainage/publication
guards and Product assembly run without replacement. It keeps all 8,584 exact
body targets, zero ordinary-block height, the original 0.0799999982 m plateau,
and matching physical rims and queries.

The later immutable 117a qualification constructs the chart and Product through
the normal maintained algorithms, without supplied meshes, cuts or conditions.
Its new general owned-prefix branch selects role0 after the older branches
fail: 19,046 pavement and 2,572 grass quads, 7,721 exact original body targets,
zero minima and 0.008895 entrance slope. Physical queries, UVs, pavement/grass
overlap and independent release checks pass. Constructor time is 339.65 seconds
including those older failed attempts, not an interactive-performance claim.
The actual coordinator Task and complete native city gates remain outstanding.
Runtime activation is now implemented in candidate 26, as described below; that
code change alone is not a successful runtime or city qualification. See [paired-prefix charts](STREETS_PAIRED_PREFIX_CHARTS.md)
and [source-body height guides](STREETS_SOURCE_BODY_HEIGHTS.md).

Bergamo surface 821 also has a planar-only diagnostic with 7,168 pavement and
3,136 grass quads. Exact source/rim/profile points, complete source-width proof,
stored quad quality, shared seams, UVs and overlap checks pass. The height and
publication path remains unqualified, and this recovery is not admitted into
normal production tasks.

Bergamo Bundle 9's source ledger and normal task pass. Separator 21 remains
incomplete: its earlier joint-gap prototype passes, but subsequent source
preparation evidence changes that construction's premise. Three retained
original Sections carry `placement=right_of:1`, previously ignored by Blosm.
The general constant-placement interpretation is integrated and passes 46
focused native/policy regressions. Both original extracts pass normal
preparation and original-source boundary support; all IDs, controls, widths and
topology remain exact to the prior capture. Only Bergamo Sections 78, 79 and 85
change their offset from zero to the specified 1.75 m; Beijing inputs are exact.
The Bundle must be reassessed from all corrected source offsets before adding
another body family. See [OSM way placement](STREETS_OSM_PLACEMENT.md).
The [OSM placement specification](https://wiki.openstreetmap.org/wiki/Proposed_features/placement)
defines the way's position within the carriageway; this investigation must retain
original nodes and direction rather than invent a geometric displacement.

The next acceptance gates remain qualified meshes for these difficult cases,
fresh complete imports of both original extracts, and progressive/blocking
parity, cancellation, publication cleanup and final visual/material checks.

### Current runtime candidate 26

`frozen_validation_source26_runtime86` is a coherent disposable copy of the
maintained repository: 1,530 files, 203,429,865 bytes, with identical source
hashes before and after copying. Both actual compound Surface Task paths now
request the original-source height policy. The public constructor still has
an explicit opt-in for older comparison tests; runtime generation does not use
that comparison default.

The candidate includes finite owned-connector contact reservations, with
separate logical corner and actual native contact identities. Their policy
and local geometric signature participate in ownership, Surface, layout and
transition-seam cache invalidation. Six projected-policy tests pass from this
exact copy, including all four new native-parameter/bank/arc/connector epochs.

The native source-height rules declare exact original curve-piece breaks.
The boundary field fits separate internal polynomial pieces while retaining
the original physical-edge projection and blend weights. At the difficult
Bergamo 351 boundary this produces 395 prediction pieces on the unchanged
392 physical edges, with maximum fit residual 9.31e-11 below the unchanged
1e-8 threshold. The degree-20 and 20,000-piece limits remain unchanged. This
is a boundary-field pass; the complete Product and actual runtime tasks have
separate gates.

Fresh checks from this candidate are the Beijing 122 actual Surface Task
dependency closure and the complete 1,255-Section Bergamo source PLAN. Neither
injects an old group, mesh, product or fingerprint. Both currently fail normal
PLAN preparation: Beijing at full bank-arc exposure and Bergamo at a repeated
retained-parent contact. The first is a near-tangent generated-endpoint
roundoff artifact; the second is a polygonal tangent alias inside an explicitly
consumed original parent interval. Bounded fixes have separate local gates;
neither failed fresh run is a runtime or city pass. The target driver saves
actual runtime failure context or completed
numerical products for diagnosis before releasing buffers, without marking a
failed or unfinished update as reusable or accepted.

The normal Bergamo 351 Product now passes from this exact full candidate in
239.12 seconds, with no supplied geometry, cuts, conditions or replacement
builders. The finite prefix admission rejects two infeasible charts and selects
the existing source-364-axis candidate: 1,236 pavement quads, 69 exact body
targets, 51 exact full-profile samples, no minima, stored edge ratio 7.62701 and
corner sine 0.05651. Physical/stored maximum slopes are 0.16025/0.16046 against
the actual 0.17600 limit; incoming longitudinal slopes are about 1.34%.
This is a normal constructor pass, not an actual city-level Surface Task pass.

A broader 130-method numerical regression attempt exposed three real failures
in existing continuation/returning fixtures after source-height activation.
They are retained as failures until their fixes pass their normal backend
tests. An additional loader error came from running a Blender-dependent module
outside Blender: that module was then run correctly in isolated factory Blender
and all eight dimension/placement methods pass. The fixed original dimensional
plane correction described in the
[ordinary-mouth contract](STREETS_ORDINARY_MOUTH_CLEARANCE.md) passes its normal
backend regression and thirteen focused tests. A complete fresh regression run
is still required after the boundary-rule fixes settle.

### Runtime candidate 27 and grass topology

The next coherent full copy, `frozen_validation_source27_runtime88`, contains
1,534 files (203,505,574 bytes), with no changed source hashes during copying.
It includes source-height epoch 2, faithful original lead/core pieces, original
evaluation-batch references and coupled retained-fillet radius reservations.
The broad numerical core passes all 129 actual methods; the eight
Blender-dependent dimension/placement methods have separate native evidence.
All 31 focused source-rule/returning/continuation methods pass. The full
67-method Role/body/jet attempt has no unrelated failures: one grass method
reports its three forward/reverse/transformed subcase errors.

Those grass errors are a local quad-ear constraint, not excessive physical
boundary slopes. The maintained direct-quad repair has a separate declared
source-variant qualification: 42 methods pass, including normal gapped-source
generation in all three configurations, exact physical boundaries, unchanged
slope/drainage limits and existing publication behavior. See
[fixed grass quad ears](STREETS_FIXED_GRASS_QUAD_EARS.md). This repair still needs
the next complete runtime qualification.

Fresh candidate-27 source PLANs do not yet qualify both cities. The 133-Section
Bergamo component stops at an original circle incorrectly assumed to occupy a
whole prepared curve; explicit circle-subinterval provenance is being
implemented. The 1,255-Section component advances through the earlier retained
fillet failures, then stops at a separate analytical bank-contact branch.
Both runs retain unchanged original inputs and source hashes and save their
actual fresh prefix for diagnosis. Neither is a complete mesh or native import
pass. The fresh 688-Section PLAN also stops, after 1,870.37 seconds, at a retained
circle parent's radius/incidence-authority check. Its original inputs and copied
source hashes remain unchanged. The actual prefix is saved for a normal local
diagnosis; it is not a qualified source PLAN or complete component.

The original Karl-Marx-Allee OSM file also passes a fresh full native candidate-27
import: both ground components commit, both report complete coverage, and the
native import reports success with zero errors. The copied source hashes stay
unchanged. Its saved scene and reports are under
`kma_native27_full_regression/` in the task artifact directory. This is native
coverage/regression evidence with the source-height runtime enabled, not an
independent whole-scene source-height audit, a progressive-UI responsiveness
qualification or a timing benchmark.

### Normal owned-connector Product qualification

Genuine Beijing owner 185 now rebuilds through the maintained source-cycle
preparer, with separate logical and actual native contact identities, sealed
finite contact ranges and original evaluation-batch references. Its normal
`CompoundProduct.build(..., source_heights=True)` and ordinary chart constructor
pass without a supplied mesh, chart, cuts, conditions or replacement builder:
28,625 pavement quads and 108 grass quads in 329.6 seconds, with about 321 MiB
peak resident memory. This is a local constructor qualification, not a complete
native runtime import or a performance target.

The independent stored-buffer audit also passes: maximum edge ratio 7.84808,
minimum corner sine 0.050135, physical/stored maximum slopes 0.128618/0.128992
below the source-derived 0.130370 limit, maximum incoming longitudinal slope
0.011716, no minima and no drainage fill. All 204 incoming-profile samples have
zero measured XY/height error. The original 3,280-cell source-width proof stays
within 2 mm; UVs, pavement/grass and self-overlap, shared queries and independent
geometry release checks pass. The remaining Beijing 691 geometry/height gate
and full city-level acceptance are separate requirements.

### Finite original-bank reconciliation

A local Beijing 1132/1156 diagnostic now qualifies two original-bank infills
with the existing 1 m STYLE radius and unchanged original profile planes.
Shared native source partitions retain every old prepared source vertex and
avoid false holes from independently sampled bank closures. The normal Product
and chart constructors pass with 2,600 pavement quads; independent stored
ratio/sine are 7.6871/0.04795, stored physical slope is 0.172819 below its
source-derived 0.176000 limit, longitudinal joins are 0.013001 and drainage
fill is zero. Incoming XY/heights, complete 604-cell source-width proof, UVs,
self-overlap and queries pass. The generic finite-bank API reproduces both
circles exactly and its four maintained tests pass.

This remains an explicit source-boundary diagnostic rather than automatic PLAN
admission or native-city qualification. See
[finite source-bank rounding](STREETS_SOURCE_BANK_ROUNDING.md) for orientation,
ownership, finite support and full-arc exposure requirements. Other source-curb
reconciliation cases, including Beijing 691 and 185, remain separate gates.

### Identical original-bank analytical contacts

The coherent cold 23 Bergamo source plan reaches a source-union contact failure
at Source 125. A focused actual-site reproduction identifies alternating short
ownership runs on the identical original Path and native bank: an analytical
OffsetCurb and a prepared chord representation. This is a singular coincident
contact, not a missing original connection or a proximity-based junction.

A conditional native-bank representation now evaluates only that proven bank
on its original source curve, retaining all original physical vertices, exact
Route/Exact knots, finite original parameter ranges and the original declared
source-chord tolerance. The existing canonicalizer removes the two redundant
ownership runs with zero endpoint discrepancy. Twenty-five focused regressions
pass; the actual local ordinary source-union constructor passes in about
4.5 seconds, retaining its complete original mouths and 278-cell, 2 mm source
proof. A fresh immutable cold24 plan advances to the later Section351 gate in
357.28 seconds, using all 5,209 genuine selected tasks and no injected products
or digests. It is not a passing complete plan. The newly exposed old profile cut
on Section364 is genuinely incomplete; a typed record of the existing outward
lead search now passes the local source-union and width checks. Its subsequent
mesh-grading gate remains separate. Complete plan, source-height and city-mesh
qualification remain pending. See
[source curb parameters](STREETS_SOURCE_CURB_PARAMETERS.md) and
[ordinary mouth clearance](STREETS_ORDINARY_MOUTH_CLEARANCE.md).

### Source-labelled curb heights: reopened qualification

The original `BoundaryConditions` policy assigned crown height to any port-free
hole. A compound can enclose an ordinary road block as well as a Bundle median.
The former's full-profile outer curbs must retain their original zero height;
the latter may retain an inner plateau and real cap transitions to lower curbs.
A loop being a hole does not establish which height policy applies.

Read-only source-provenance audits identify 436 wrongly raised full-crown edges
on the accepted 318 product, 136 on 483, and 218 across the completed historical
133-Section component. Beijing's first hole has ordinary zero curbs; its second
hole has a real median plateau but mixed-height caps. Exposed outer curb chains
also lose source heights when one long interpolation spans consumed roads.
Previously passing height-envelope, slope, drainage, incoming-profile and
overlap checks do not prove this missing boundary requirement.

The correction must derive conditions from actual strip sides, original
approach profiles and original corner/contact provenance, including explicit
source-span mappings on normalized carriers. Original incoming profiles,
analytical curbs, source IDs and real median plateaus remain authoritative.
No boundary may acquire a height merely because its loop is closed. All affected
products require source-labelled height/jet checks and normal height/publication
requalification after the shared fix; a planar mesh pass alone is insufficient.
The [source curb and grass-height document](STREETS_SOURCE_CURB_HEIGHTS.md)
describes the implementation, physical query behavior, formulas and tests.

The optional source-provenance implementation is not enabled globally. Some
actual compound exteriors contain non-tangent joins between exposed original
corner curves. For Beijing owner 185, one real-bank rounded replacement passes
its source/profile checks and exposes the full preferred 2 m circle, but five
other genuine arc-to-arc reflex joins remain. Exact analytical turns range from
about 28 to 99 degrees; they are not sampling phase errors. A nonzero inherited
crossfall cannot be assigned unchanged as a smooth zero-height jet at two
nonparallel curbs. A defined finite junction-only crossfall transition may
reconcile the jets while preserving the boundary heights, full incoming profiles
and contact-free source bodies. It does not round the XY boundary. Actual
geometric rounding and height transitions are separate requirements; solver
clamping or changing genuine median peaks would not satisfy them.

### Latest ordinary-source and regression checkpoint

Complete source-width checks now run on remaining ordinary junctions after
compound ownership has been established. They distinguish disjoint owned
prefixes from an unowned Section middle and retain the original source sampling
budget. Failed style defaults may use existing bounded radius fitting; explicit
user radii remain strict. Source-backed singleton recovery uses real original
mouths and source identities rather than fabricating a graph cycle.

The bounded three-arm root recovery passes ten focused tests, an additional
strict/untyped-context check, and a byte-exact comparison with the earlier
default implementation. Its fixed-mouth variant retains all original mouths;
the separately tested advanced-mouth variant is not admitted. These are local
implementation checks, not a successful complete source plan.

Both original preparation and actual Blender source-object transport checks pass.
Candidate 21's Bergamo 1,255-Section cold source plan progresses past the previous
ordinary failures and stops at another width transition. Its source-backed
candidate needs the normal wider-side lead reserve to expose a complete profile.
The local normal Product constructor and independent checks pass in 8.41 seconds:
2,357 pavement quads, stored ratio 7.795, complete 2 mm width proof, no isolated
minima, entrance slope 1.46%, exact boundary heights and no self-overlap.
Rebuilding the same width curve at its old profile plane fails the strict
partial-profile check; the reserve is not a redundant station representation.
The narrowly admitted production rule and its tests are documented in
[ordinary mouth clearance](STREETS_ORDINARY_MOUTH_CLEARANCE.md). A fresh complete
source plan and subsequent normal tasks remain required.

Candidate 21 also exposed a Karl-Marx-Allee source-plan rejection at Section 241.
The actual complete-width failure is real and remains corrected by its existing
source-backed singleton. The subsequent rejection was instead a redundant
reversed-Route distance differing by approximately 1e-15 m, while every actual
native trim, complete profile frame, offset Jacobian and physical gradient
remained byte-exact. The representation-only fix uses no numeric epsilon and
admits no real source movement. Coherent candidate
`frozen_validation_source22_native_mouth_representations` passes all 827 cold
source tasks, producing six groups in 9.844 seconds with unchanged candidate
hashes. `cold22_kma253_source_report.json` and its actual context retain that
evidence. This is a source-plan pass, not a complete mesh/native-gate pass.

Sixty-five focused height, provenance, dependency and progressive-dispatch
regressions pass against coherent candidate 21. They cover the named modules
only; whole-city geometry, height semantics and GUI callback timings remain
separate acceptance gates.

Coherent candidate `frozen_validation_source23_dimensional_reserve` includes
the distinct station-alias and wider-side reserve rules. All 14 affected
regressions pass with unchanged candidate hashes, including full network
construction, reversed/transformed controls, STRICT radii, borrowed-record
rejection, exact native-frame authority and mandatory fresh source proof.
Karl-Marx-Allee again passes all 827 cold source tasks under this policy epoch.
The loaded concurrent run is not a performance benchmark. The corresponding
Bergamo source-plan rerun progresses past the wider-side reserve case, but
terminates at a different singleton: original Section 125 fails complete source
width, its source-backed raw union passes the fresh prefix proof, and canonical
sampling then raises `AnalyticalContactError`. The actual pending context and
unchanged hashes are saved in `cold23_1255_dimensional_reserve_report.json` and
its companion context. No PLAN is ready. Its 715.99-second loaded run peaks at
301.6 MiB; these are diagnostic resource observations, not a production benchmark.
The curve pair and contact ancestry require diagnosis before another recovery
is introduced. Neither the local width-transition constructor nor its regression
fixture proves this remaining source-plan gate.

After memory headroom recovered, candidate 23 also passes the fresh full native
Karl-Marx-Allee import with its normal temporary performance preference enabled:
two components complete, zero errors, 57.29 seconds for Import and 59.74 seconds
including the wrapper/scene save. The process stays below 693.1 MiB sampled
private memory and the candidate hashes remain unchanged. Reports and the
inspection scene are under `kma_native23_full_regression`. This verifies the
earlier source/width corrections against the normal native ground gate. The
source-curb/body policy is still opt-in, so its complete height semantics are
not qualified by this pass.

The subsequent [canonical curb parameter handoff](STREETS_SOURCE_CURB_PARAMETERS.md)
passes five focused tests, including Exact/smooth source maps, reversed native
and multi-Section Route orientation, unchanged physical vertices/labels, invalid
ancestry rejection and real coordinator epoch invalidation. It changes metadata
only and was added after candidate 23; the native result above is not relabeled
as a native qualification of that later change.

The separate Beijing 185 artifact now adds explicitly derived tangent fillets
between retained original corner arcs. Exact finite-parent parameter handoffs
remove the shallow polygon-contact phase without relaxing the original normal
approximation budget. All seven inserted arcs and the earlier original-bank
replacement are fully exposed; no genuine analytical reflex switch remains.
The single median keeps positive area and its actual inner-peak authority, while
the fresh full-width proof passes 3,280 cells. Original mouths, Profiles, source
identities, claims and original parent radii remain unchanged. Participation by
a STRICT parent is explicitly declined by this STYLE-only prototype. Its
surface/height construction and normal runtime integration remain unqualified;
the planar artifact is not a completed Task or native import.

## Evidence location

Generated diagnostics, scenes and logs belong under
`C:/dev/workspace/artifacts/streets/beijing_bergamo_fixes_20261114/`.
The initial audit and partial-import outputs remain in their earlier artifact
folders. Current work is in progress; this file is not a completion claim.

## Ground-only scope decision

The user confirmed that elevated roads remain excluded in this task. At a real
shared source node, ground approaches form an explicit ground-only subgraph.
One remaining ground approach is an open scope boundary, not a dead end. Two or
more ground approaches keep their real connectivity and generate the supported
ground junction. Every omitted elevated incidence is recorded, including its
original Section/end and OSM way IDs. Invalid or missing source records do not
qualify as intentional exclusions. No elevated road is placed at ground level.

Reports distinguish overall `coverage` (still partial when elevated roads are
excluded) from `ground_coverage`. Ground completion requires every retained ground
component to pass construction/publication and no blocked or invalid input.
It does not mean that elevated meshes, bridges or tunnels have been implemented.

The `preparation_03` run has zero blocked ground Sections in both originals. It
defers 11 elevated Sections in Beijing and 90 in Bergamo, and records one and six
ground-only mixed junctions, respectively. These counts are preparation evidence;
full mesh qualification remains in progress.

## Implemented corrections under qualification

- Original OSM node identity controls connectivity. Crossing lines and distinct
  nodes sharing XY no longer create unintended graph connections. Source points
  are not displaced to distinguish their identities.
- Returning Sections retain two endpoint incidences but have one path/profile
  resource. A bounded secondary contact search rejects a coincident zero-sweep
  root and can recover distinct exact circular contacts on a returning road.
- Failed smooth paths can retry with certified evaluation controls. Source
  vertices remain unchanged; the complete curve must still satisfy the original
  source-distance envelope and offset-fold checks. See
  [dense-control smoothing](STREETS_DENSE_CONTROL_SMOOTHING.md).
- Bundle continuation at a branching terminal uses unique original-way
  continuity instead of container order. Explicit roundabout/circular arcs are
  excluded from parallel-carriageway grouping, while remaining ordinary roads.
- Heuristic overlap with another street group does not by itself make a short
  source segment a transverse Bundle splitter. A unique directed continuation
  at the original shared endpoint, with at most a 20-degree local turn, retains
  that segment in its lateral chain. Ambiguous forks, transverse links and
  merely nearby endpoints do not pass this test. Original points and linked
  connectors remain unchanged. This prevents one-member pseudo-Bundles created
  by cutting an otherwise continuous carriageway into separate candidates.
- Close source junctions can defer individual mouths until their shared physical
  exterior is constructed. Consumed source spans require explicit ownership and
  containment; no rejected surface is published as an accepted replacement.
- Width-continuation curves use a bounded refined Bernstein curvature
  certificate when the coarse global certificate is inconclusive. This can
  accept the original curve without moving its controls or reducing its radius.
  If the original curve genuinely needs more reach, a separate finite-budget
  construction extends the incoming leads along their original paths. Profiles
  stay complete at the new mouths, circular contacts stay fixed, and real
  source/scope endpoints may not be extrapolated. Two dedicated tests check the
  resulting radius, profiles, stored-coordinate quad quality and finite-end
  rejection.
- Preferred style radii retain their ordinary fitting range first. Long acute
  common-node forks can subsequently fit the rounded gore nose down to the
  ordinary 0.5 m heuristic. A narrow source-proven median has a separate
  failed-only scale: its nose radius is capped by the genuine hole's area divided
  by perimeter, with a finite smaller-radius family. This preserves the median
  rather than forcing an overlarge cap across it. A gore nose is distinct from
  a vehicle turning radius. Explicit user radii remain strict; fitted values and
  the fitting branch are reported.

These are general geometry/topology policies, not conditions on fixture names or
OSM identifiers. A successful preparation pass is not a successful mesh import.

New shared-boundary policies and direct quad partitions are described in
[compound source ownership](STREETS_COMPOUND_SOURCE_OWNERSHIP.md) and
[direct quad partitions](STREETS_DIRECT_QUAD_PARTITIONS.md). These documents
describe construction contracts; the full-city acceptance result remains the
separate gate below.

Microscopic remnants of consumed internal profiles are covered by
[compound boundary sampling](STREETS_COMPOUND_BOUNDARY_SAMPLING.md). The correction
retains original analytical curbs and complete external profiles, and requires
continuous complete-width containment. Real clipping and larger fragments remain
construction failures.

Additional finite-boundary policies are described in
[merged ground scope ports](STREETS_FINITE_SCOPE_PORTS.md) and
[junction-owned free pads](STREETS_FINITE_FREE_PADS.md). Both retain original
source identities while making explicit which physical boundary profiles are
reconciled inside their shared surface.

## Qualification checkpoint

Earlier full-traversal checkpoints are `full_mesh16_beijing_current` and
`full_mesh_13_bergamo`. Both pass original-source boundary support and have no
unmapped-Bundle preparation errors. Their strict ground gates fail: Beijing
publishes nine of ten components, and Bergamo publishes 32 of 35. The later
interrupted Bergamo attempt is recorded separately below.
These component counts are **not** coverage percentages: the failed large
components contain many streets. Corrected ownership can also change component
counts, so counts from earlier runs are not directly comparable.

The empty-member Bundle guard remains mandatory. Source-chain normalization now
maps the genuine corridors rather than disguising them as independent full-crown
roads. Ground portions disconnected by excluded bridges retain explicit scope
handoffs; multiple ground extents keep their original Bundle lineage. See
[Bundle source chains](STREETS_BUNDLE_SOURCE_CHAINS.md).

The following table records earlier focused checkpoints. Later complete-width
source checks supersede their ownership-plan acceptance; the current remaining
gates are recorded in the continuation sections below.

| Construction | Current evidence | Remaining integration |
| --- | --- | --- |
| Beijing expanded source cycle, including partial and whole Bundle intervals | Current float32-source replay passes with 25,017 total quads; original rims, height/drainage and joins pass | Complete ownership plan, remaining Bundle outputs and fresh native import |
| Bergamo source cycle at 903, retaining both real grass holes | Complete surface passes with 15,049 pavement and 510 grass quads; independently stored coordinates, UVs, seams and overlap pass | Fresh native import and complete progressive parity |
| Bergamo 1,255-Section component | Fresh normal ownership plan passes in 374.87 s. The finite C2 root case passes with 8,715 pavement and 1,167 grass quads. Subsequent normal-backend continuations pass the formerly failing radial-budget fork, long lead/short arc, ordinary drainage and three-junction center-chart cases | Remaining component surfaces, Bundle outputs, native import and progressive parity |
| Bergamo finite scope composition at 1628/1656 | Complete surface passes with 2,936 quads, permitted 1 m corner, stored ratio 7.79, zero minima and joins below 3% | Fresh complete component and native import |
| Bergamo 133-Section component | All 48 normal junction surfaces pass in the complete backend; four formerly failing surfaces now have independent exact-profile and stored-coordinate checks | Complete backend reaches an overlapping Bundle separator; shared-profile body integration and native import remain |
| Beijing finite terminal source cycle | Complete normal surface passes with 3,566 quads, stored ratio 7.96, zero minima and joins below 3% | Complete connecting Bundle body and native import |
| Multi-Section overlapping Bundle runs | The Beijing central group and all four side-junction products now pass in one canonical original-Path graph, including the genuine side-junction grass hole. Finite body and rounded external-cap consumers pass numerical fixtures | Complete real Bundle assembly, runtime discovery and fresh native import remain unqualified |

At that earlier checkpoint, the Beijing 1,072-Section component's ownership plan
passed in 344.85 s. These are preparation-only results from the unchanged native input
records. The completed preparations can be reused for diagnostic normal-backend
continuations; final initial-import qualification still requires a fresh native
run. An earlier `*_failure.json` can predate a later accepted plan, so check the
corresponding run's stdout and accepted-plan artifact rather than treating an
old failure file as the current result.

The Bergamo 688-Section component's seven-site source pool with an acute
unequal-width corner now passes full-width source preparation and independent
surface checks. It keeps both 279.353- and 16.008-square-metre grass holes, with
19,025 pavement and 6,373 grass quads. Its bounded adjacent-curve envelope keeps
the existing two circular fillets, original source routes and zero-height outer
curbs. Stored-coordinate ratios, slopes, drainage, profile joins and exact whole
grass/pavement edges pass. The subsequent fresh complete ownership-plan run
advances to later groups; the focused seven-site pass is not a complete
component pass.

The subsequent 1625/1651/1653/1747 pool now passes with 8,569 pavement and 3,096
grass quads, preserving its 32.729-square-metre hole and the real half-profile
plateau. The retained 1641/1648/1649 group also passes, with 3,224 pavement quads.
Both have independent complete-profile, stored-quality, drainage, slope and
overlap audits. The latter uses source-union preparation only after a typed
local-corner contact failure; full source coverage proves the failed internal
corner is hidden rather than clipping or relabeling it. Complete 688-Section
preparation and subsequent body outputs remain separate gates.

The later 468/1683/1684/1751 pool also passes its complete normal surface:
6,905 pavement quads, stored ratio 7.802, minimum corner sine 0.0518, no
interior minima and entrance longitudinal slope below 1.70 percent. Its bounded
finite-profile search preserves all original paths and complete external
profiles. An adjacent original C2 curve is adjusted only after the established
circular/contact families cannot make the actual zero-height outer cap
internal. The unchanged minimum radius and original plateau remain constraints.
This pass does not settle a separate millimetre-scale road remnant between
distinct physical owners, or qualify the whole component.

Further normal-backend continuations pass the elongated 236 fork (13,046
quads), ordinary four-arm junction 9 (3,743 quads), and the 261 compound with
both real grass holes (12,054 pavement and 4,602 grass quads). The 856 compound
keeps its 402.619-square-metre separator with 11,011 pavement and 2,363 grass
quads. These checks include continuous source-width coverage, stored-coordinate
grading, exact profiles, drainage, slopes and forbidden overlap. Failed-only
computational cuts or grass bands preserve successful earlier chart families.

The multi-Section Beijing Bundle has a fresh ownership plan built from its
native input with all four additional side junctions installed before planning.
It has one owner for each original junction. All five resulting junction
products, four internal body pieces and both external body/cap pieces pass
focused construction. The external pieces have 3,304 pavement/768 grass and
4,634 pavement/1,662 grass quads, respectively. All four pavement halves beside
the two external separator openings pass independent fixed-boundary shape,
height, slope and drainage checks. The complete Bundle plan qualifies with
14,970 pavement and 2,430 grass quads. Its global body/junction overlap audit,
all original source intervals, query ceiling, finite material UVs and geometry
release checks pass. Normal-provider integration and regression checks are now
enabled; a fresh native import remains required.

The normal native-input continuation now rebuilds those side junctions and
generates the 14,970-quad Bundle pavement through the actual provider. It then
reaches a later, separate Bundle whose old separator-cap assumption rejects a
boundary interrupted by another physical mouth. Thus the new Bundle's normal
execution is confirmed, while complete Beijing publication still fails.

Further focused surfaces include the 504-quad asymmetric service junction,
where explicit fixed-chart entrance constraints keep longitudinal joins below
three percent, and the 3,390-quad two-junction compound recovered with the
existing source-certified reflex atlas. Both retain exact incoming profiles,
stored quality, original source coverage and the elevation envelope. The later
Bergamo thin median qualifies with 2,942 pavement and 867 grass quads, retaining
4.899 square metres of grass (91.39 percent of its original source hole).
The subsequent three-junction contact case qualifies 2,865 quads using the
existing one-micrometre source-rim normalization only after a typed analytical
contact failure; contact tolerances are unchanged and all three full original
strip widths are recertified before and after sampling.

The new two-junction body consumer has passing Exact and Smooth numerical
fixtures with actual accepted junction products, real and derived registry
entries, complete physical fronts, and continuous source containment. Its
external-end counterpart also passes multi-Section fixtures with actual rounded
remote caps, including reversed and rotated source traversal. Pavement and its
profile joins retain the three-percent longitudinal bound; the grass cap retains
its existing ten-percent physical slope and drainage limits. These consumers
are used by the jointly prepared multi-Section source path. Complete real
Bundle assembly and joint sampling have focused qualification, while whole-city
coverage and publication remain separate unfinished gates.

The long-lead fork's piecewise analytical sampling qualifies 8,455 quads without
changing its prescribed curves or profile rows. The ordinary drainage recovery
qualifies 3,211 quads by trying coupled spill-floor and slope reconciliation only
after the existing refinement family fails. The three-junction radial core case
qualifies 4,609 quads using bounded fixed-rim charts and complete assembled/stored
quality checks. Each has an independent exact-profile, slope, drainage and
overlap audit. They are focused surface evidence, not final import timings or
proof that subsequent Bundle outputs will succeed.

The current focused Beijing source-cycle surface takes about 106 seconds; its
older replay took about 403 seconds. The Bergamo 903 surface took about 181
seconds in its recorded diagnostic run. These establish geometry
correctness, **not** satisfactory interactive responsiveness. A synchronous
surface builder remains an indivisible progressive-import callback; final
progressive testing must report its blocking time honestly.

The source-footprint union was independently optimized without changing its
boundary or source ownership. A captured 46-region, 11,016-segment Beijing union
dropped from 138.06 to 8.63 seconds using bounded midpoint batches and conservative
spatial rejection. Indexed and exhaustive results retain identical coordinates
and source spans. This isolated timing is not a full-import speedup claim.

`inspect_street_construction_failures.py` is an offline diagnostic only. Raw
preparation captures have not yet passed the normal original-source support step
at the import rectangle. Consequently, some contact failures in that inventory
are expected selection-boundary cases that normal import can recover. Use fresh
native `--full-mesh --require-ground-coverage` results for acceptance.

Ground acceptance also requires a successful overall import: no caught error,
operator result exactly `FINISHED`, outcome `success` with no import errors,
and nonempty `COMPLETE` ground reports containing only `READY` surfaces. The
audit records this as `SUCCESSFUL_IMPORT_AND_COMPLETE_GROUND_V2`. An older
`ground_qualification_passed` field alone is insufficient: a controlled native
late-finalization failure left every surface READY while the import was partial,
and an early second-input failure could retain the preceding input's reports.
Recheck the complete audit record and successful process exit for older runs.
Explicit elevated exclusions may still leave overall coverage `PARTIAL` while
the retained ground coverage is `COMPLETE` and the import succeeds.

No change is ready for synchronization merely because its focused test passes.
The final gate also includes stored-coordinate mesh quality, original supported
network coverage, existing fixtures, progressive publication and visual review.

The `kma_regression_01` native run passes strict ground coverage for the original
Karl-Marx-Allee extract (two components, no import errors). Its saved inspection
scene is in that artifact folder. `regression_01` passes 192 numerical tests
covering paths, profiles, source support, ownership, Bundle preparation and
progressive product events. These checks do not substitute for final Beijing
and Bergamo qualification.

Additional regression evidence includes `bundle_regression_01` (88 tests),
`region_union_tie_validation.json` (eight actual-owner, transformed/indexed union
replays), and focused finite-scope/free-pad tests. All evidence is under the task
artifact directory. The original `full_mesh_03` process has no completed import
report and is excluded from qualification evidence.

For an offline inventory that includes the normal import-boundary preparation,
run `check_street_issue_audit.py` with `--export-inputs --source-support` and feed
the resulting `supported_networks.json` to
`inspect_street_construction_failures.py`. Check each component's
`source_support_error` before interpreting that inventory. Neither diagnostic
publishes or accepts a partial mesh.

The audit harness defaults to disabling the temporary import performance request
for reproducible lifecycle tests. Add `--performance-request` to exercise the
normal preference when measuring interactive-import timing. Record that flag and
other running work before comparing elapsed times: the `kma_regression_01` and
`kma_regression_02` inputs/settings are identical, but the second run's unchanged
path tasks also took several times longer. Those two concurrent-development
runs are coverage evidence, not a controlled performance comparison.

The later `regression_02` passes 198 numerical checks. Focused finite-boundary
checks also cover the open ground-scope spur and a narrow free stub inside a
wider bent through road. Each runs full containment, complete external profile,
drainage, longitudinal-join and stored-coordinate quality checks. They remain
local evidence until the original cities pass the full-import gate.

`kma_regression_03` is another successful strict native import after the
longitudinal-link grouping fix. A read-only comparison of its saved scene with
`kma_regression_02` finds all 328 generated mesh outputs exactly unchanged,
including positions, topology, UVs, normals, shading and material assignments.
Source records/settings also match. `progressive_lifecycle_01` passes six native
test methods covering real-adapter parity, cancellation before construction,
hidden staging cleanup, preservation of accepted results and user edits, and
the extended authored fixture matrix's exact blocking/progressive parity.

## November 15 continuation: source runs and standalone grass

These are additional **focused qualifications**, not final full-city passes:

- The finite merged-terminal consumer is routed through the normal Bundle task.
  It qualifies 20,529 pavement and 3,169 grass quads while retaining the complete
  separated original road remainder. The task declares actual crossroad source
  inputs and physical-port readiness rather than reading cached road meshes.
- Source-cycle side preparation now runs inside the normal ownership stage,
  after all source groups are resolved and before Products. A current Beijing
  source-task replay yields 43 groups over 100 original sites, with one complete
  Bundle source resource and the previously missing actual side owner exactly
  once. Fourteen dependency, exact phase and terminal regressions pass.
- The resulting complete source-owned Bundle passes the normal constructor from
  that current source-task witness, with freshly built physical junction
  Products: 15,920 pavement and 318 grass quads, full original interval coverage,
  exact profile/grass seams, stored shape/height checks, required-owner overlaps,
  UVs, height queries and buffer release. This is not yet a native import pass.
- The long Bergamo separator now passes its normal provider with 21,432 grass
  quads. All 6,756 boundary vertices and whole adjacent pavement edges remain
  exact. Stored ratio is 7.85360 with positive corners; height is its original
  flat plateau. Exact topology, material, UV and 59,620 overlap candidate checks
  pass. Fifteen focused band/patch/old-separator regressions pass. A subsequent
  normal component continuation advances past this separator and stops at a
  different partial end-cap case; the component is not yet qualified.
- Thirty-four existing free, angular, dimensional and ground-scope Bundle tests
  pass after the new dependency routing. The syntax check currently covers 258
  changed/new Python files, without errors.

The new exact-rim band fallback and direct two-/three-quad seed patches are
documented in [direct quad partitions](STREETS_DIRECT_QUAD_PARTITIONS.md).
Source-cycle resources, finite terminal evaluation cuts and their cache inputs
are described in [one-sided paired junctions](STREETS_ONE_SIDED_PAIRED_JUNCTIONS.md).
Both constructions retain the established successful families first.

Remaining gates include the actual remaining Bergamo compound/body/end-cap
cases, fresh complete-city native ground coverage, independent final overlaps,
progressive publication/cancellation and scene review. The ten-site original
source exterior now qualifies its complete Product with 20,854 pavement and
1,888 grass quads. Larger source pools have a separate measured-work admission:
at most sixteen original physical sites, 256 finite source operands, 20,000
source segments and containment cells, sixteen profiles and eight real holes.
All final shape, source containment, height and allocation limits still apply.
The ordinary successful eight-site policy runs first. This local qualification
does not qualify the remaining surfaces in its connected component.

The fresh `full_mesh_15_beijing` run advances through all junction/road products
and the repaired complete Bundles, then stops at a later Bundle whose inner
curbs cross twice over an approximately 61 m interval. Its single grass loop
is therefore invalid; the body needs a real local pavement union. The report
still has 9/10 networks READY. This is progress through different geometry,
not improved whole-city coverage. Its approximately 1,507 s elapsed and 271 s
maximum callback were measured during concurrent development; the blocking run
has no progressive previews and does not establish satisfactory responsiveness.

The subsequent native source investigation confirms a separate underlying
dimension defect related to #121: correction propagation can replace a genuine
two-lane width with an earlier three-lane width. The new propagation boundary
fix and 28 native regressions retain original dimensions while preserving valid
correction runs. `preparation_16_source_widths` reports no structural findings
in either original extract. It changes only three Bergamo widths and none in
Beijing. The corrected-width normal component continuation advances through
the new finite grass closure and reaches a further unequal-visible-range
separator case. It remains unqualified; no final global-overlap success has
been claimed for that component.

### Complete-width guard and failed-run captures

An accepted raw junction-footprint union must contain every fully consumed
original carriageway, not only its centerline. A genuine continuous source-strip
escape now triggers preparation of the full original-source exterior. A proof
budget or unsupported curve certificate is not treated as clipping and cannot
silently invoke that recovery. Both consumed roads and incoming original
prefixes must pass the two-millimeter containment certificate afterward.

This guard found two further Bergamo source-width losses and a Beijing Bundle
passage requiring additional source preparation. A surviving pair of consecutive
original member traversals may retain its complete half-profile ports even when
another member interval is consumed. Neither traversal in that surviving pair,
nor any original piece in either incoming route, may itself be consumed. This
permission creates no additional source connection or ownership claim.

For costly fresh native runs, the audit harness accepts
`--capture-failed-products`. On failure it saves the actual pending numerical
network, products and source/task digests under the selected artifact directory.
It never changes evaluation order, manufactures digests, substitutes products,
commits the failed update or suppresses the error. Such captures support focused
diagnostics; a continuation from them is not a new complete-city acceptance run.

The controlled native capture check deliberately fails three synthetic road
builders. All three real pending contexts are saved with sixteen validated
upstream products apiece, and the incomplete import still fails its ground gate.
Only keys in `validated_product_keys`/`task_digests` have completed evaluation in
that attempt: other cached products in a failed context may not have been
revalidated. Never reconstruct those digests or promote the capture to an
accepted update.
Source/task digests do not authorize reuse after a changed geometry policy or
task definition. Check the implementation and source epoch; diagnostics using
an earlier source checkpoint must say so. Final qualification uses fresh native
imports, not reconstructed fingerprints or relabeled cached products.

For long fresh runs, `tests/trace_street_issue_audit.py` accepts the same arguments
as `tests/check_street_issue_audit.py`. Launch it with the same isolated factory
profile and resource guard. It forwards the original progress callback, logs
normal task starts to `native_task_trace.jsonl`, and writes periodic diagnostic
Python stacks to `native_stack_trace.log`. It changes no builders, task digests,
evaluation order, products or acceptance rules. `source_epoch.json` records
Python source hashes before and after the run; a changed epoch must be disclosed
and cannot support a frozen-implementation qualification claim. Its controlled
three-component native check passes all ground gates and records sixty progress
events with no source changes during execution. These traced runs are functional
diagnostics, not controlled speed benchmarks.

### Original median width changes and displaced prefix contacts

A newly checked Beijing three-site junction retains a genuine widened Bundle
member. Its inner cross-section previously left an 821 mm preparation seam in a
grass hole. A failed-only width support now joins the two consecutive original
member traversals at their real shared junction. It requires unconsumed source
routes, centered unequal widths, STYLE radii, nearly opposite tangents, identical
positive facing crown heights and zero facing profile derivatives. It retains
the original minimum curvature radius and spends only finite original leads.
Outer circular corners and source profiles stay unchanged. Zero-length spans,
outer-boundary roots, outer halves, STRICT overrides and unrelated incidences
cannot authorize this construction.

The normal source exterior and Product qualify 9,056 pavement and 1,277 grass
quads. Independent stored coordinates pass the shape, slope, profile and UV
checks; complete two-millimeter source containment and changed-part overlap
checks against 2,634,111 existing quads pass. This is a focused junction result.

An already owned short dimensional prefix can place an approach's actual
external frame away from the original source node. Source-cycle fillets now
close through that exact original prefix and external route rather than drawing
a bridge between displaced origins. If two internal closure seams cross once,
their entirely internal lobe can be removed only after continuous containment
inside an actual finite original source strip is proved, with analytical curb
clearance. Every physical curb sample is retained. Multiple crossings, removal
of a curb-bearing lobe, insufficient source coverage or exhausted proof budgets
remain failures. The complete resulting junction must still pass source,
profile, physical-rim, height, mesh and publication checks.

### Current regression evidence

`core_regressions_current02` passes 166 numerical tests, including ordinary
roads, Bundles, source ownership, cached updates, original prefix contacts and
progressive product dispatch. The latest syntax pass covers 313 changed/new
Python files. The documentation link audit covers twenty-one changed documents
and finds no broken local links. Later additions require these checks again.

`shared_sampling_regressions_final01` passes 31 tests covering generated
contact-phase admission, unchanged profile/hole/area semantics, failed-chart
dispatch, certified large source pools, and existing region-union behavior.

`kma_regression04_source_widths` passes the original Karl-Marx-Allee strict native
ground gate with both components READY and all 328 mesh outputs present. The
saved-scene comparison finds 288 outputs byte-identical to `kma_regression_03`
and forty changed. This is **not** an unchanged-geometry optimization claim:
four source Sections on Lichtenberger Straße now retain the actual three-lane
style width, 7.5 m, instead of inheriting 5 m from an earlier two-lane segment.
The fresh source snapshot verifies the original lane tags and resolved widths;
it has no structural findings. Other source geometry, connectivity, settings
and Bundle semantics are unchanged. The resulting total is 354,292 quads,
versus the previous 350,623. The run's 58.46 s elapsed during concurrent
development is functional evidence rather than a controlled speed comparison.

`progressive_lifecycle02` passes eight native methods, including the extended
fixture matrix's exact blocking/progressive parity, normal OSM adaptation,
cancellation at different construction/publication stages, hidden staging
cleanup, and restoration after injected precommit/postcommit failures. This
qualifies those lifecycle contracts. It does not establish a short blocking
callback or complete-city parity for the still unfinished original imports.

### Fresh native continuation after the complete-width guard

`full_mesh16_beijing_current` runs the original OSM through normal preparation
and construction with no substituted products. Its ownership task succeeds, but
the large component stops at `junction_surface:122` on a direct-quad layer
budget. The nine smaller components are READY; the overall ground gate still
fails. Its 235.49 s elapsed is an earlier failure point, not an import speedup.
The saved scene is explicitly named `beijing_center_unqualified.blend`.

This run saves the actual failed coordinator context, including 4,328 validated
products and their real source/task digests. The numerical diagnostic loader
translates only the Python package alias and verifies the original source
fingerprints; it does not reconstruct task digests or change source identities.

The all-group comparison against the older captured context identifies eight
changed existing physical exteriors, plus the separately prepared side owner.
Consequently, qualification of the three difficult local source cases does not
qualify every other changed exterior. The initial additional gate was the
122/126 group, followed by the remaining changed groups in the normal task
order. These are additional validation gates, not intentionally
excluded roads. Elevated-road exclusions remain the original explicit scope.

The focused 122/126 diagnostic subsequently qualifies 2,470 pavement quads after
normalizing only the generated exterior sampling, with at most one micrometer
of projection error. All 68 original profile samples remain exact, the original
area operands and holes remain authoritative, and independent stored geometry,
drainage, slope, UV, query and local-overlap checks pass. That first diagnostic
used an explicit sampling-admission override. The genuine native16 coordinator
continuation subsequently passes the normal 122 task in 239.89 s; its complete
mesh, physical rim, peak map and original profile map are byte-identical to the
independently audited diagnostic. The saved task digest is the real runtime
digest, not a reconstructed value. This qualifies that repaired task and its
dispatch, not the remaining changed exteriors or complete Beijing import.

The current Bergamo checks also qualify the complete 1,255-Section ownership
plan, with thirty source groups, and the normal 213 pavement product, with
8,663 quads. The wholly owned/free-member Bundle 6 separator passes a separate
constructor diagnostic with 771 grass quads. Its report deliberately records
`task_qualified=false`: private constructor inputs are not a substitute for a
real dependency task or native import. The later complete-OSM attempt is recorded
below; neither city has passed the final complete-ground qualification.

### Source-helper review corrections

A follow-up review fixes three bounded helper cases: all free-end-plane
crossings are counted before requiring a unique oriented separator closure;
an eligible owned/free group without a prepared physical layout returns the
existing typed unavailability; and canonical source sampling enforces its
aggregate allocation before creating arrays. Rendering-only source-contact
operands may sample just their owned intervals while retaining compulsory cuts.
The shared source-cycle grids must sample the complete original finite paths:
future mouth candidates and remaining-strip proofs also consume those grids.

Twenty-six focused tests pass, including full versus partial owned intervals,
reversed/smooth sources, the ten-kilometer unused-tail witness and existing
source-cycle/contact/companion behavior. Captured 213 and 199 source operands
remain byte-identical. The Bundle 6 diagnostic still has exactly the same 771
quads and XYZ data, with its stored-quality and overlap checks passing.
These are local review checks; a native run started before the changed source
epoch is not automatically a qualification of the corrected implementation.

A later captured-source ablation exposed a contract error in the first
owned-interval optimization. A sparse future tail represented 65.99 m by one
chord and falsely requested the far owner. The complete historical grid found
the last actual overlap at 52.496 m and a valid mouth at 52.694 m. Complete
future sampling is restored under the same 20,000-row allocation bound;
rendering-only maps explicitly reject use by future-source searches. The
ownership-plan task has a separate policy version so cached plans from the
incorrect sampling epoch rebuild while raw paths, profiles and local layouts
remain reusable. Input records were identical across the compared captures;
the changed ownership closure was a sampling-contract regression.

The additional changed Beijing groups 486, 629, 784 and 1132 also pass direct
normal surface construction and independent stored/physical grading, slope,
drainage, exact-profile, source-area, UV, query and local-overlap audits. Together
with the earlier focused cases, this covers all nine changed source groups in
the comparison. Genuine saved coordinator outputs for 122, 486, 629 and 784
match their independently audited geometry. It remains necessary to finish the
coordinator, all Bundle bodies, complete-update overlap checks and fresh native
publication. A local pass is not a whole-component acceptance.

### Bergamo attempt 17 and coherent test candidates

`full_mesh17_bergamo_current` finishes after 3,215.45 s with a partial import.
The actual 688-Section component stops at `junction_surface:483` on assembled
quad grading. The following 1,255-Section component is interrupted by an import
of `_sampling_intervals`: the process had loaded the older `compound_cycles`
module before the review edit, then loaded the newer contact helper afterward.
This is a mixed test epoch, not a missing function in the current coherent
source tree. Three small components are READY and only four component records
are published before the interruption. These counts cannot be compared with a
complete thirty-five-component traversal as a coverage percentage.

The V2 acceptance predicate rejects the full record. Its 986.67 s maximum
blocking callback and absence of previews do not qualify progressive user
experience. The saved scene is `airport_bergamo_unqualified.blend`.

The real completed products before 483 are useful independent evidence:
285 passes with 2,102 pavement quads; 318 passes with 20,504 pavement and 1,940
grass quads; 468 passes with 6,905 pavement quads. Each passes fresh profile,
complete-source, stored/physical quality, slope, drainage, UV and local-overlap
checks. These are actual coordinator products, not substitutions from an older
source graph. Their metrics do not qualify the failed component or whole city.

Long future native checks should launch a frozen working-tree candidate made by
`tests/freeze_street_test_source.py` and verified with its source manifest.
Maintain source in this repository; the generated candidate is disposable test
evidence under street artifacts, not a recovery checkpoint or installed addon.
The copied `tests/trace_street_issue_audit.py` supplies the same native test
plus task traces and source-epoch checks, using the chat's isolated factory
profile and existing resource guard. Never edit the candidate while it runs.

### Attempt 18 and the corrected source-sampling contract

The coherent `full_mesh18_bergamo_frozen` traversal finishes all thirty-five
components in 31,590.01 s: thirty-two READY, three failed. Its source manifest
and trace show no source changes during execution. This is still a failed
whole-city import, and its duration does not qualify progressive responsiveness.
The failures are the 688-Section component's surface 318, the 1,255-Section
component's Bundle 9 source inventory, and the 133-Section component's source
plan. The candidate predates the corrected complete-future sampling policy.

Comparing actual upstream records from attempts 17 and 18 found identical raw
Paths, Profiles and local physical intersections. The sparse future grid alone
changed derived source ownership in the 688-Section component. Restore the
complete canonical finite-source grid before diagnosing its enlarged compounds;
those enlarged groups are not evidence for changing original OSM geometry.
The policy version is included in projected layout, surface and seam signatures
as well as the PLAN task, so same-membership stale geometry cannot survive a
normal update merely because a projected signature omits the parent PLAN hash.

Normal coordinator source-plan rebuilds with the corrected policy now pass for
Beijing (256.59 s), Bergamo's 1,255-Section component (345.77 s), and Bergamo's
133-Section component (159.54 s). Only PLAN was rebuilt in these checks; genuine
raw input identities and task digests were preserved. The thirty groups in the
1,255-Section plan match the coherent attempt-18 geometry exactly. These are
source-plan checks, not complete mesh or native-import qualifications.

Current remaining mesh cases include Beijing's overlapping Bundle 18 body and
separator; Bergamo Bundle 9's missing accounting for a whole Section genuinely
owned by an ordinary routed junction prefix; and continued compound-surface
qualification. The Bundle 9 source-width proof passes: the missing original
visit belongs in the ownership ledger, not in an invented source connection.

The corrected 133-Section continuation stops at surface 995 after 751.34 s.
Its original physical XY and provenance spans match the previously audited
source exactly. The immediate failure is inside an optional normalized-rim
candidate: finite hole-cut exhaustion was untyped and escaped before the raw
construction alternatives could run. A distinct `PhysicalHoleCutError` now
allows that candidate to decline under the existing one-attempt budget.
Invalid original source proofs and unrelated implementation errors remain
fatal. Twenty-seven focused sampling, normalization, and physical-hole tests
pass; a real coordinator continuation from the genuine prefix is running from
the coherent `frozen_validation_source19a` candidate. That continuation now
passes the complete 133-Section numerical backend in 1,426.73 s, including
surface 995, both Bundle plans/separators, road outputs, UV fallback and the
final changed-part overlap check over 244,030 quads. It reuses genuine validated
upstreams, substitutes no mesh, and records identical candidate source hashes
before and after. This is coherent pre-coupling-policy component evidence,
not a fresh whole-city native import or a performance benchmark.

Independent surface-995 checks pass on the actual coordinator product: 18,828
pavement quads, stored ratio 7.98123, stored corner sine 0.04286, maximum physical
slope 0.09048 against 0.1, entrance longitudinal slope 0.00669, and no interior
minimum or drainage fill. Complete 2 mm source containment, exact full profiles,
physical/float32 grass seams, UVs and local pavement/grass overlap checks pass.
A separate all/all check covers all 153 published parts, including pairs where
both parts came from the genuine validated prefix: all 244,030 physical and
float32 quads pass, with 682,419 and 682,150 candidate pairs respectively.
Every retained junction's physical elevation stays within its datum/crown
envelope; the rounded buffer stays within the corresponding stored envelope.
A static Blender inspection scene and top/oblique/control-overlay images are
under `native19_component133_visual/` in the task artifact directory. Its object
height exaggeration is explicitly ten; mesh coordinates retain physical height.
The scene visualizes accepted numerical buffers and does not claim native city
publication or an editable source network.

The corrected 688-Section source-plan task also returns normally in 2,287.0 s,
restoring thirty-three groups, the original ten-site group at 318, and separate
ownership at 468 without acquiring 472. All 2,818 genuine raw inputs retain
their identities and digests. Every published float32 rim matches attempt 17
byte-for-byte; the largest double-coordinate difference is 42 nm, with unchanged
memberships, claims, row counts and source-span identities. A strict diagnostic
double-byte comparison rejects roundoff, while the actual source-plan task and
its normal guards pass. Its numerical continuation uses the matching immutable
19a source epoch; it does not relabel that saved PLAN as the later coupling policy.

### Local update dependencies and owned-prefix accounting

The Bundle planner's broad hydration closure must not become a declaration
that every potential distant road changes the local separator. Construction
still waits for those inputs, but hydration-only route layouts project onto
the actual mouths reaching local Sections. Local Path/Profile inputs, actual
physical-owner surface digests, and complete visible port stamps stay in the
signature. The separator's own input Sections exclude the extra hydration
closure. Eight dependency and partial-cache tests pass, including a real
shortened continuation where the remote mouth does enter the local crossroad:
that edit rebuilds its Bundle plan and road, and all incremental meshes equal
fresh meshes. The unchanged distant-mouth case retains the original road,
separator, and furniture product identities.

The Bundle 9 ownership-ledger correction records whole and partial original
routed prefixes without creating body visits, source Ends, or SharedClaims.
Twenty-two focused inventory, validity, runtime and body tests pass, including
reversed/rotated sources and explicit duplicate, missing, wrong-owner,
wrong-Path, order, and proof rejection. A normal 1,255-Section backend
continuation uses the genuine corrected PLAN context and current task registry;
its actual Bundle 9 task now passes in 12.39 s with 6,695 pavement and 1,499
grass quads. All seven original Sections remain ordered in the ledger, including
the whole 271 and partial 356 prefixes owned by junction 541. Independent stored
pavement ratio/sine are 7.57051/0.04804; grass ratio is 7.58797 with positive
corners. Bundle 6's actual normal separator task also passes with 771 quads and
stored ratio 6.56328. UV fallback and global overlap remain pending because the
component subsequently fails at separator 21 after 3,624.77 s.

Separator 21 reaches two different original crossroad Streets, rather than the
single connecting Section required by the established cap branch. The original
incidences and existing physical ownership are under review; this is not yet
evidence that concatenating a same-Street chain would solve that cap. No
coordinate-based connection or replacement cap has been accepted. Post-launch
source changes in the 1,255-Section diagnostic are recorded separately from its
actual compiled task registry and genuine digests; final current-epoch native
qualification remains required.

The combined integration suite passes 63 methods with no source changes during
the run. A subsequent coupling-proof cache refinement has twelve focused
methods passing: its dependency scope includes all original Bundle Path/Profile
records used by the admission witness, unions inherited scopes, and versions
the PLAN policy so an earlier proof-bearing plan cannot silently hide the scope.
Sampling policy and source geometry remain unchanged by that cache refinement.

### Source-only Beijing side/terminal coupling

A new bounded admission family uses original source data before any junction or
road mesh exists. A unique short original continuation only selects a candidate;
it is not a geometric proof. Admission additionally requires a visible original
carriageway overlap reinverted on the analytical curbs, and a real side-junction
pavement overlap beyond the combined chord allowances. The ordinary exterior
builder then supplies all finite-mouth, complete-width, original-incidence,
rounded-boundary and existing resource-budget proofs. Repeated admissions must
reduce the original owner forest and rediscover candidates under the new cuts.

On Beijing Bundle 18, this builds exactly the five original sites
122/125/126/137/821, consuming original Sections 52/53/54/62/69. The full future
grid finds finite cuts on 61 and 307; neither the distant head group nor the
whole 307 is acquired. The actual normal PLAN task passes in 384.95 s with
forty-two groups, reusing all 4,290 genuine upstream records. It rebuilds PLAN
only and preserves the existing prepared Bundle resource.

The corresponding actual surface task still fails after 578.60 s at bounded
physical-hole/visibility partitioning. Source preparation is therefore qualified
for that recorded epoch, while the five-site mesh, height field, body and
separator remain unqualified. Existing paired-interval computational extensions
are being tested on the exact fixed rim before introducing another chart family.

Those extensions alone do not solve the partition. The recorded failing leaf
lies on the original Section 69 bank: an internal cut normal to unpaired
crossroad 52 creates a long thin wedge. The existing partial-member cuts also
ignore off-center pavement intervals because only whole Bundle claims are
marked paired. An artifact-only policy retaining one validated partial-member
axis, all actual pavement intervals on its normals, and the existing hole
bridges partitions both original roles successfully in about 2.3 s, with exact
canonical rims/claims/holes and the original angle, clearance and cut budgets.
The subsequent raw chart trials still fail the cell or quad-quality guards;
partition success is not a qualified surface. Generated-phase normalization and
bounded computational cadence are the next checks; no physical source changes
or work-budget increases have been promoted.

### Bergamo 688-Section continuation: strict height feasibility

The coherent 19a continuation passes surfaces 23, 27, 135, 147 and 285, then
fails at 295. This is a raw-chart height-projection case, not an optional
normalized candidate. Its last bounded iterate correctly fails: physical slope
0.176000685 exceeds the inherited 0.175999996 limit, the stored slope is worse,
and a 1.206 mm drainage spill remains. Small iterate movement is not permission
to accept that geometry.

An artifact-only failed-caller retry uses the already existing strict feasible
mode: target 0.16950 above the immutable fixed-gradient floor of 0.15000,
acceptance at the unchanged original limit, and the same 1,000-iteration/96 MiB
limits. It finds a feasible intermediate mesh in 107 iterations. The complete
normal Product diagnostic and independent audit pass with 6,666 pavement quads,
stored ratio/sine 7.85393/0.05305, stored slope 0.17548, joins 0.00482, no minima
or drainage fill, exact boundary/XY/quads, complete source-width checks and UVs.
The failed-only integration is now present, with 24 height/projection tests
passing. A coherent immutable source variant containing just the approved height
change passes the genuine coordinator task in 31.15 seconds; its independent
product audit also passes. The subsequent ten-site surface 318 fails a separate
planar mesh-grading check. No original limit or validation tolerance is relaxed.

### Failed-only computational cut endpoint reuse

The corrected ten-site surface 318 exposes a cut-phase issue. Its remaining
short hole-rim chords are generated intersections inside existing analytical
boundary chords, rather than compulsory original source or profile points.
Changing the height solver cannot repair those planar cells. Blanket reuse of
all nearest boundary vertices also fails the existing grading checks.

An explicit `PhysicalHoleCuts.near_endpoint_fraction` option instead reuses the
nearest existing chord endpoint only when a generated contact lies strictly
within one eighth of that endpoint. One eighth is the reciprocal of the
unchanged 8:1 edge guard; this avoids creating a very short cut fragment without
moving the physical boundary. Exact endpoint contacts and every original
physical/profile vertex remain fixed. The option accepts finite fractions from
zero through one eighth; zero retains the original cut references and output.

The source-cycle atlas tries all three existing original-corridor candidates
first. Only their failure admits the same three corridors, in the same order,
with endpoint reuse. The complete family therefore has at most six candidates;
each keeps the existing cut, source, memory and 50,000-cell limits. Full source
width, holes, physical/stored quad quality and all height/publication checks
remain mandatory. This changes computational cuts, never road incidence,
widths, analytical curbs or corner radii.

The captured 318 normal Product diagnostic passes with 19,677 pavement and
1,588 grass quads. Independent stored pavement ratio/sine are
7.76941/0.04868; stored slope is 0.16540 under its inherited 0.17600 limit,
entrance joins are 0.01150, and no minima or drainage fill remain. All 3,176
compulsory canonical rim points are byte-exact; extra points lie on those same
chords within 1.6e-13 m. Original area operands remain the same objects, all
three holes remain, complete 2 mm source/prefix checks and 72,032 overlap pairs
pass. The matching immutable-source coordinator task also passes in 914.87 s,
selecting the fifth candidate after all three default corridors and the first
endpoint-reuse corridor fail. It saves the genuine 2,885-key validated prefix;
an independent audit of that actual task Product passes with byte-identical
geometry to the constructor diagnostic. This remains numerical component
evidence, not a complete current-epoch native import.

A subsequent source-labelled height audit exposes a separate existing defect:
the port-free-hole rule assigns the crown ceiling to full-crown outer curbs
whose original heights are zero. The accepted 318 task has 436 such checked
edges at +0.08 m; the accepted 483 task has 136. Their earlier mesh, source-width,
profile, slope and drainage checks do not establish that missing source-curb
height contract. The cut repair is qualified for planar geometry and complete
source coverage; complete height qualification remains pending the shared
source-labelled-hole correction. True Bundle inner plateaus must keep their
original heights.
