# Experimental Geometry Nodes streets

**Development checkpoint, 19 October 2026:** the current GN sources, exact
private catalogs, inspection scenes and outstanding limitations are documented
in the [GN handover and restart guide](../../prochitecture_gn/docs/streets/blosm_integration/GN_HANDOVER_2026_10_19.md)
in the sibling development repository. The latest private catalog is not the
bundled library described below. A future **Python / Geometry Nodes** generation
choice will share the prepared inputs and geometry contracts; the new Python
solver is not implemented yet and is distinct from legacy generation.

The `buildings_gn` integration targets Blender **5.2**. It imports actual prepared
Blosm Streets into clean editable meshes and uses one shared GN implementation
for carriageways, supported intersections and two-member Bundle separators.
**Qualification is incomplete.** Whole maps can contain rejected or unsupported
parts; inspect the coverage report and native red diagnostics.

## Import

1. Install/update this Blosm checkout and restart Blender after replacing Python
   modules. The distribution must include `assets/prochitecture_streets_v8.blend`
   (revision 8.14), `assets/prochitecture_streets_original_spans.blend`
   (revision 1.14), and the `way/gn` Python package.
2. Enable experimental features in Blosm preferences. Use the existing realistic
   OSM import mode, enable highways, and select **Geometry Nodes (experimental)**
   in the street backend control. Existing product/UI availability rules apply.
3. Choose any OSM source supported by Blosm. For an initial flat-ground test,
   disable buildings, forests, aviation and other unrelated components. A
   streets-only native import needs no external building/vegetation assets folder.
4. Start at mesh resolution **1** and height exaggeration **10**. Exaggeration
   **1** displays physical profile rise. The current flat-ground profile rise
   defaults to 0.08 m. Terrain placement is a later implementation stage.
5. Inspect the `Blosm / Geometry Nodes streets` collection and its import report
   in Blender's Text Editor. Reimport creates a separate owned collection.

An optional street-library path can select a compatible external `.blend`.
It uses the same [general asset loader](ASSET_LIBRARIES.md) as bundled Base/Pro
libraries. This does not change the required street interface or asset revision.

## Edit and inspect

### Connected-junction preview

The additional **Street construction → Connected junctions (preview)** choice
loads `assets/prochitecture_streets_original_spans.blend` (revision 1.14) through
the same universal asset service. The existing Current choice and 8.14 library
remain available. Both use the normal OSM import operation; no fixture filename
or junction ID is embedded in the construction policy.

The newer provider keeps complete editable Streets, addresses their original
Section ends and loads one shared GN root for its components. Its Streets
sidebar exposes cell size, resolution, exaggeration, pavement and grass materials, and
width/offset edits for a selected Section. Saved vertex and width edits update
natively, including with script execution disabled.

Revision **1.14** adds road runs bounded by a change to an excluded category,
including runs with no intersection at either end. The complete source Street
remains editable; its excluded footway Sections retain their original IDs and
are not assigned carriageway surfaces. No junction or free endpoint is invented.
Dimensions must be constant across the selected run; Bundle membership and
dimensional transitions remain explicitly pending in this branch.

The source-built release passes all 24 asset checks, seven universal-loader
tests and isolated normal imports of two different OSM inputs. Exact saved-mesh
comparisons preserve the qualified candidate's coordinates, quads, heights and
ownership. The adopted 839 production files match both tested packages
byte-for-byte. Berlin now has 179,861 quads, including the previously omitted
Schulweg road Section; Nashville retains 155,433. Nine synthetic fixtures and
the new real road pass native vertex/width edits, profile bounds and quad
topology checks at resolutions 1 and 2 in the qualified candidate.

Revision 1.13 added separators between two complete Bundle members attached to
the same intersection at one end and genuinely free at the other. It follows
their existing inner pavement boundaries and retains the previous pavement
exactly. The release passes all 20 asset checks and isolated imports of two
different OSM inputs, with exact geometry against the qualified candidate.
Berlin has 178,133 quads and 201 owners represented by faces, plus one checked
matching-profile seam with no area; Nashville has 155,433 quads and 178 owners.
The free-separator branch passes native edits and recovery at resolutions 1
and 2. Partial Bundle extents, interior islands and general overlapping-member
union remain unfinished. The existing combined gallery still awaits migration.

Revision 1.12 added original two-ended connections, bounded transitions between
unequal profiles, and supported road/category handoffs. It retains the original
rounded corners, profiles and native source-footprint guards. Stable evaluation
of the quad energy difference repairs the measured coordinate-frame sensitivity;
it intentionally corrects some interior layouts while preserving fixed boundaries.
Both candidate installations reproduce the qualified geometry exactly and pass
independent R2 profile checks. Native Edit Mode checks pass at R1/R2. The final
1.12 library passed 18 contract checks. Both final bundled installations passed
normal import without candidate overrides, and independent reopening preserves
all 200 Berlin and 178 Nashville candidate surface owners exactly at R1. These
results do not establish complete Bundle coverage or interactive performance;
the existing combined gallery still awaits migration.

The original-span feature now opts into ground-level building-passage admission:
`tunnel=building_passage`, no bridge, zero or absent layer/level, no indoor
designation, and surface or absent location. Original tags are preserved. The
legacy carriageway feature keeps its earlier scope. Revision 1.12 supplies the
native seam and cross-Street category-boundary construction needed for these
newly connected approaches. Passage admission alone does not establish that
every newly selected connection has a generated surface; inspect coverage.

Revision 1.8 adds native local failure containment, complete unattached Streets,
and road junctions with excluded pedestrian/path attachments. An excluded
attachment keeps its original identity and is reported separately; its surface
is not generated as a carriageway. At least three in-scope road arms are still
required. Revision 1.9 also retains complete Streets containing road and path
Sections, while generating only their selected carriageway runs. Interpolation
separates at original category boundaries, so an excluded footway bend does not
invalidate its neighboring road. Ordinary all-road Streets retain one curve.
Category changes at native junctions, unattached category-bounded runs, wider
Bundle coverage and some logical connections remain unsupported.

At a duplicated Section seam, move both recorded source endpoint copies
together. With the street tools enabled, select either endpoint in Edit Mode,
then choose **Select Section Joint Copies** in the Streets sidebar or F3 search.
Press **G** to move the selected pair normally. The command uses the explicit
Section binding; unrelated coincident vertices stay independent. If the other
Section belongs to a separate visible source object, that object joins
multi-object Edit Mode. A hidden or uneditable counterpart is reported first.
Selection is undoable and changes no coordinates, attributes or topology.
There is no background selection handler. Moving only one copy without this
selection step still fails the solver's continuity check; automatic propagation
is not implemented. Exclusive footway control edits leave road geometry
unchanged; supported road controls update the pavement natively.

At revision 1.9 publication, the bundled library passed ten interface checks and
two isolated street-capable installations passed normal import with Python reads from the development
checkouts forbidden. Independent reopening with street scripts disabled confirms
that all five generated component meshes preserve the qualified candidate's exact
coordinates, faces, physical heights and surface ownership. These are correctness
checks, not responsiveness measurements or full-map qualification. Those
installation records precede the subsequent source-preparation repair below;
they do not establish package equivalence for later Python changes.

The latest preparation repair preserves genuine free member ends through
logical Bundle aggregation. It first fixes Street concatenation overwriting an
outer connector, then retains immutable original endpoint identities before
aggregation. A missing or unmapped junction is never treated as evidence of a
free end. This uses the existing revision 1.9 GN approach branch. The normal
Berlin import now contains 178,292 R1 quads and four plus one unresolved routes;
all 190 preceding surface owners and all three Nashville meshes remain exact.
The restored free cap passes profile, actual endpoint editing, local failure
and recovery checks at R1/R2. Broader Bundle and separator coverage remains open.

Revision 1.5 introduced automatic construction of eligible closed two-member Bundle
separators, preserving their surrounding pavement. The separator declaration,
count and grass-material inputs are required by the renderer's asset preflight.
Missing inputs or a mismatched revision reject the library before source
objects are authored. The default original-span asset build includes this branch.
Two numerical repairs retain valid endpoint corners and recover four previously
rejected high-resolution height projections. Selective separator-boundary welding
preserves the geometry and reduces the measured large-network edit cost; it does
not yet meet the interactive performance target.

The independent provider arrangement introduced in revision 1.4 remains in use.
All patch providers share one GN tree; all corridor/free-approach providers share
another; all network outputs share the assembly tree. A strip references its
existing Street and its adjacent live patch providers. Editable Streets are not
duplicated. Common numeric controls propagate through native property drivers,
and the final output assigns pavement material. The universal asset loader
requests those shared trees from the same library, plus the common diagnostic
display and material when local failure support is enabled. Grouped input objects
remain supported but can invalidate more providers when edited. Global control
changes still require rebuilding all affected providers, and the final assembly
still checks the full component.

This preview currently requires **Smooth interpolation**. Managed topology
commands are not yet connected to its new binding compiler. Local geometric
failures preserve independently accepted surfaces and show native diagnostics;
missing or inconsistent ownership still rejects the affected component. Bundle separator coverage and interactive
performance remain incomplete. Do not interpret successful parsing/import as
proof that the whole map was generated.

A bounded preparation rule assigns eligible short original spans to the existing
two-junction construction. The candidate threshold is original control-path
length less than the span's largest width plus its two prescribed terminal
radii. Already owned logical sites are excluded; competing candidates are
ordered deterministically. This selects an owner, not a solved surface: all GN
geometry checks remain mandatory. It never changes source controls or radii,
and ownership remains stable during ordinary coordinate edits.

The selected three-site real import passes independent saved-file editing and
geometry checks, with 3,512 quads and exactly the same XYZ/faces as its earlier
reference. The revision 1.9 isolated Berlin import produces a 70-site component
and one smaller component, totaling 177,877 quads at resolution 1. Their
unresolved-route counts are 5 and 1 respectively; whole-network coverage
remains partial. The same construction was independently qualified for all 216
entrance profiles at resolutions 1 and 2, including crown limits, oriented
boundaries and positive quads. The R2 candidate contains 442,656 quads.
The isolated Nashville import generates three components: 149,695 and 5,404
quads for its connected networks, plus 396 for an unattached Street. Its 95
selected supported Streets are generated and 18 Streets are deliberately
excluded. Neither sample filename nor any of these IDs selects a production
construction branch.

The additional compound branch uses a direct rounded exterior only when the
original retained-corner construction fails and the connecting Street has two
controls, compatible full profiles and prescribed radii. It checks the entire
absorbed carriageway footprint. It does not shrink radii or alter source
controls. Revision 1.2 batches directional clearance queries while retaining
the same samples and checks. Both Nashville meshes, all three Berlin meshes
and six dense regression networks remain exactly unchanged. A fresh normal
Nashville import also matches the independently checked reference after
reopening without street scripts. Broader Bundle coverage, quad quality and
responsive performance are unfinished. This preview is not full-map-qualified.

### Current construction

Each complete Street normally has a source mesh object. Its Sections are directed
edge chains with `st_section_id`, `st_width`, `st_offset` and `st_interpolation`
attributes. Sources have no isolated metadata points. Move source vertices in
Edit Mode: the saved GN graph updates without rerunning a Python script. Enable
Blosm for its Streets sidebar, parameter controls and managed topology commands.
Do not register an old `OPTIONAL - Street editing tools.py` at the same time.

The default interpolation is the accepted optimized cubic; the other menu option
is Exact polyline. Invalid bends and mismatched mouths produce local diagnostics.
Undo the edit or adjust the inputs to recover. Arbitrary topology edits require
the managed street commands/rebinding. Source object scaling is not the normal
editing workflow; change vertices and physical parameters instead.

The sidebar's **Last explicit check** is the result of the last import or managed
refresh. Native red diagnostics update live. Ordinary GN rejections keep live
editing enabled; only broken topology bindings need a one-shot Python recovery.
The registered observer does not scan/rebuild v8 geometry after each vertex move.

Each connected output component references the **same GN datablock**, loaded from
the asset library. Separate controller objects permit independent source
collections and local diagnostics. Source objects and their controller share a
nearby local frame. The saved projection origin and source-coordinate offset
recover the original import coordinates; reversible gallery placement is separate.
Widths and profiles are metres, including when the scene uses a different unit
scale. Source coordinates are converted to Blender units at import.

## Coverage and limitations

The report contains prepared IDs, OSM way IDs, original node-pair segment
intervals, width-origin metadata, deferred items, native diagnostics, and per-Street
coverage. Segment intervals describe contributing preparation segments; they are
provenance, not new editable GN indices. Later preparation splits may retain the
ancestor interval list.

Coverage measures **prepared input centerline length**, not surfaced area or the
length remaining after mouth trimming. `GENERATED` requires every Section and its
active junctions and explicit transitions to be present without relevant native errors. `PARTIAL` includes
local rejection and explicitly clipped boundaries. Deferred categories and test
selection exclusions are recorded separately. A successful import operation does
not mean that every street was generated.

Unsupported junctions no longer suppress all neighboring roads. Their adjacent
Sections stop before an explicit exclusion area, using the incoming widths and
radius policy. Authored points remain intact. These are visible gaps awaiting a
proper junction solution, not dead ends inferred from OSM. They remain partial
coverage. Short/bent approaches can themselves fail the native trim checks.

Current restrictions include undeclared/general lane transitions, partial Bundle
membership and internal Bundle junctions, general curved nonparallel divided
mouths and outer arm gaps above 210 degrees. Convex
rounded corners on bent through roads are supported up to that limit, subject
to the unchanged mesh validity checks. Grade-separated features
retain metadata but are not terrain/bridge surfaces. An audit records legacy
XY connections without shared OSM endpoint nodes. Affected carriageway Streets
remain unsupported source context, rather than being published with unverified
connections. Events involving deferred components or unmapped provenance are
reported separately. This policy does **not** replace the graph with a full
grade-aware connectivity algorithm.

The parallel-Street detector preserves closed paths and repeated control visits
instead of deleting nonconsecutive duplicate coordinates. Such roads remain in
the input inventory but cannot become open Bundle members; closed carriageway
sweeps still require their own provider. Bundle reports include endpoint
locations and all owned Streets, distinguishing the ordered parallel members
from internal attached paths.

Lane-transition preparation plans chains against the original graph before
changing Section ownership. It preserves source controls and dimensions while
assigning each Section once and remapping each junction incidence to its own
completed Street. The saved `ownership_audit` lists any duplicate Section IDs;
ambiguous source owners are retained with diagnostics rather than silently
deduplicated or generated as overlapping roads.

Minor-junction cycles are detected before linear Street concatenation. Their
original segments and complete connector rings remain explicit junctions, so
preparation cannot create an infinite item list. The report records these
reclassifications under `preparation_audit.minor_cycles`. A complete Street
returning to a major junction remains inspectable with a closed-path limitation;
it is not silently omitted. This does not add a general roundabout provider.

Complete Bundle members may contain several compatible Sections. Their profile
and interpolation mode must agree along each member. Width/offset differences
between consecutive Sections require a declared finite transition; the two
members can also have different widths. GN joins their existing cubic spans
and retains each Section ID on its surface faces. Both endpoint copies at a
Section seam must move together when editing that seam. Incompatible edits are
diagnosed and recover natively when restored. Explicitly clipped unattached
Bundles can use a common cut plane; ordinary free endpoints remain unchanged.

Ordinary carriageway `SymLane` and `SideLane` seams now have an explicit width
transition provider. Width and offset stay constant on each authored Section;
the generated surface uses finite quintic tapers around the seam, retaining
the full crown and zero edge height. Live width changes update both halves and
their extents. Too little clear length produces diagnostic 65. Both copies of
the shared Section endpoint must move together. Revision 8.11 applies these
declared transitions within complete two-member Bundles and derives the live
separator and pavement overlap from their varying physical boundaries. Free
terminal cuts retain their centerline-based axis when the envelope width changes.
Prepared transition items with their own centerline still need an extent provider.
General curved attached approaches remain limited by the straight-lead checks.

When the through carriageway changes width or offset at a deferred pedestrian
or cycleway attachment, preparation adds the same finite taper declaration to
the original `Intersection` item. Its `surface_treatment` remains
`DEFERRED_ATTACHMENT`: attachment IDs, location, side and through Section-end
references are preserved for the future crossing provider. The declaration's
`source_kind='Intersection'` and `source_id` retain that provenance. No new
Street, Section or junction identity is invented, and no crossing surface is
claimed by generating the taper. Native carriageway junctions own their own
dimensional joins and do not receive this longitudinal declaration.

This policy applies only to consecutive connected Sections at the recorded
attachment location. The initial limits are 1.5 m minimum half-length and 0.3
maximum lateral dimensional change per metre, using the existing quintic law.
GN evaluates the available length and rejects overlapping or infeasible tapers.
Equal-dimension attachments need no new declaration. Width changes at other
undeclared Section seams still require a managed transition declaration.
Transition IDs in native tables are distinct from original item IDs, whose
namespaces can overlap. The original kind and ID remain separate provenance.

Parallel Bundle members may have different widths at a junction. The shared
plateau follows the actual separator: its center can differ from the envelope
center, and each carriageway retains its own parabolic half-profile. This updates
live when widths or signed offsets change. Revision 8.6 checks include reversed
member directions and R1/R2/R4. Revision 8.7 adds exact oblique cuts for complete
straight nonparallel members, preserving physical widths, corner radii and zero
outer-edge height. Both ends must be compatible attached mouths or explicitly
clipped boundaries; arbitrary free endpoints are not moved silently. Revision
8.8 permits curved interiors with certified straight terminal leads while
retaining every authored control. Incompatible cuts or insufficient straight
room remain diagnosed. A real divided intersection is independently qualified;
this is a bounded selection, not unrestricted Bundle support. The test runner's
`--open-bundle-clearance` is an explicit boundary policy and defaults to zero.

The selected connected three-junction test passes R1/R2/R4 at exaggerations 1
and 10. Patch seams now use integer chart correspondence: nearby interior
vertices cannot collapse merely because the cells are small. A physical-position
agreement check remains active before a shared identity is accepted. R4 produces
157,568 quads with two closed boundary loops and no nonmanifold boundary vertices.
Grouped, per-Street and mixed imported inputs preserve the same output; native
vertex edits and managed reversal undo/redo pass. The earlier Street 36
divided-mouth edit limitation remains unresolved. No automatic radius reduction
is used.

## Bundle validation

Revision 8.12 adds a linear-time perimeter certificate for ordered, convex
Bundle cross-sections. Cases that cannot satisfy this sufficient proof use the
existing exhaustive perimeter check. Geometry sampling and all other validity
checks are preserved; every one of the 123 combined-gallery meshes retains its
exact coordinates and topology. On the measured curved width-transition pair,
mean native edits at R1/R2/R4 take 0.052 / 0.060 / 0.187 seconds. These are
fixture measurements, not full-map performance guarantees.

Revision 8.13 extends generated straight leads when an explicitly clipped free
Bundle needs more room for its common cut, retaining every authored control.
The boundary lookup preserves authored knot identities and canonical cut stations
within bounded float precision. The linear perimeter proof has its own larger
work budget; the exhaustive fallback retains its old limit. A real changing-width
Bundle passes R1/R2/R4 and native edit/recovery tests. Its deliberate boundary
cuts still mean partial coverage. General curved mouths and short approaches
remain limited.

## Extend and rebuild

`StreetRenderer` delegates to the selected renderer. `GeometryNodesStreetRenderer`
creates an import context; `StreetFeature` implementations consume its prepared
snapshot and shared asset service. The carriageway provider owns connected-component
selection, explicit unresolved boundaries and native output evidence. Future
sidewalk/marking/rail providers can reuse that context without owning asset loading.

Canonical GN builders/runtime live in `prochitecture_gn/experiments/streets/`.
`blosm/way/gn` is a generated, qualified-import distribution: edit the canonical
source and rerun `experiments/streets/blosm_integration/package_runtime.py`.
Add `--check` with the same `--blosm <addon directory>` argument to verify every
generated module and its source hash without changing files. A nonzero exit
reports missing or out-of-date modules. This checks packaging consistency;
the generated library and installed behavior still require qualification.
Build the library in a fresh background Blender process with `build_assets.py`.
Use `--blosm <addon directory> --original-spans` to build the current original-span
library using the catalog revision. Closed separators, native providers, local
failure containment and standalone Streets are included by default. Omitting
current capabilities requires an explicit comparison revision; qualify such a
build at a separate output path before replacing a bundled asset.
Do not build it on every import or depend on a saved fixture scene.

The integration runners, tests and combined inspection scene are under
`prochitecture_gn/experiments/streets/blosm_integration/`. The development plan,
progress record and illustrated guide are under `prochitecture_gn/docs/streets/`.
Generated Blender files and JSON reports are local artifacts and are not staged
automatically.


Revision 8.14 adds a curvature-checked shorter generated lead when the original
lead leaves insufficient bending room after an explicitly clipped free Bundle
cut. The retry uses conservative per-span width/offset bounds; it excludes a
distant taper only when a chord-distance lower bound proves it cannot reach the
span. Existing valid curves keep their original construction. Source vertices,
physical widths and radii remain unchanged; failed retries still report errors.
This is bounded free-terminal support, not unrestricted curved junction mouths.

The preparation adapter additionally records explicit Bundle member spans in
each Bundle's `topology` report, retaining source Section IDs and internal
junction/attachment references. `MAPPED` describes addressing only; it does not
mean the native solver supports that surface. Internal Bundle junctions remain
rejected until their member-span patches are implemented. Prepared tail ordering
uses the actual tail locations, and minor-attachment owner corrections are
recorded in `ownership_audit.remapped_attachments`. These adapter changes do not
change the 8.14 GN library or the existing saved inspection scene.

`BundleJunctionPlanner` additionally exports `junction_patches` with semantic
inner-edge roles and adjacent span incidences. `PLANNED` is a local capability
record, not generated geometry. A candidate half-profile kernel is qualified
in the development checkout; it is not yet packaged into the 8.14 asset. Its
native separator bridge is qualified locally against independently sampled
straight opposite Sections. A subsequent candidate uses native curved spans
from original Section controls, with uniform quad rows and a checked boundary
approximation; four prepared local assemblies pass at R1/R2/R4. Complete
member-span assembly across Section boundaries and public
registry bindings remain necessary before enabling internal Bundle junctions
in the renderer. A further native candidate joins complete incident Sections
to those half-profile junctions and extends the separator along both approaches.
Three connected synthetic arrangements pass R1/R2/R4 and saved-file native
editing. These candidate groups are not yet in the 8.14 library. Local candidate
meshes do not count as imported network coverage.

The adapter also exports `incident_networks`: supported member-side patches
grouped through shared semantic Sections, with one Section record and independent
start/end junction references. Unbound here means outside that component, not
necessarily a free end in the complete network. A canonical candidate compiler
and shared native group now qualify seven components at R1/R2/R4, including a
three-junction network with two doubly attached Sections. Saved Edit Mode moves
work without street Python. Full Bundle extents, opposite members, consecutive
Section seams and terminal connections still require provider integration;
these records do not bypass the existing unsupported-Bundle report.

The optional span-aware planner additionally exports `member_networks`, with
explicit consecutive Section seams and original authored order. A candidate
network joins those rows and uses declared finite width/offset transitions.
Six real member components pass R1/R2/R4 and saved-file Edit Mode checks, including
four connected intersections with ten once-owned Sections. This candidate is
not yet packaged; dense performance remains open, and opposite members, separator and paired terminal mouths must
be connected through the public provider before enabling full Bundle coverage.

Street editing registration now defers its initial source index while Blender
exposes restricted scene data during add-on startup. The one-time initialization
is removed on unregister; subsequent file loading retains the existing explicit
reindexing. This runtime fix does not change the GN library.

## Retained physical source metadata

Physical aggregation can remove an original Bundle from the active logical
inventory while retaining its original Streets. The adapter now exports those
Bundles separately as `retained_physical_bundles`, using the original Bundle
objects referenced by retained Streets. Member order, head/tail identities and
locations remain available for physical-mouth profile preparation. No source
IDs are replaced and membership is not inferred from spatial proximity.

`PhysicalStreetInventory.from_report` restores both retained Street and retained
Bundle records into a physical preparation view. Duplicate active/retained IDs
are rejected; incomplete or inconsistent member records remain explicit
diagnostics. Retained Bundle coverage is `PENDING_PHYSICAL_PROVIDER` until an
actual provider supplies its surfaces. This metadata does not enable unsupported
Bundle rendering or alter the packaged 8.14 asset.
