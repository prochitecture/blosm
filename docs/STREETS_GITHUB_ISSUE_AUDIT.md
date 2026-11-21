# GitHub street issues versus current buildings_gn

Audit date: **13 November 2026**. Audited production revision:
`386abdb15808ce397149b6c2c951b6bb65e6c87f` (`buildings_gn` and this chat's
`streets_development` initially matched). Later uncommitted large-area test and
documentation work is preserved; it does not alter the production importer.

**Implementation follow-up:** the findings below describe that audited revision.
Subsequent fixes and whole-import qualification are tracked in
[Beijing and Bergamo import qualification](STREETS_BEIJING_BERGAMO_FIXES.md).
Focused fixes on the development branch do not by themselves establish that
either complete city import is qualified.

## Findings

The issue tracker contains useful evidence about the preparation layer. However,
**an open GitHub issue does not establish a current defect**, and a closed issue
does not establish that every related geometry case is supported.

The audit found two confirmed preparation defects and several current adapter or
mesh-construction limitations. The old Berlin ownership/list bugs did not recur
in the inspected inputs. The most important remaining problems are:

| Finding | Layer | Current evidence | Priority |
| --- | --- | --- | --- |
| F1. XY crossings can create connections without shared OSM nodes, including grade-separated routes | Shared OSM preparation | Reproduced in Beijing, Bergamo, Tokyo and other extracts | High |
| F2. A few unsupported approaches suppress a large connected network | Python adapter/publication scope | Previous full imports: 14/1,139 Sections accepted in Beijing; 31/2,264 in Bergamo | High |
| F3. A returning/closed Section duplicates a numerical task dependency | Python task declaration | Duplicate profile dependencies reproduced in Beijing, Friedrichstrasse and Bergamo | High |
| F4. Automatic traffic-side detection selects right-hand traffic in Tokyo | Shared preparation/configuration | Tokyo sample reports `right_hand_traffic=True`; Japan is missing from the lookup | High |
| F5. Some prepared Bundle configurations cannot be mapped to the current surface contract | Shared-preparation-to-backend boundary | End-order changes, unresolved membership, and missing opposite-member crossing ports | Medium/high |
| F6. Rounded-contact, short-bend and boundary-support construction rejects additional inputs | Python geometry/source support | Reproduced failures in the prior full Bergamo import; further full-import checks below | Medium/high |
| F7. Width correction crosses a genuine minor-junction lane/width change | Shared source preparation | Confirmed in Bergamo during the November 15 mesh investigation; nominal-dimension boundary fix passes focused/native preparation checks | High |

F2, F5 and F6 include **unsupported scope and conservative rejection**, not proof
that every rejected input has a valid solution under the current constraints.
The report distinguishes that from the confirmed topology, traffic-side and
dependency defects. No production fix is included in this audit.

## Method and evidence limits

The GitHub inventory includes **123 issues** (67 open, 56 closed), fetched from
both open and closed collections. Titles/bodies were screened across the whole
inventory; **1,663 comments from 90 street-related discussions** were retrieved
to check fixes, changed requirements and sample references. Large threads were
paginated. The [issue-by-issue matrix](STREETS_GITHUB_ISSUE_MATRIX.md) records
relevance and the current disposition for all 123 retrieved issues.

Source evidence comes from the current local branch, not the repository's
default branch or the historical source links in comments. Source issue IDs and
links are retained, but transient Street/Section IDs from old posts are not
assumed to match a newly prepared network. Reproductions use the original local
OSM files and retain OSM way identities where possible.

The preparation audit ran **26 configurations over 22 distinct extracts** in
isolated factory Blender 5.2 processes. Kuzmice and Toronto were also run with
railways enabled. It checks the actual prepared objects before the Python mesh
backend, including:

- Active Street item ownership and unique Section ownership.
- Head/tail reachability, predecessor/successor reciprocity and chain cycles.
- Major connector-ring closure, reciprocity and owner identity.
- Endpoint connector types and availability of referenced Bundle members.
- Width/offset propagation after SideLane/SymLane through minor intersections.
- Existing provenance-based XY connection diagnostics, snapshot scope results,
  numerical adaptation and task dependency declarations.

The preparation-only mode replaces the renderer **inside the isolated test
process** after normal parsing/generation. It does not construct meshes. Its
`success: 0 completed` operator message means the observation hook completed;
it is not full street-generation qualification. The full-mesh mode runs the
normal renderer and records component outcomes separately.

The checks found **zero raw invariant violations** across those 26 preparation
runs, including **202 transition-propagation checks**. This is evidence against
recurrence of the tested ownership/list defects, not a proof about every possible
OSM graph or screenshot-only geometry report.

All generated evidence is under:

`C:/dev/workspace/artifacts/streets/github_issue_audit_20261113/`

Raw GitHub snapshots are `issues.json` and `comments.json`. Each Blender run has
its source hashes, observations and `process/` logs/resource report. Code and
this documentation remain in the repository. The user's main Blosm checkout,
normal Blender profile and installed addon were not modified.

## F1. Geometric crossings are used as topology

**Confirmed preparation defect.** `StreetGenerator.findSelfIntersections` sends
planar segment pairs to the sweep intersector. `createWaySectionNetwork` then
splits/joins using those positions. OSM node identity and grade compatibility do
not gate that graph construction. `ConnectivityAudit` records suspicious events
but explicitly does not repair the graph.

Relevant source:

- [action/generate_streets.py](../action/generate_streets.py),
  `findSelfIntersections` and `createWaySectionNetwork`.
- [way/connectivity_audit.py](../way/connectivity_audit.py), `ConnectivityAudit.read`.
- [way/street_network/connectivity.py](../way/street_network/connectivity.py),
  later surface-connectivity review.

| Original extract | Different-grade XY events | Other inferred XY connections |
| --- | ---: | ---: |
| Beijing center | 20 | 36 |
| Bergamo | 103 | 0 |
| Tokyo Shibuya | 115 | 29 |
| Moscow Leningradsky Prospekt | 8 | 0 |
| Taipei | 6 | 3 |
| Milano | 3 | 1 |
| Karl-Marx-Allee | 0 | 2 |

These counts are **audit events without a shared OSM endpoint node**. Different
tag tuples identify a grade warning, not a derived numeric elevation. Same-grade
events can also reflect mapping defects; the importer must diagnose them rather
than silently claim that a physical connection exists.

Concrete examples:

- Beijing ways `121903964` (bridge footway) and `121903991` (primary road),
  crossing near local XY `(1913.663, -551.419)` without a shared endpoint node.
- Bergamo ways `49419548` and `54680038`, a bridge/ground tertiary/motorway pair,
  near `(933.711, -1018.887)`.

This can alter minor/major intersections, Street aggregation and Bundle
membership before any new mesh formula runs. Filtering out elevated surfaces
later does not undo those topological decisions.

GitHub relationship: the older workflow and topology discussions
[#72](https://github.com/prochitecture/blosm/issues/72),
[#80](https://github.com/prochitecture/blosm/issues/80),
[#104](https://github.com/prochitecture/blosm/issues/104) provide context.
The old Beijing division-by-zero report
[#82](https://github.com/prochitecture/blosm/issues/82) is **not** the same defect:
that exception did not recur in either Beijing or Tokyo preparation.

Suggested work: preserve OSM node/endpoint incidence and grade semantics before
building connectivity; retain explicit diagnostics for malformed source data.
Do not fix this by arbitrary XY tolerances or by deleting all bridge-tagged
segments after they have already affected the graph.

## F2. Component-wide rejection amplifies local unsupported input

**Confirmed current limitation at the adapter/renderer boundary.** The ground
surface policy excludes indicated elevation, bridges and most tunnels. Missing
required approaches and invalid Bundle members add blocked Sections. The Python
renderer rejects a component whenever its Section set intersects that blocked
set.

Relevant source:

- [way/street_network/physical_ports.py](../way/street_network/physical_ports.py),
  `GroundSurfacePolicy.reason`.
- [way/street_network/python_backend/adapter.py](../way/street_network/python_backend/adapter.py),
  `PreparedPythonAdapter.read` and `components`.
- [way/street_network/python_backend/renderer.py](../way/street_network/python_backend/renderer.py),
  `_render_steps`.

| Full import | Accepted/prepared Sections | Largest rejected component |
| --- | ---: | ---: |
| Beijing | 14 / 1,139 | 1,122 Sections |
| Bergamo | 31 / 2,264 | 824 Sections, plus other large rejected components |
| Piešťany north | 9 / 526 | 513 Sections |

These are Section counts, not road-length percentages. The full-import evidence
and blocking/progressive Beijing parity are in the
[large-area assessment](STREETS_PROGRESSIVE_LARGE_AREA_TEST.md). The new
preparation audit independently reproduced blocked sets of 1,123, 2,031 and 517
Sections, respectively. Additional rejected Sections fail for other reasons.

Conservatively refusing an incomplete junction is appropriate. Suppressing most
of an otherwise useful city network is the functional limitation to address.
The dependency concerns in [#85](https://github.com/prochitecture/blosm/issues/85)
are related design context, not a reproduced dictionary-mutation exception.

Suggested work: define explicit supported boundaries around unsupported
incidences, preserve diagnostics and consistent profile/ownership seams, then
render the qualified subnetwork. Simply ignoring a missing approach is unsafe
for correctness. Implementing elevated roads and improving failure isolation are
separate deliverables.

## F3. Closed/returning Sections duplicate task dependencies

**Confirmed Python backend defect; related to, but distinct from, the old
connector-key collision in [#54](https://github.com/prochitecture/blosm/issues/54).**

`IntersectionProvider.tasks` correctly retains separate endpoint incidences, but
also appends a profile dependency once per incidence. When both ends reference
one Section, the same profile key occurs twice. `DependencyIndex` rejects that
task with `Repeated task dependency`.

| Extract | Observed task | Repeated dependency |
| --- | --- | --- |
| Beijing | `physical_intersection:3081` | `profile:1541` |
| Friedrichstrasse | `physical_intersection:107` | `profile:54` |
| Bergamo | `physical_intersection:3599` | `profile:1800` |

These local IDs belong to the recorded audit runs. Beijing's stable source
example is closed service way `724548765`, attached to way `724548762`.
The shared snapshot can already flag a closed Section as requiring dedicated
treatment, but Python conversion still admits it to task declaration.

Sources:
[intersections.py](../way/street_network/python_backend/intersections.py),
`IntersectionProvider.tasks`;
[updates.py](../way/street_network/python_backend/updates.py), `DependencyIndex`;
[shared adapter](../way/street_network/adapter.py), `street_support`;
[Python adapter](../way/street_network/python_backend/adapter.py).

The previous full Beijing import reports this exact exception. In Bergamo it is
also a latent declaration failure; the normal full import may encounter a path
rejection earlier. A preparation declaration failure is not automatically the
first error the full renderer will report.

Suggested work: keep distinct `(Section, end)` incidences while deduplicating
resource dependencies. Then qualify the closed-loop geometry itself. Removing
the duplicate-dependency guard alone would not prove that a valid returning-road
surface can be constructed.

## F4. Traffic-side lookup is incomplete

**Confirmed preparation/configuration defect.**
[#126](https://github.com/prochitecture/blosm/issues/126) requested correction of
the old `leftHandTraffic` naming/default. Its comments describe a committed
rename. That part exists: `WayManager` now owns `rightHandTraffic`.

However, the current automatic lookup uses a small collection of latitude/
longitude rectangles and defaults to right-hand traffic outside them. Japan is
not in the table. The unmodified Tokyo Shibuya import therefore reports
`right_hand_traffic=True`; Melbourne correctly reports `False` in the same audit.
Japan uses left-side traffic; see the
[Japan Automobile Federation's traffic rules](https://english.jaf.or.jp/safe-driving/traffic-rules-in-japan).

Sources:
[way/manager.py](../way/manager.py), `isRightHandTraffic`;
[way/way_properties.py](../way/way_properties.py), `leftHandTrafficAreas` and
`lanePattern`;
[setup/__init__.py](../setup/__init__.py), `getWayManager`.
Setup currently forces the automatic mode, so an override policy also needs
review. The old compatibility name remains inside StreetGenerator, but the
observed Tokyo misclassification is a functional defect beyond naming.

Suggested work: make traffic-side resolution explicit, cover all supported
regions and respect intentional overrides. Add a Tokyo regression case and
tests at the lookup/override boundaries before assessing asymmetric lane/offset
behavior there.

## F5. Bundle topology has unhandled configurations

**Confirmed contract/feature limitations; not proof that every upstream Bundle
is wrong.** Relevant discussions are
[#80](https://github.com/prochitecture/blosm/issues/80),
[#104](https://github.com/prochitecture/blosm/issues/104),
[#113](https://github.com/prochitecture/blosm/issues/113) and
[#115](https://github.com/prochitecture/blosm/issues/115).

The current Bundle topology mapper explicitly diagnoses member order changes,
incomplete membership, ambiguous chain joins and insufficient distinct members.
For example, the existing Bergamo full-import snapshot has seven order-change,
five unresolved-membership and two ambiguous-chain-join diagnostics. Beijing has
three order-change and three unresolved-membership diagnostics.

Task declaration also rejects prepared Taipei and Milano networks with
`Bundle interior Streets need two explicit opposite-member crossing ports`.
These are declared-topology checks, not performance failures. Full-render
outcomes are recorded separately below because retained boundary support can
change what reaches final declaration.

Sources:
[bundle_topology.py](../way/street_network/bundle_topology.py),
`PreparedBundleTopology`;
[Python Bundle adaptation](../way/street_network/python_backend/adapter.py);
[Bundle crossing code](../way/street_network/python_backend/bundle_regions.py).

The original Piešťany bug [#114](https://github.com/prochitecture/blosm/issues/114)
was reported fixed. Its current source passes the structural checks, but a
matching linked-list structure alone cannot verify every yellow road shown in
the old screenshot. That visual membership claim remains distinct from the
confirmed present-day mapper limitations.

Suggested work: use source-backed examples to distinguish an incorrect
preparation group/order from a valid configuration requiring an expanded backend
contract. Do not reorder or pair members by proximity just to satisfy a solver.

## F6. Geometry and boundary-support rejection remains

The previous full Bergamo import also rejects rounded-contact fitting, short
smooth bends, possible offset folding and incomplete outside-boundary incidence
support. These remain actual unsupported/rejected inputs in the current branch.
See [the large-area report](STREETS_PROGRESSIVE_LARGE_AREA_TEST.md) for the exact
messages and affected component counts.

Related older topics include
[#65](https://github.com/prochitecture/blosm/issues/65) (overlaps),
[#83](https://github.com/prochitecture/blosm/issues/83) (fillets),
[#109](https://github.com/prochitecture/blosm/issues/109) and
[#120](https://github.com/prochitecture/blosm/issues/120) (nearby junction groups).
The current backend already has compound-junction planning; these issues cannot
be summarized as “no grouping exists.” Nor should a present-day curved-contact
rejection be labeled as the identical old polygon bug without matching evidence.

The offset-collapse question [#111](https://github.com/prochitecture/blosm/issues/111)
does not supply a precise local reproduction file/coordinate pair. Its exact
outer-edge event was not independently reproduced. The current smooth-offset
fold rejection is related geometric risk, not confirmation of that same event.

## Historical fixes that should not be reintroduced as current bugs

The relevant comments explain several fixes even while issues remain open:

- [#122](https://github.com/prochitecture/blosm/issues/122): minor Intersection
  predecessor/successor behavior was patched; the issue stayed open for a cleaner
  class structure. `Intersection.pred/succ` now derive the through-Street
  Sections. Reciprocal-list checks passed on the audited inputs.
- [#121](https://github.com/prochitecture/blosm/issues/121): width/offset
  propagation and configurable pixel correction were implemented. Current
  `StreetGenerator.postProcessing` propagates dimensions through minor
  intersections; 202 propagation checks passed. Pixel correction defaults off
  for the current geometric backend. A later related over-propagation defect
  is described in the November 15 section below; the original check did not
  exercise a genuine source lane-count boundary.
- [#105](https://github.com/prochitecture/blosm/issues/105),
  [#110](https://github.com/prochitecture/blosm/issues/110),
  [#123](https://github.com/prochitecture/blosm/issues/123): no active-chain owner
  mismatch or duplicate active Section ownership was found.
- [#124](https://github.com/prochitecture/blosm/issues/124),
  [#125](https://github.com/prochitecture/blosm/issues/125),
  [#106](https://github.com/prochitecture/blosm/issues/106): no unavailable
  connector Street/Bundle terminal member was found. The current backend also
  has an explicit retained-physical-source inventory for absorbed logical roles;
  a retained source is not automatically a ghost Street.
- [#117](https://github.com/prochitecture/blosm/issues/117),
  [#118](https://github.com/prochitecture/blosm/issues/118): the tested endpoint
  types and ring closure/reciprocity are valid, including rail-enabled inputs.
- [#82](https://github.com/prochitecture/blosm/issues/82),
  [#115](https://github.com/prochitecture/blosm/issues/115): the original division
  errors did not occur in current preparation of Beijing, Tokyo, Kuzmice or
  Toronto. Rail-enabled preparation is not a claim that the Python carriageway
  backend now generates rail or tram meshes.

Many older reports reference `IntersectionArea`, `TrimmedWaySection`,
`TransitionSideLane`, graph-cycle/PyGEOS triangulation or the original GN asset
renderer. Those are not the active Python street-rendering pipeline. Shared
geometric concepts remain relevant; obsolete private attributes and renderer
exceptions should not be presented as current bugs merely because a ticket is
still open.

## Reproduction

Run `tests/check_street_issue_audit.py` with the recorded production revision and
this chat's isolated factory profile. Every launch must retain
`--factory-startup --disable-autoexec` and a resource guard. After Blender's `--`,
pass `--artifacts <output-directory>` and repeat `--osm <original-file>`.

- Default: production parsing/preparation and adapter/task-declaration audit;
  no mesh construction. `--railways` includes rail/tram input.
- `--full-mesh`: normal renderer and import reports; no preparation replacement.

Run folders `preparation_01`, `preparation_02` and `rail_preparation_01` contain
the 26 preparation configurations. Resource peaks were below 0.4 GiB resident
for those runs under the 3 GiB guard and 1 GiB system reserve. Sample names,
source SHA-256 and detailed findings are in each `preparation_audit.json`.

Do not infer correctness from an old transient numeric object ID, an open/closed
GitHub state, a successful preparation-only operator message or the absence of
an exception. Full output coverage, geometry and source-backed incidence need
their own acceptance checks.

## Later confirmed dimensional propagation defect (November 15)

The initial #121 check proved that width/offset propagation existed; it did not
prove that the propagation stopped at every legitimate dimensional boundary.
The Bergamo full-mesh investigation finds an actual related defect: a two-lane
source road (`304685645`) reaches the renderer at 10.5 m, inherited from an
earlier three-lane Section through a minor side junction. Its original lane
count and category/style resolve to 7 m. The excessive width overlaps an
otherwise separate nearby junction. No geometric Bundle inference is needed
to explain or correct this mismatch.

`Section.sourceWidth` now records the original resolved width before optional
pixel correction. `TransitionDimensionPropagation` transfers a prepared
correction only through the same original width/lane-count run. A true width
or lane-count change ends the propagation; its existing native junction owns
the transition. Equal-source-width runs still inherit the original pixel and,
for SideLane, offset corrections. No source centerline, OSM incidence, lane
count, width tag or established successful correction is rewritten.

Twenty-eight native source/transition/Bundle regressions pass. Fresh native
preparation of both original extracts reports no structural findings. Comparing
all supported inputs finds three corrected Bergamo widths (10.5 to 7 m, 3.5 to
7 m and 7 to 3.5 m), and no Beijing dimensional change. These are source-stage
results; the complete corrected imports, meshes and global overlaps remain
separate qualification gates. Evidence is `preparation_16_source_widths` and
`dimension_source_regressions01` under the current task artifact directory.

## Additional full-import results (initial audit)

Six original extracts were also passed through the **unmodified normal mesh
renderer**, after the preparation-only audit. These are functional observations,
not timing benchmarks: the blocking test does not enable the GUI performance
request. All runs used the same isolated profile and resource guard.

| Extract | Full import outcome | Main observed result |
| --- | --- | --- |
| Taipei | Partial: 3/5 components READY | 237-Section component rejects shared profile for incomplete lateral members; another component exceeds the stored quad edge-ratio limit |
| Friedrichstrasse | Partial: 1/2 READY | 102-Section component rejects a closed Section requiring a dedicated loop provider |
| Versmannstrasse | Failed: 0/2 READY | Rounded contact lacks full-profile entrance space; another component is blocked by unsupported source context |
| Nabrezie | Failed: 0/1 READY | Compound clipping retains a partial incoming profile |
| Nashville Wyoming Avenue | Complete: 2/2 READY | All 138 prepared Sections generated; 254,547 quads across the two components |
| Milano | Partial: 1/3 READY | 334-Section component rejects insufficient length for a smooth bend; another component exceeds the stored quad edge-ratio limit |

The successful Nashville result is direct current evidence that close-junction
grouping is implemented for the sample in #120. The current Nabrezie and
Versmannstrasse failures are reproducible but do not establish the exact old
visual symptoms of #68/#69: the old screenshots cannot be requalified while
the current renderer rejects those components.

Evidence is in `full_mesh_01/preparation_audit.json` and
`full_mesh_02/preparation_audit.json`. The second run captures the terminal
reports for the two fully failed imports, for which Blender raises an operator
error. An expected recorded import failure does not count as successful mesh
qualification merely because the evidence-collection process exits normally.
