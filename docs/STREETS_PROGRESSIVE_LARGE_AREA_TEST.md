# Progressive street import on larger OSM extracts

Assessment: **13 November 2026**, production code at
`386abdb15808ce397149b6c2c951b6bb65e6c87f`, in the isolated
`C:/dev/workspace/blosm_streets_dev` checkout. The production importer was not
changed for these tests. Test instrumentation adds a configurable observation
interval/time limit and exports the importer report even when coverage fails.

## Conclusion

The current implementation **cannot yet be qualified as a satisfactory general
large-area street import**. The scheduler continued returning to Blender, with
no multi-second freeze in these three runs. However, all three imports produced
partial coverage: a few unsupported or invalid approaches blocked large connected
components containing most of the prepared Sections.

These short runtimes do not establish the cost or responsiveness of generating
the complete larger networks. Most of that construction was rejected before it
could run. Earlier successful Karl-Marx-Allee measurements remain valid for that
extract; they cannot be extrapolated to these larger areas.

## Inputs and configuration

The selection includes a dense city center, a larger airport/surrounding road
area and a smaller town area. Full original extracts were imported, without
removing difficult streets or modifying the OSM files.

| Extract under `C:/dev/workspace/osm_extracts` | Approximate bounds area | Carriageway-tagged OSM ways |
| --- | ---: | ---: |
| `facade_visibility/beijing_center.osm` | 6.41 km² | 520 |
| `aviation/airport_bergamo.osm` | 10.49 km² | 1,282 |
| `streets/piestany_north.osm` | 2.26 km² | 222 |
| Previous reference: `facade_visibility/berlin_karl_marx_allee.osm` | 0.50 km² | 146 |

Areas use the file's bounds rectangle with longitude corrected at its middle
latitude. Way counts include `motorway`, `trunk`, `primary`, `secondary`,
`tertiary`, their links, `unclassified`, `residential`, `living_street` and
`service`. A way count is input context, not the number of final Blender objects
or independent construction units.

Each run used:

- Blender 5.2.0 LTS, a fresh factory process, isolated chat profile and disabled
  automatic execution; no normal Blender session or installation was changed.
- The actual production GUI timer and progressive street preview implementation.
- Realistic 3D mode, Python streets only; buildings, forests, railways, water and
  aviation disabled. Airport roads were tested, not runway/taxiway generation.
- Cell size 0.7 m, resolution 1, exaggeration 1, profile rise 0.08 m, default
  smooth interpolation, no terrain and the normal shared street materials.
- The scoped Import performance request enabled, four Blender threads, one
  numerical-library thread, a 3 GiB process cap and 1 GiB system-memory reserve.
- Sequential runs, local OSM input and missing-member downloading disabled.
- A 50 ms observer timer, without simulated continuous viewport navigation.
  This interval only changes the test observer, not the import scheduler.

The separate GUI process was hidden. The heartbeat measures opportunities for
Blender to run UI callbacks and whether output objects are visible in its active
view layer. This is not a direct measurement of on-screen frames or a manual
navigation/Cancel trial on the full larger network.

## Measured results

| Extract | First visible preview | Import elapsed | Longest UI heartbeat gap | Components READY / total | Accepted Sections / prepared Sections |
| --- | ---: | ---: | ---: | ---: | ---: |
| Beijing | 8.37 s | 13.83 s | 0.66 s | 8 / 11 | **14 / 1,139** |
| Bergamo | 22.57 s | 55.87 s | 1.19 s | 22 / 75 | **31 / 2,264** |
| Piešťany north | 6.63 s | 9.46 s | 0.53 s | 5 / 8 | **9 / 526** |

Prepared Sections are those represented in the attempted numerical components.
The accepted fraction is approximately 1.2%, 1.4% and 1.7%, respectively. This is
a count of Section records, **not** road-length, road-area or all-OSM-way coverage.
Component counts alone disguise the severity because one rejected component can
contain far more streets than many small accepted components combined.

| Extract | 95th-percentile heartbeat gap | Longest import callback | Peak resident / sampled private memory | Published quads |
| --- | ---: | ---: | ---: | ---: |
| Beijing | 0.195 s | 0.594 s | 1.036 / 1.380 GiB | 17,828 |
| Bergamo | 0.273 s | 1.175 s | 1.249 / 1.612 GiB | 87,655 |
| Piešťany north | 0.133 s | 0.485 s | 0.897 / 1.339 GiB | 11,086 |

There were no heartbeat gaps over five seconds; only Bergamo had a gap over one
second. In all three cases, the longest callback was in the initial OSM feature
preparation/crossing phase. Each process stayed within its memory cap and exited
normally. Nevertheless, the strict GUI qualification is **failed**, because its
required complete import ended as `partial`. Process completion is not a passing
coverage result. Memory measurements also cover partial output, not the missing
complete meshes.

## Failure causes

The Beijing blocking/cooperative comparison passed in
`coverage_parity_03/checks.json`: both modes produced the same outcomes, reported
the same rejected source OSM ways, and published exactly matching accepted
geometry. The three Beijing component failures therefore also occur without
progressive previews. This comparison was not repeated for the other two extracts.

### A blocked approach rejects its whole connected component

In `PreparedPythonAdapter.read`, grade/elevation exclusions, missing required
approaches and invalid Bundle members can mark otherwise eligible Sections as
blocked. `components()` groups by the complete declared junction/connection and
Bundle topology. `PythonStreetRenderer._render_steps` rejects a component if its
Section set intersects `adaptation.blocked_sections`.

This prevents silently treating an unsupported junction approach as a valid free
road end. It also means a local unsupported case can suppress hundreds of other
ground-level Sections in the same connected network.

Observed adapter issues:

| Extract | Grade/elevation exclusions | Missing required approaches | Excluded/invalid Bundle membership |
| --- | ---: | ---: | ---: |
| Beijing | 11 | 3 | 0 |
| Bergamo | 175 | 103 | 3 |
| Piešťany north | 3 | 3 | 0 |

These are counts of reported issue records, not the number of rejected components.
The largest rejected components contain 1,122 Sections in Beijing, 824 in Bergamo
(with further large rejected components of 688 and 320), and 513 in Piešťany.

### Other numerical/topological rejections

Beijing also reports `Repeated task dependency` for a two-Section service-road
component. One source is a closed loop, OSM way `724548765`, attached to service
way `724548762`. This needs a dedicated closed-loop dependency investigation;
the current assessment does not bypass the dependency guard.

Bergamo additionally reports failures to fit rounded corner contacts, insufficient
source length for smooth bends, a potential smooth-offset fold, and incomplete
original incidence support at an outside-boundary junction. Preview publication
does not make those rejected surfaces valid.

## Next work suggested by the evidence

1. Define safe local boundaries around unsupported features so one excluded
   approach does not automatically suppress an entire connected urban network.
   Preserve explicit unsupported-junction diagnostics and valid seam/profile
   contracts; simply dropping an approach or disabling validation is insufficient.
2. Correct the closed-loop dependency case and add its small OSM-derived
   regression fixture.
3. Improve corner fitting, short-bend handling and retained boundary incidence
   using the reported source identities as reproducible cases.
4. Repeat full-coverage large-area GUI measurements after those changes. Only
   then assess the memory and longest-builder pauses for a truly generated large
   connected network, followed by a manual navigation/Cancel trial.

The preview scheduling mechanism is useful, but broad construction support and
failure isolation are the immediate limitation exposed by this assessment.

## Evidence and reproduction

All generated evidence is under:

`C:/dev/workspace/artifacts/streets/progressive_large_area_20261113/`

The `beijing_01`, `bergamo_01` and `piestany_north_01` folders contain
`ui_checks.json`, progress snapshots and `process/` resource reports/logs.
`ui_checks.json` includes the original source SHA-256, runtime source hashes,
heartbeat history, component outcomes, complete importer diagnostics and counts.
`summary.json` contains the reduced measurements and failure breakdown.

Use `tests/check_street_progressive_ui.py` in a resource-guarded isolated Blender
GUI, with `--artifacts <run-folder> --osm <original-extract>
--performance-request --no-navigation --timeout-seconds 1800
--heartbeat-seconds .05`. Keep the production code at the recorded revision to
reproduce these observations. The script deliberately records a failed check
when output coverage is incomplete.

Regenerate the aggregate summary with a normal Python interpreter:

```powershell
python tests/summarize_street_progressive.py --artifacts C:/dev/workspace/artifacts/streets/progressive_large_area_20261113
```

The additional `test_real_extracts_coverage_parity` in
`tests/check_street_progressive.py` compares blocking and cooperative imports,
including accepted geometry and failed source OSM identities. It supports partial
coverage as an explicit outcome; passing parity does not qualify that coverage
as satisfactory.

The Beijing parity check passed in a separate bounded background Blender process
(12.58 s for the complete test, peak resident memory 0.481 GiB). Initial harness
attempts were corrected to select the current collection's report and compare
stable OSM way identities rather than transient Section IDs. Their failed
assertions remain in the earlier `coverage_parity_01`/`02` logs; the passing
record is `coverage_parity_03`.
