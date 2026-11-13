# Python street generation: profiling and implementation plan

Measured on 11 November 2026, runtime revision
`b9517de05889a3a7f0d7bc8107d923f7af871f41`, Python backend `0.1.0.dev78`.
This document supersedes the earlier timing investigation as the optimization
plan. [STREETS_PERFORMANCE.md](STREETS_PERFORMANCE.md) retains historical measurements;
[IMPORT_PERFORMANCE.md](IMPORT_PERFORMANCE.md) describes the already implemented
Windows scheduling preference.

## Objective and conclusion

The requested stretch target is a **tenfold reduction in generation time**:
approximately **72 seconds to 7 seconds** on the current Karl-Marx-Allee test.
It applies to a fresh generation, not only a cached Update. The implementation
must work from geometry and input contracts, without rules for particular OSM
files, Section IDs or junction IDs.

**Preserving the visual properties of the generated streets and intersections
is a hard acceptance requirement, and takes precedence over the 10× target.**
Optimization may change internal computation or mesh connectivity, but must
preserve the visible surface shape, smoothness, silhouettes, terrain fit,
separator appearance and material/texture mapping. A faster result that visibly
changes these properties does not qualify, even if it satisfies the existing
minimum geometric quality limits.

The profile gives concrete opportunities, but **does not yet demonstrate that
10× is achievable** at the same settings and quality. Small allocation/caching
changes are useful first deliveries; reaching 10× requires a substantially
cheaper construction path for expensive junctions/separators, together with
faster terrain queries and validation. Existing generic construction must remain
available for inputs that do not qualify for a faster path.

The baseline profiling task changed only benchmark/inspection tools and
documentation. Implementation began afterward; measured changes and the current
resume point are in [STREETS_PERFORMANCE_IMPLEMENTATION.md](STREETS_PERFORMANCE_IMPLEMENTATION.md).
The tenfold target has not yet been achieved. Development stays
in this chat's `blosm_streets_dev` checkout; the user's main checkout and normal
Blender profile remain untouched.

### Implementation checkpoint — 12 November 2026

The initial optimization pass is implemented in runtime `0.1.0.dev81`.
Final alternating reference/candidate controls measured **70.278 → 41.768
seconds**, an observed **1.68× speedup / 40.6% time reduction**. An earlier five-run
candidate batch measured 37.418 seconds median, compared with the original
72.133-second baseline. The [implementation report](STREETS_PERFORMANCE_IMPLEMENTATION.md)
keeps both batches and explains the timing drift; it uses the more conservative
alternating control as the headline.

The final saved Berlin output, and additional Canal Street and Nashville
outputs, retain exact geometry, normals, UVs and material assignments. All 983
numerical regressions and native lifecycle/Update checks pass. All eight
controlled Berlin render pairs are pixel-identical. The report records each
stage's selected changes, unselected experiments, cache/Update qualification,
runtime hashes, artifacts and the remaining construction bottlenecks. The
tenfold target is still unmet; no lower resolution or weaker checks are used.

## 1. Measurement method and baseline

Three separate factory-startup Blender processes executed the normal **3D
realistic Import**, with streets as the only imported feature. They read the
effective settings and original terrain from
`C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend`. The production
performance preference was enabled. No affinity or working-set throttle was
applied. The OS file cache was not flushed.

| Run | Street generation | Whole Import action | Numerical preparation | Publication, including terrain |
| --- | ---: | ---: | ---: | ---: |
| 1 | 70.926 s | 71.418 s | 56.567 s | 11.431 s |
| 2 | 72.133 s | 72.663 s | 57.005 s | 12.313 s |
| 3 | 73.616 s | 74.111 s | 58.365 s | 12.319 s |

Median street generation: **72.133 s**. Median whole Import: **72.663 s**.
Three observations establish a practical baseline, not a statistical latency
guarantee. Saving, Blender startup, addon registration, benchmark assertions,
rendering and interactive UI redraw are outside the timed Import interval.
The generator's own production validation is inside it.

Conditions:

- Blender 5.2.0 LTS, build `fbe6228777e7`; bundled NumPy 2.3.4.
- Intel Core i5-12500H; 16 Blender threads, one numerical-library thread.
- Resolution 1; effective cell size 0.6999999881 m; smooth interpolation;
  crown rise 0.0799999982 m; exaggeration 1.
- Python profile-preserving terrain placement, sample spacing 1 m, clearance
  0.02 m, maximum reference grade 10%, maximum lift 25 m.
- Library-default pavement/grass materials. This is not a textured render test.
- 3 GiB child memory cap and 1 GiB system reserve; all runs completed normally.
  Maximum resident memory was **601 MiB**, maximum sampled private memory
  **673 MiB**. Extra memory is not the first optimization to pursue.
- Input: `osm_extracts/facade_visibility/berlin_karl_marx_allee.osm`, SHA-256
  `d22cf78ef0bd8daa2473dea0a85687f5483f1bd8567c2d7d8deebf9e6c2583dc`.
- Source scene SHA-256:
  `4a2865752b8be826f01e1944a8baadd195026030a455e9faded554dc2ec57e17`.

All runs reported COMPLETE coverage, two networks, 328 output meshes and
350,571 quads. Independent inspection compared the saved third ordinary run
against the separate profiled run: normalized sources, prepared and final
float32 positions, face indices, UVs, material names and construction metadata
matched exactly. Matching counts alone was not used as that proof.

### Where the ordinary run spends time

The following rows partition **run 2**, the median-duration street run. They
are non-overlapping; using one actual run avoids adding unrelated medians.

| Work | Time | Share of street generation |
| --- | ---: | ---: |
| Junction surfaces | 30.895 s | 42.8% |
| Separators | 16.244 s | 22.5% |
| Remaining numerical work, including overlap checks | 9.865 s | 13.7% |
| Terrain preparation and fitting | 7.408 s | 10.3% |
| Remaining publication, including stored-precision validation | 4.905 s | 6.8% |
| Adaptation, source objects, materials and other renderer work | 2.815 s | 3.9% |
| **Street generation** | **72.133 s** | **100%** |

OSM preparation and the remaining operator work add approximately 0.53 s.
Initial path preparation takes 0.560 s and ordinary road mesh tasks 0.304 s
in this run. However, later construction repeatedly queries path frames;
those calls are charged to their calling surface tasks, not to path preparation.

### Workload size

| Output kind | Objects | Quads |
| --- | ---: | ---: |
| Ordinary road parts | 237 | 238,720 |
| Junction pavement | 70 | 66,699 |
| Separators | 9 | 32,357 |
| Bundle pavement | 2 | 7,360 |
| Section seam surfaces | 9 | 5,248 |
| Junction grass | 1 | 187 |
| **Total** | **328** | **350,571** |

Ordinary roads contain **68% of the quads**, although their direct construction
is inexpensive. They therefore matter for downstream terrain sampling and
overlap validation. Optimizing only intersection face counts misses most of
the scene's faces.

## 2. Function-level profile

A separate run recorded numerical and publication profiles independently.
It took 94.253 s for generation, including instrumentation overhead, and is
**not** a fourth baseline sample. The numerical profile contains 67.76 million
calls; the publication profile contains 6.36 million.

Times below include called functions. **Nested rows must not be added** or
treated as predicted unprofiled savings.

| Routine / operation | Calls | Profiled cumulative time | Implication |
| --- | ---: | ---: | --- |
| `CompoundQuadChart.build` | 43 | 29.530 s | Generic compound partition and chart construction dominate many junctions |
| `FixedRimSeparator.build` | 4 | 21.716 s | Four difficult separator regions account for almost all separator work |
| `visibility_center` | 1,440 | 17.373 s | Repeated half-plane clipping uses many tiny Python/NumPy operations |
| `VisibilitySplitter.split` | 125 | 15.945 s | Candidate polygon cuts repeatedly invoke the kernel calculation |
| `FixedRimQuality.improve` | 164 | 13.784 s | Repeated shape repair and line searches are costly |
| `FourSidedQuadChart.build` | 11 | 13.579 s | Multiple atlas proposals and their repairs amplify work |
| `PatternLibrary.quad` | 155 | 8.035 s | Cold pattern construction is significant despite an existing cache |
| `Path.frame` | 9,399 | 5.853 s | Repeated downstream station/frame evaluation deserves batching/reuse |
| Numerical `validate_overlaps` | 4 | 4.578 s | Spatial validation is significant even before publication |
| `SlopeProjection.solve` | 45 | 3.935 s | Bounded iterative slope correction is another repair cost |
| `BoundaryJetField` setup / evaluation | 43 / 43 | 1.740 / 2.244 s | Height-field work is measurable, but not the first dominant target |
| Terrain sampler | 1,238 batches | 6.310 s | Python dispatch around millions of rays is costly |
| Stored-precision `validate_overlaps` | 2 | 4.287 s | Publication is largely validation, not Blender datablock creation |
| Blender mesh creation, `SceneWriter._mesh` | 328 | 0.388 s | Bulk mesh publication already works well |

Important details:

- There were **3,541,148 terrain rays**. Only 1.604 profiled seconds were inside
  the BVH `ray_cast` method; the sampler's own Python work took 4.705 s.
  Per-face sampling includes at least nine points, and the final lift check
  queries mesh vertices again.
- Four-sided mapping calls account for 7.559 of the 8.035 s in `quad`.
  Its longest-boundary-edge search performs **1,463,469 scoring calls**.
  The general NumPy norm routine was called 2.69 million times across the
  numerical backend. Large vector operations are not the only workload here.
- The main network's existing pattern cache recorded 1,102 hits and 390 misses
  and held about 7.87 MiB against its 8 MiB limit. A miss count does not reveal
  which entries were rebuilt after eviction; add that distinction before
  changing the memory budget.
- The accepted-output metadata identifies four partitioned separator charts.
  Runtime failed attempts are included in the profile, but a complete
  attempt-by-attempt rejection ledger is not yet recorded.
- Blender bulk `foreach_set` calls totalled only 0.213 profiled seconds.
  Replacing publication with another backend, merging all objects, or adding
  more Blender threads is not the first response supported by these numbers.

## 3. What 10× requires

Junction and separator tasks total 47.140 s in run 2. Even removing them
entirely would leave **24.992 s**, at best a 2.89× speedup. Accelerating them
alone by 10× would still leave **29.706 s** overall, about 2.43×.

The target must therefore be end-to-end. The following is an **engineering
budget to test**, not an achieved result or a prediction:

| Work package | Stretch budget |
| --- | ---: |
| OSM preparation and other operator work | 0.5 s |
| Adaptation, source objects, assets and renderer overhead | 0.7 s |
| Path/layout/ownership preparation | 0.8 s |
| Junction, separator and ordinary surface construction | 2.5 s |
| Both numerical and stored-precision validation | 1.0 s |
| Terrain preparation and fitting | 1.0 s |
| Remaining Blender publication | 0.5 s |
| **Whole Import target** | **7.0 s** |

This budget requires roughly an order-of-magnitude reduction in expensive
surface construction, a major reduction in repeated geometry checks and
terrain queries, and improvement in the unprofiled adaptation remainder.
It may prove unattainable for all general inputs. Report the fast path's
**share of execution time eliminated**, not only the percentage of junctions
using it: one expensive fallback can dominate the entire import.

## 4. Development stages and concrete tasks

### Stage A — reusable performance and correctness evidence

**Deliverable:** a repeatable optimization harness and a small replay set of
the most expensive geometric classes.

1. Retain the three-run baseline and separate numerical/publication profiles
   added for this task. Continue using the normal Import and the production
   scheduling preference.
2. Extend diagnostics with chart family, attempt success/failure and cause,
   boundary size, output faces, repair iterations, kernel calls, exact cache
   hits/evictions and peak workspace. Keep instrumentation optional and separate
   from acceptance timings.
3. Replay dependency closures for expensive junctions and separators from
   serialized numerical inputs. Use these for rapid development, then run the
   complete Import. IDs select test cases only; runtime dispatch uses geometry.
4. Also measure the remaining adaptation/source-object/material phase before
   attempting its 0.7 s budget. The present function profile covers preparation
   and publication, not that entire remainder.
5. Establish separate measurements for fresh Import, repeated Import in one
   process, no-op Update, vertex edit, width edit, and terrain edit. No cached
   Update result counts as a first-generation speedup.
6. Freeze reference scenes and repeatable camera/lighting/material setups for
   the visual-preservation comparisons in section 5. Establish surface and
   normal comparison tolerances before evaluating a changed construction;
   do not widen them afterward to accept a faster candidate.

`check_python_street_update_profile.py` now explicitly loads this checkout,
registers inside the isolated process and opens scenes with `use_scripts=False`.
It supports separate vertex/width measurements and a cold/warm terrain-cache
comparison. Keep this explicit bootstrap when extending the harness; do not
resolve a sibling main checkout through `from blosm` alone.

**Exit gate:** stable receipts, selected checkout verified, identical saved
baseline/profile results, and reproducible expensive local workloads.
The full-Import baseline and independent equality check are complete. The
implementation now also has an optional nested attempt ledger, captured kernel
and separator replays, and initial separate vertex/width/no-op Update measurements.
See the implementation record for qualification status and outstanding work.

### Stage B — remove repeated work while preserving construction behavior

Implement separately benchmarked patches, starting with the first two items.

1. **Pattern construction:** replace the repeated scan for the longest boundary
   edge in `PatternLibrary.quad` with an indexed priority queue. Use lazy
   invalidation for edges replaced by a split and deterministic tie ordering
   corresponding to the current owner insertion order. Avoid calling a general
   vector norm on each two-component candidate. Test exact pattern UVs,
   connectivity, orientation and boundary order across equal/unequal counts,
   parity transitions and cache eviction. This attacks the measured 1.46 million
   scoring calls. If preserving tie behavior is impossible, qualify the result
   explicitly as a construction change.
2. **Visibility kernels:** vectorize each half-plane clipping step, reuse edge
   directions/lengths, and cache exact repeated polygon queries within a build.
   Keep clipping order, epsilon semantics and the final clearance test. Store
   proven failure as well as success. Cache keys include exact coordinates,
   winding and tolerance. The splitter already caches index tuples within one
   split; extend reuse only where new measurements show repetition. Removing
   redundant constraints requires proof against all original half-planes;
   physical curved boundary samples must remain unchanged.
3. **Topology workspaces:** share immutable edge incidence, boundary masks,
   adjacency and face-color groups between quality and slope stages. Invalidate
   after edge rotations, refinement or changed connectivity. Cache by topology
   identity/version, not a mutable `Mesh` object address.
4. **Repeated path queries:** batch station requests and reuse exact frame
   results within footprint/boundary construction, keyed by path version,
   station and orientation. Keep the optimized cubic, continuous curvature
   certificate and source-deviation limits. Do not reduce interpolation safety
   sampling merely because path preparation is a small top-level task.
5. **Repair inner loops:** pre-gather constant color/corner arrays outside slope
   iterations; retain update order, line search, stopping limits and final tests.
   Quad energy scratch buffers and `bincount` gradients already exist.
6. **Boundary rules and payloads:** supply exact coefficients for constant
   plateau rules and known polynomial curb rules; retain generic fitting and
   residual checks for other rules. Share explicitly immutable UV buffers
   through Z-only placement, keeping defensive copies for external inputs.
   These are secondary to the measured kernel and pattern costs.

**Exit gate:** identical sources, topology, published positions, UVs, materials
and geometric construction results for changes that preserve operation order;
any numerical reordering gets a declared tolerance and adversarial comparison.
No geometry validation is removed. A **55–60 s** whole-Import checkpoint is an
initial target to evaluate, not a guaranteed sum of estimated savings.

### Stage C — construct suitable quad meshes directly, avoiding search and repair

**Deliverable:** general faster construction families with explicit eligibility
checks and the existing generator as fallback.

The code already has direct corner bands, Coons-style maps, a Bundle atlas and
cached quad patterns. The new work is extending their qualified coverage to
the costly cases currently reaching generic partition/search/repair.

1. Add a small strategy interface: `classify(inputs)`, `construct(inputs,
   workspace)`, `validate(result)`. A failed eligibility or quality test returns
   a structured fallback reason. Classification uses arm count, angles, width
   ratios, ownership, holes, rim ordering and monotonicity, never fixture IDs.
2. **Paired-rim separators first:** derive correspondence from source-side order
   and stations. Construct a central ribbon and small end caps; reconcile unequal
   boundary counts with direct quad transitions. Require non-crossing pairing,
   exact pavement-adjacent rim XYZ, correct plateau and valid convex cells.
   The four expensive partitioned examples and their perturbed variants are
   the first experiments. Do not force this branch on non-monotone or ambiguous
   regions; preserve their current fallback.
3. **Ordinary intersections:** use the existing semantic mouth and rounded
   corner boundaries to assemble adjacent half-profile sectors and a structured
   central patch. Support T/general three-arm, skew, four-arm and larger arm
   counts, unequal widths and divided approaches. Preserve boundary sample
   identities. Prove each mapping's positive orientation and qualify its quality
   before bypassing iterative repair.
4. **Compounds and holes:** reuse the physical union and source ownership, but
   first partition along known mouth, corner and island features. Use the same
   structured families for suitable resulting blocks. Generic visibility
   decomposition remains for unresolved blocks.
5. Precompute more than face connectivity: semantic boundary maps, transition
   partitions, incidence data and admissibility tests. Optional generated pattern
   packs need versioned keys, a reproducible generator and a bounded runtime
   cache. They must not bake a particular street's world coordinates or heights.
6. For Update, try the previous accepted family/atlas first, but map it onto the
   edited inputs and validate again. Never reuse stale geometry solely because
   the arm count stayed the same.

The objective is **avoiding** hundreds of repair steps on an unsuitable seed,
not accepting the same unsuitable seed after fewer checks. Record both time
saved and fallback time. Reordering successful families can change interior
topology; this stage uses geometric equivalence gates rather than claiming
bit-identical meshes. The visual-preservation gate in section 5 also applies:
retaining the boundary and passing slope limits alone does not establish that
the interior surface looks the same. Use the existing construction wherever
the faster family cannot preserve the reference appearance.

**Exit gate:** one new family gives a substantial measured improvement on an
expensive geometric class and passes generalization tests before adding more.
Evaluate a **25 s** whole-Import checkpoint after broader coverage, keeping
the complete result and existing fallback. This target remains experimental.

### Stage D — terrain query batching and validation throughput

This stage can be developed independently of the new chart families and is
necessary for the 10× target.

**Terrain:**

1. Deduplicate exact XY samples in bounded batches or component-local caches.
   Reuse heights in the final vertex lift check when coordinates are identical.
   Measure lookup/sorting overhead as well as rays avoided. Do not quantize XY
   or reduce the sample set for this first change.
2. Prototype a batch height-query provider for an evaluated triangle terrain:
   spatially bin candidate terrain triangles, evaluate XY barycentric queries
   in NumPy batches, and select the same highest valid hit below the ray origin.
   For boundaries, degeneracy, vertical faces, multiple layers or ambiguous
   tolerances, use the existing BVH ray as the reference/fallback. Accept only
   after comparison on realistic terrains and adversarial query positions.
   Build cost and memory count against the result; use BVH for small workloads
   when batching does not pay.
3. Separate terrain signature/evaluated-mesh preparation from BVH construction,
   creating the sampler lazily when fitting is actually required. For Update,
   invalidate by evaluated terrain, transform, controller frame and units.
4. Keep sampled no-cut coverage, complete profiles, shared junction minima,
   grade and maximum-lift checks. The current method is a sampled constraint;
   do not claim a continuous terrain-clearance proof.

**Validation:**

1. Measure the cost of broad phase, pair filtering and SAT independently. The
   implementation already has a spatial CSR grid, deduplicated pairs and
   bounded vector batches; merely proposing a spatial index is not new work.
2. Tune bounded batch sizes and reuse immutable face/index workspaces. Consider
   two-level part/patch and face bounds so unrelated regions avoid face work.
3. Give certified structured charts a reusable validation certificate covering
   their exact connectivity and coordinates. Internal pairs may be excluded
   only where injectivity/non-overlap is actually established; positive local
   cell Jacobians alone are not sufficient for arbitrary global geometry.
   General or altered charts keep the full checks.
4. Reuse conservative candidate-pair supersets across numerical and publication
   phases where valid. Still test the actual float32 stored geometry with its
   own tolerance. Rounding and large coordinate frames can create overlaps
   absent in float64. Do not simply delete the second overlap validation.

**Exit gate:** equality for exact sample reuse; independent BVH and all-pairs
overlap oracles for changed algorithms; unchanged failure behavior on missing
terrain, collapsed/inverted faces and real overlap. No unbounded whole-network
query cache. Terrain and both validation stages must approach their respective
1 s budgets to support the final target.

### Stage E — optional adaptive density and cheaper height construction

These are separately qualified construction changes, not silent reductions of
the user's selected resolution.

**Adaptive density:** 68% of faces lie in ordinary road parts. Use longitudinal
curvature, profile approximation error, terrain variation and marking/UV needs
to propose spacing, with compatible transitions at contacts. Long straight
regions can use fewer samples only when the profile, terrain and cell-quality
bounds permit it. Keep rounded boundaries at equally accurate or better
sampling, exact mouth partitions, gradual size changes and the current maximum
edge-ratio limit. Avoid long thin quads and visible size jumps.

First measure performance at today's sampling policy. Report a new adaptive
policy separately, with face counts and errors; do not claim a same-settings
speedup by moving the resolution slider or weakening tolerances. Define and
review the policy's error contract before making it a default. Reduced density
qualifies only if the visual-preservation comparisons pass, including close-up
silhouettes, surface shading, texture distortion and terrain contact. If they
fail, retain the original density for that region.

**Height construction:** only if it remains important after the chart work,
prototype a low-degree patch basis with exact incoming profile collars and
boundary gradients, zero outer curb relative to the common datum, and the
highest incoming crown as ceiling. Bernstein coefficient bounds can certify
elevation range; slope and continuity require additional derivative/Jacobian
checks. Positive blending weights alone do not prove gentle slopes, smooth
joins or drainage. Retain the current height/repair provider for ineligible
cases and preserve the sampled mesh checks.
Matching the crown ceiling and slope bounds is not sufficient: the new height
basis must also preserve the reference surface's visible curvature and shading.

**Exit gate:** documented geometric error bounds, original rounded outlines,
profiles, slopes and UV continuity; no regression in previously fixed ditch,
bump, skew-width or centerline-fit cases. Consider **12–15 s** as an intermediate
whole-Import target to test, not an established expectation.

### Stage F — bounded Update reuse and final integration

1. Cache terrain constraints and boundary topology for unchanged components.
   Refit the affected connected elevation-reference graph, including remote
   street/furniture outputs when grades or shared datums propagate. A lowering
   edit requires recomputation, not only an increasing-height propagation.
2. Keep the explicit Update workflow, source identity, unchanged output identity
   and transactional publication/rollback. Do not retain a second unlimited
   copy of every published mesh merely to make updates faster.
3. Resolve the remaining adaptation/material/source-object costs using their
   measurements. Material-only changes already have a mesh-preserving path.
4. Test the 7 s budget with all changes combined. If a small fallback class
   dominates, address that class next rather than reporting fast-case timings
   as the complete import.
5. Consider compiled numerical kernels only after the Python/NumPy architecture
   experiments locate a remaining unavoidable hotspot. This would need a
   separate Blender-extension packaging/maintenance decision. Do not assume
   SciPy, Numba, a compiler, or extra packages in Blender's installation.
   Multiprocess Blender imports are not the initial approach: they multiply
   memory, complicate deterministic ownership and cannot safely publish bpy
   data concurrently.

**Exit gate:** at least five ordinary runs with a median whole Import near
7.0 s (street generation at or below 7.2 s), complete coverage and all quality
gates. Publish the range, memory and fallback contributions. If the result
remains slower, report the actual improvement and remaining cost; the target
is not permission to drop output or checks.

## 5. Correctness and generalization contract

### Visual preservation is mandatory

Compare against the current accepted output for identical source inputs and
settings. Preserve the following properties:

- Street course, widths, outer silhouettes and rounded corner curves.
- Carriageway crown/crossfall, intersection curvature and smooth entrance joins;
  no new depressions, bumps, faceting or shading seams.
- Separator outlines, plateau and surface heights, and joins to the pavement.
- Terrain fit, clearances and visible contact with the original terrain.
- Material assignment and appearance; UV orientation, real-world texture scale,
  continuity and distortion, including diffuse and PBR asphalt materials.
- Conforming direct quads with the agreed reasonably uniform dimensions and
  gradual transitions; speed does not justify elongated or irregular cells.

The internal vertex count and connectivity need not be identical for a new
construction family. Its visible geometry and appearance still must be
preserved. Passing existing absolute quality limits alone is insufficient:
two different-looking surfaces can both lie below the same crown ceiling.

For every construction-changing optimization, require both kinds of evidence:

1. **Geometric comparison:** compare outlines and cross-sections, and compare
   surface heights, slopes, normals and texture coordinates at corresponding
   physical locations. Use shared spatial samples when vertex counts differ;
   include boundaries, entrances, interior extrema and terrain contacts. Keep
   analytical/continuous checks where already available. Record maximum as
   well as typical differences so a localized defect cannot be averaged away.
2. **Controlled visual comparison:** use identical cameras, lighting, materials,
   smoothing, normal maps, render settings, random seed and color management.
   Check tangent-space normal-map appearance as well as base color. Compare top, side and
   oblique views, whole fixtures and close-ups; use neutral smooth shading,
   grazing light, textured asphalt and a wireframe inspection. Include
   exaggeration 1 for normal appearance and 10 to expose surface defects.
   Retain side-by-side images and difference views. Image comparison complements
   the geometry checks and must not hide defects through downsampling or
   render noise.

Also compare small source edits followed by Update on both sides of a
construction-family boundary. Changing the selected algorithm must not introduce
new visible popping or abrupt surface redistribution.

Keep generated comparison scenes/renders under `artifacts/streets`; illustrations
included in documentation stay in the repository. These visual comparisons are
planned acceptance work, not evidence already produced by the profiling runs.
If a candidate changes the visible result, fix it or retain the existing
construction for that input. Report the smaller achieved speedup rather than
trading away appearance to meet the time budget.

### Numerical and generalization checks

Two categories of optimization require different acceptance tests:

| Category | Required evidence |
| --- | --- |
| Reuse / allocation / exact batching | Identical source data, ownership, face topology, float32 prepared/final coordinates, UVs, materials and accepted construction results; include cache invalidation and negative tests |
| New chart, height or sampling family | Exact authoritative boundary identities; unchanged geometric tolerances; preserved visible surface and texture mapping; full independent topology, quality, slopes, coverage, terrain and UV checks; controlled visual comparisons at exaggeration 1 and 10 |

Across both categories retain:

- Rounded corners, exact/smooth centerline behavior and the recent sparse-span
  alignment fix; no OSM-specific runtime exceptions.
- Complete entrance profiles, zero outer edges relative to the shared datum,
  highest-crown ceiling and existing longitudinal slope/depression guards.
- Direct quads, conforming seams, manifold topology, convex positive cells,
  existing edge-ratio limits and reasonably uniform transitions.
- Bundle overlap ownership, separator boundaries/plateau, retained holes and
  correct boundary-junction import behavior.
- Terrain sample coverage/clearance, grade and lift limits; correct material,
  UV and future furniture dependency behavior.
- Existing 1 GiB system reserve and bounded buffers. Keep the current Berlin
  workload near its approximately 0.6–0.7 GiB footprint; any material increase
  needs an explicit measured tradeoff, not a requirement for more user RAM.

Use all six intersection cases: T, general three-arm, non-right-angle, five or
more arms, unequal widths, divided carriageways/separators. Include four-arm
crosses, sharp dimensional transitions, short streets, near-collinear/acute
angles, holes, overlapping carriageways, translated/rotated/scaled synthetic
copies, small vertex edits and adversarial invalid inputs.

Karl-Marx-Allee is the principal benchmark, not the product definition. Reuse
the existing additional real OSM fixtures in the test corpus before acceptance.
The original function profiles in section 2 cover **Karl-Marx-Allee only**.
Subsequent real-network and manual-Update qualification is recorded in the
[implementation report](STREETS_PERFORMANCE_IMPLEMENTATION.md).

Relevant suites include `test_python_street_local_smoothing`, transitions and
mixed/smooth transitions, quad charts/four-side charts/rotations/workspaces,
sampling, overlap batches, street terrain and UV tests. Use native Blender
saved-scene checks in addition to pure numerical tests. Test vertex/width edits,
material-only Update, terrain changes, undo/reopen and failure rollback.

## 6. Recommended implementation order

1. Add per-attempt diagnostics and replay the expensive geometry classes.
2. Deliver the longest-edge priority queue and exact visibility-kernel
   optimizations independently, with before/after full Import measurements.
3. Start the **paired-rim separator fast-path experiment early**, alongside
   exact terrain-query reuse. Do not spend the whole effort on micro-optimizations
   if the required algorithmic reduction fails to materialize.
4. Expand successful structured construction to expensive junction/compound
   classes; optimize validation and batch terrain evaluation in parallel work
   streams, measuring Blender runs serially to avoid contention.
5. Decide from those results whether adaptive density or new height bases are
   necessary. These have larger geometric risks than exact reuse.
6. Integrate, test other networks and Update behavior, and assess the 10× target
   against complete unprofiled imports.

The first tangible delivery should include the unchanged-output optimizations,
their measured improvement, and a Blender comparison for the first new fast
construction family. Do not wait for the full tenfold target to deliver a
qualified improvement.

## 7. Code map and evidence

Paths below are relative to this Blosm checkout.

| Concern | Main source files under `way/street_network/python_backend/` |
| --- | --- |
| Pattern generation/cache | `patterns.py` |
| Visibility and compound decomposition | `quad_charts.py`, `decomposition.py`, `partition.py` |
| Separator/Bundle chart dispatch | `separator_surfaces.py`, `four_side_charts.py`, `bundles.py`, `bundle_compounds.py` |
| Junction ownership/direct bands | `ownership.py`, `intersections.py` |
| Repair workspaces/height fields | `quad_workspace.py`, `quality_correction.py`, `slope_projection.py`, `surface_fields.py`, `conditions.py` |
| Path/frame queries | `paths.py`, `routes.py` |
| Terrain | `terrain/alignment.py`, `terrain/corridor.py`, `terrain/blender.py` |
| Numerical/stored-precision validation | `validation.py`, `backend.py`, `mesh.py`, `scene.py` |
| Update dependencies and publication | `updates.py`, `scene.py`, `materials.py` |

Generated evidence root:
`C:/dev/workspace/artifacts/streets/performance/optimization_plan_20261111/`.

- `baseline/summary.json`, `baseline/run_01..03/timings.json`: three unprofiled
  receipts, per-owner task durations and bounded-process resource records.
- `baseline/run_03/street_benchmark.blend`: complete current output.
- `profile/run_01/numerical.prof`, `publication.prof` and their `_profile.txt`
  companions: separate Python call profiles; `timings.json` labels overhead.
- `profile/run_01/street_benchmark.blend`: profiled output.
- `equivalence/report.json`: independent saved-output comparison.
- `workload/report.json`: mesh-kind counts and accepted chart metadata.

The baseline launcher began before the profiler-scope option was added to the
benchmark tools. Its final summary initially rejected the mixed old/new receipt
schema, although all three child runs completed successfully. The summarizer
now normalizes a missing unprofiled scope to `null`; rerunning summary succeeded
without changing raw timings or rerunning generation. Runtime files stayed at
the same revision. This bookkeeping issue is not a slow or failed import.

### Reproduction

Run from this chat's development checkout. Use fresh output folders for every
new batch. The host tool configures the isolated profile and resource guard;
do not launch these tests inside the user's running Blender.

```powershell
& 'C:/Users/vladi/AppData/Local/Programs/Python/Python310/python.exe' `
  tools/benchmark_streets.py `
  --osm C:/dev/workspace/osm_extracts/facade_visibility/berlin_karl_marx_allee.osm `
  --scene C:/prochitecture/setup/bldg2_berlin_karl_marx_allee.blend `
  --operator --import-performance --runs 3 --threads 16 --numeric-threads 1 `
  --profile-root C:/dev/workspace/artifacts/blender_profiles/01a0807c-abda-7822-8855-02127742633c `
  --output C:/dev/workspace/artifacts/streets/performance/next_baseline
```

For a separate profile, use a fresh output directory and replace `--runs 3`
with `--runs 1 --diagnostic-profile --profile-scope all`. The default profiler
scope remains `numerical`; `publication` is also supported. Do not mix profiled
and ordinary reports in one timing summary. `tests/inspect_street_workload.py`
reads a saved scene without regenerating geometry; the independent equality
checker is `tests/check_street_timing_equivalence.py`.
