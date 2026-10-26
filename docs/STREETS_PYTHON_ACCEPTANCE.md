# Python streets: deliverable and verification record

Reviewed 25 October 2026. Target: **Blender 5.2**, branch **`streets_python`**,
checkout `C:/dev/workspace/streets_dev/blosm`. The buildings checkout at
`C:/dev/workspace/blosm` remains on `buildings_gn`.

**Current development:** dev72, fingerprint
`1066f7fbcd448e7e77e358acc99ab2d7cd3a8caaa7b0653bfe7c8610856fca5d`.
The bounded matching-crown mixed-join fallback keeps source XY/curbs/direct
quads, whole mouths plus neighbor rows and the entire crest fixed. It projects
only interior heights, retains the ordinary physical/entrance slope bounds,
and requires transverse downhill paths and a basin-free boundary-spill check.
The [complete final-source regression](../../experiments/streets/python_v1/suite_dev72_audit_01.json)
now passes **982 tests / 114 modules** in one uninterrupted run, with no errors,
failures or skips. This includes the formerly failed short/sharp cases at all
three resolutions. The preceding [93 focused checks](../../experiments/streets/python_v1/mixed_height_focused_dev72_01/focused_checks.json)
retain their earlier pre-catalogue snapshot scope.

The [112-fixture dev72 gallery](../../experiments/streets/python_v1/gallery_generated_dev72_01/generation_checks.json)
passes authored generation with 726,698 direct quads, no GN modifiers and
nonoverlapping cells. [All 110 predecessor fixture sets remain exact](../../experiments/streets/python_v1/gallery_capture_dev72_01/catalogue_checks.json).
Both [saved-buffer quality](../../experiments/streets/python_v1/gallery_quality_dev72_01/import_quality_checks.json)
and [independent complete-gradient checks](../../experiments/streets/python_v1/slopes_gallery_dev72_01/slope_checks.json)
pass, with explicit float32 uncertainty in the full-gallery gradient proof.
Both [new native editing cases](../../experiments/streets/python_v1/mixed_height_native_dev72_01/native_checks.json)
pass cold saved identity, real Edit Mode/width/Update, mode changes, local/full
mesh and furniture, remote identity, restoration and R2/E0/R4/E10. The two
actual saved joins also pass a [strict measured incoming-profile slope check](../../experiments/streets/python_v1/mixed_height_published_dev72_02/published_checks.json)
without an uncertainty allowance, exact zero curbs and drainage. Fixed crest
rounding is explicitly bounded by its native half ULP; the earlier verifier
requiring bitwise constant exaggerated crest floats is retained as failed.
The current isolated installation and delivery checks are recorded below.

All [112 actual cold Updates](../../experiments/streets/python_v1/gallery_cold_dev72_01/cold_update_checks.json)
now pass with exact saved source/mesh/furniture/layout and retained object/data
identities. The [seven-job current gallery audit](../../experiments/streets/python_v1/gallery_dev72_audit_01.json)
reconciles terminal results, fixture inventories, source/scene/script hashes,
caps, one numerical thread and the 1 GiB reserve. Current full regression is
independently audited against every ordered test event, module inventory,
source/script hashes and terminal resources. It takes 1,072.52 seconds and peaks
at 189,472,768 sampled private bytes under 0.25 GiB. No installed-copy
qualification is inherited from the older release.

The [current real-input audit](../../experiments/streets/python_v1/dev72_real_receipt_audit_01.json)
also passes all nine serial normal-import/quality/reopen-edit jobs for Berlin,
Canal Street and Nashville. It checks complete coverage, prior prepared
identities/parameters and exact reference mesh buffers, native seam/width edits,
incremental/full equality, restoration and unrelated output retention. Local
Update measurements are 1.66 / 7.63 / 0.94 seconds, with at most 715,980,800
sampled private bytes under 0.75 GiB. These are three individual measurements.

All [31 current native asset tests](../../experiments/streets/python_v1/native_assets_dev72_01/native_asset_checks.json)
and [five current interaction jobs](../../experiments/streets/python_v1/native_closure_dev72_checks_01.json)
pass. They cover Base/Pro/external datablocks and exact-revision recovery,
saved editing bootstrap, late-publication rollback, Undo/Redo, source grouping,
cache release and absence of a continuously running geometry observer. The
installable candidate now passes the isolated-copy queue below.

**Qualified installation candidate:** [dev72 ZIP](../../experiments/streets/python_v1/packages/dev72_01/blosm_streets_python_dev72.zip),
SHA256 `25a7b4891bc119b285191e5eafa67bab5e97f4fc4d4c61b344d4a710b6ca0aa6`.
The [installed-copy audit](../../experiments/streets/python_v1/package_dev72_qualification_01.json)
reconciles 13 terminal native jobs and 43 installed fixture checks, every
archive/installed payload and all 1,245 Python files. Full addon registration,
backend selection, external material loading, manual Update, normal import,
quality and reopening/edit/restoration pass with installed-only runtime modules.
Original preferences were not saved. No public publication or Git commit was
performed.

The [actual registered Import operator](../../experiments/streets/python_v1/installed_ui_import_dev72_01/operator_import_checks.json)
also passes from that installation with normal menu settings and bundled street
assets, without setup/renderer substitution. Its [saved quality check](../../experiments/streets/python_v1/installed_ui_quality_dev72_01/import_quality_checks.json)
and [reopening/edit/restore check](../../experiments/streets/python_v1/installed_ui_reopen_dev72_01/import_reopen_checks.json)
pass separately. UI float properties produce a distinct 21,838-quad input run;
the 21,818-quad direct-pipeline run retains its exact-reference qualification.

| Initial plan stage | Current evidence |
| --- | --- |
| Shared contracts/classes | Prepared identities and parameters retained in all three real inputs; grouping/mixed/regrouping checks; injectable providers and shared asset service |
| Reusable patterns | Immutable direct-quad grids/strips/atlases, bounded LRU and per-part limits; complete 982-test suite and native published-mesh checks |
| Connected surfaces | All six named families in the 112-fixture gallery: T, general three-arm, oblique, five-plus, unequal widths and divided carriageways/separators; rounded/profile/height/slope/quad checks |
| Manual Update | Actual button/Import operators, source/width/mode edits, old/new interactions, local/full equality, rollback, Undo/Redo and unchanged remote output identities |
| Furniture anchors | Source-relative station/side/offset, movement and width changes, actual-surface queries, buffer release, saved-file reconstruction and invalid-anchor diagnostics |
| Integration/qualification | Normal sample-independent Blosm pipeline, three real inputs, 112 authored fixtures/cold Updates, 31 native asset tests, lifecycle/grouping checks and the qualified installable ZIP |

This review is for the original **initial** six-stage Python milestone. Stage 6
explicitly requires unsupported construction to be reported and postpones full
terrain and later transport features. The limits listed later in this document
remain unresolved; no universal all-configuration solver, completed GN network,
managed topology commands or production furniture assets are claimed. The
current guide is consolidated in [the quickstart](STREETS_PYTHON_QUICKSTART.md).
The [delivery review](../../experiments/streets/python_v1/initial_delivery_dev72_audit_01.json)
maps those six stages to actual evidence, verifies all 39 terminal native jobs,
the exact source/archive/scene identities, the three UI-path checks and current
guide links. It also records the unresolved broader limits rather than claiming
an all-configuration solver. The input snapshot confirms that the UI/direct
comparison differs only in the two components' float32 cell-size settings.

**Previous checkpoint:** dev71, fingerprint
`0aa98c264e079a3c8829b7ecb8fcabb98e050c89cd9a6cc4d87b22e949a823bc`.
Equal-profile Exact/Smooth bent seams now use the existing mixed geometry when
their endpoint frames disagree. Task declaration and evaluated layout make the
same decision; matching dimensions no longer bypass the necessary join.
[Eight current-source tests](../../experiments/streets/python_v1/uniform_focused_dev71_03/focused_checks.json)
pass 28 signed/order/length cases, mouth profiles, source loci, constant crown,
mode edits/removal/restoration, source/width edits, local/full/furniture,
neighboring-junction clearance and separate/overlapping Bundle ownership.
The [combined focused audit](../../experiments/streets/python_v1/uniform_focus_dev71_audit_01.json)
reconciles 76 passing IDs. It proves that the earlier production source is
byte-identical; only the new optional-mesh assertion changed. The two original
Bundle verifier errors are retained. This is not a single green 76-test run
or a full current-source regression claim.

The [new authored 110-fixture gallery](../../experiments/streets/python_v1/gallery_generated_dev71_01/generation_checks.json)
is terminal PASS with 718,954 direct quads and no overlapping gallery cells or
GN modifiers. Four new examples include both standalone orders and separated
and merged Bundle pavements. All
[four native cases](../../experiments/streets/python_v1/uniform_native_dev71_01/native_checks.json)
pass saved-file cold identity, actual Edit Mode vertex/width/Update, mode
changes that remove and restore owned joins, local/full mesh and furniture,
remote identity, exact restoration and R2/E0/R4/E10. Peak sampled private usage
is 356,368,384 bytes under 0.5 GiB, one thread and the fixed 1 GiB reserve.
The [independent six-job gallery audit](../../experiments/streets/python_v1/gallery_dev71_audit_01.json)
verifies actual authored generation, all 106 predecessor source/mesh/furniture/
layout sets exact, actual float32 quality, independent complete bilinear
gradient bounds with explicit float32 uncertainty, all 110 cold Updates and
the four native edit sequences. The fixture inventory, scene/source hashes,
scripts, terminal resources and caps are reconciled. Current full regression,
real inputs and installed-copy qualification remain separate gates.

**Metric-station checkpoint:** dev70, fingerprint
`19707c80d571f298841cb2e27518ccf0715c68442213c7a7fd71501f4802c94c`.
The source-relative furniture conversion now uses the source path's metric
station inverse and the original native parameter for its mixed-transition
piece. An independent probe measured up to 53.3 mm XY error before the fix;
the [same probe now reports zero](../../experiments/streets/python_v1/mixed_anchor_station_dev70_01/anchor_station_checks.json)
for its six Smooth-leg interior stations. A new permanent regression reproduces
the original failure. [64 focused tests](../../experiments/streets/python_v1/anchor_focused_dev70_01/focused_checks.json)
pass both modes, source orientations and sides, endpoint/interior stations,
local/full equality and released-buffer reconstruction. These focused checks alone do not
qualify a complete dev70 gallery, real-network import or installed archive.
The [native three-fixture check](../../experiments/streets/python_v1/anchor_native_dev70_03/native_checks.json)
also passes saved-file upgrades, real Edit Mode/width edits, actual Update,
local/full equality, restoration, R2/E0 and R4/E10. Upgrade mesh buffers, source
data and semantic keys remain exact; only the Smooth-leg furniture anchor is
corrected. Current-version cold Updates retain object/data identity, and
ordinary edits retain the remote road. Native anchor XY error is at most
0.000000111 m. The worker peaks at 346,828,800 sampled private bytes under
0.5 GiB, one thread and the same 1 GiB reserve. Failed verifier attempts are
retained: attempt 01 incorrectly demanded native identity across a VERSION
change, and attempt 02 treated an offset profile centre as its source
centerline. Attempt 03 retains the actual curve check with the documented
offset compensation; no physical mesh guard is relaxed.

The [freshly generated 106-fixture dev70 scene](../../experiments/streets/python_v1/gallery_generated_dev70_01/generation_checks.json)
is terminal PASS and contains 699,522 direct quads. The
[independent saved capture](../../experiments/streets/python_v1/gallery_capture_dev70_01/catalogue_checks.json)
verifies all 106 source/mesh/layout sets against dev69 and all 103 older
source/mesh/furniture sets against the qualified checkpoint. Only
`mixed_dimensional_3_0 / placement:85` is corrected. The
[saved mesh comparison](../../experiments/streets/python_v1/gallery_quality_dev70_01/import_quality_checks.json)
passes convexity, ratio 8, crown envelope, overlaps and exact reference buffers.
The [independent gradient check](../../experiments/streets/python_v1/slopes_gallery_dev70_01/slope_checks.json)
passes all 699,522 quads with explicit float32 uncertainty in its complete
bilinear bound. All [106 cold reconstructions](../../experiments/streets/python_v1/gallery_cold_dev70_01/cold_update_checks.json)
now pass using the actual Update operator, preserving saved source, mesh,
furniture, layout and object/data identities. The
[complete gallery audit](../../experiments/streets/python_v1/gallery_dev70_audit_01.json)
reconciles all five terminal jobs, exact scene/source and resource limits.
Current full regression, real-input and installed-copy qualification remain
separate gates; this is not an install-qualified dev70 release.

**Preceding dimensional-seam checkpoint:** dev69, fingerprint
`ae0d7fef71288fbafc64b992b1296b30e255e14ecec337c886effbc090616b16`.
Bent dimensional seams with one Exact and one Smooth source now use bounded
piecewise polynomial curbs, keeping both source centerline curves and all
ordinary surface guards. The
[37 focused tests](../../experiments/streets/python_v1/mixed_transitions_focused_dev69_04/focused_checks.json)
pass. The [regenerated 106-fixture scene](../../experiments/streets/python_v1/gallery_generated_dev69_02/generation_checks.json)
and [saved quality comparison](../../experiments/streets/python_v1/gallery_quality_dev69_02/import_quality_checks.json)
pass all 699,522 quads, with all 634 previous mesh parts exact and three named
new networks. The [new native edit checks](../../experiments/streets/python_v1/mixed_native_dev69_01/native_checks.json)
pass actual Edit Mode/width/Update/furniture/local-full/restoration and R2/E0,
R4/E10 for all three examples. The
[complete regression audit](../../experiments/streets/python_v1/suite_dev69_audit_01.json)
verifies all 968 unique tests in 112 modules, in one uninterrupted green run
with zero failures, errors or skips and exact source before/after. It uses one
numerical thread, a 0.25 GiB cap and the unchanged 1 GiB reserve. This evidence
belongs to dev69 rather than automatically qualifying the later anchor fix;
current real-input/full-gallery cold/installed qualification remains separate.
This must not be described as an install-qualified dev69 or dev70 release. The 106-input catalogue keeps
all previous fixtures and adds three mixed seams, including unequal rises.
Short sharp, same-profile bent and route-spanning mixed cases remain open.

**Qualified preview checkpoint:** dev68, Python fingerprint
`376c108b156007d8688652a21ee28ff058fe760688e71bd8082210b961a176c4`,
adds a bounded planar corner-band seed after an ordinary single-junction
grading failure. The established compound solver still owns physical heights;
source controls, actual rim partitions, full entrance profiles and radii stay
fixed. Its [45 focused tests](../../experiments/streets/python_v1/corner_seed_focused_dev68_02/focused_checks.json)
pass physical/E0/E10, connected complete interfaces, whole-edge entrance bounds,
drainage and local/full/restoration checks, together with existing core,
junction, condition and surface regression. The extra three tests verify that
stored coordinates and bounds share float32 precision across rebuilt/retained
meshes and scene units. Publication tolerance is unchanged.

The [22-case matrix](../../experiments/streets/python_v1/normal_merged_matrix_dev68_01/matrix_checks.json)
passes all variants on the preceding a9949428 source, before that publication
fix. The [actual Blender check](../../experiments/streets/python_v1/corner_seed_native_dev68_04/native_corner_seed_checks.json)
passes the new case's width edit, moving furniture, local/full equality and
restoration on the current 376c108b source. Its
[103-fixture scene](../../experiments/streets/python_v1/corner_seed_native_dev68_04/python_streets_corner_seed.blend)
uses the freshly regenerated current-source 102-fixture gallery as its parent,
then preserves those source/mesh/furniture/object sets while adding the new case.
The [independent saved quality check](../../experiments/streets/python_v1/corner_seed_quality_dev68_02/import_quality_checks.json)
passes all 688,258 published quads: convexity, ratio 8, no overlap, crown envelope
and exact previous mesh buffers. The [fresh reopening](../../experiments/streets/python_v1/corner_seed_reopen_dev68_02/reopen_checks.json)
repeats the new case's cold identity/edit/full/restoration checks and retains
the scene hash. Maximum sampled resident usage across these native jobs is
517,672,960 bytes under 0.75 GiB, one numerical thread and the same 1 GiB reserve.

The [current regression audit](../../experiments/streets/python_v1/suite_dev68_reconciled_audit_01.json)
reconciles all **960 unique tests in 111 modules** with passing evidence on the
same source. The original main run has 950 passes and ten asset-test failures
caused by sandbox temporary-directory rename permissions. All 21 tests in the
affected asset modules pass with the standard workspace temporary directory.
Original failures are retained: this is combined coverage, not one uninterrupted
green run. No street geometry failure was waived.

The [regeneration audit](../../experiments/streets/python_v1/gallery_dev68_regeneration_audit_01.json)
and [saved-gallery audit](../../experiments/streets/python_v1/gallery_dev68_saved_audit_01.json)
verify all **102 authored fixtures** generated from catalogue inputs, exact
reference source/mesh/furniture, saved float32 quality and all 102 native cold
Updates with retained object and mesh identity. The parent gallery has 677,596
quads across the complete scene; this is not a per-Bundle count.

The [real-input audit](../../experiments/streets/python_v1/dev68_real_receipt_audit_01.json)
verifies all nine fresh import/quality/reopen jobs for Berlin, Canal Street and
Nashville. Prepared identities and saved meshes remain exact; source seam/width
edits pass incremental/full equality, restoration and unrelated identity.
Observed local Updates take 1.64 / 7.51 / 0.92 seconds respectively. These are
individual observations, not controlled speedup measurements. Maximum sampled
private usage is 712,732,672 bytes under a 0.75 GiB cap, with one numerical thread
and the unchanged 1 GiB reserve.

The [native interaction closure](../../experiments/streets/python_v1/native_closure_dev68_checks_01.json)
passes five jobs covering bootstrap, failed-publication rollback, Undo/Redo,
cache lifecycle and grouped/per-Street/mixed inputs. Undo/publication checks
load their target fixture selectively rather than duplicating the whole gallery
into history. All assertions in the original tests are retained.
The [native asset check](../../experiments/streets/python_v1/native_assets_dev68_01/native_asset_checks.json)
passes all 31 tests for shared libraries, asset kinds, external identity,
Base/Pro policy and saved-library recovery.

The [new addon ZIP](../../experiments/streets/python_v1/packages/dev68_01/package_checks.json)
has verified payloads and current Python sources. The
[independent installed-copy audit](../../experiments/streets/python_v1/package_dev68_qualification_01.json)
now verifies twelve terminal passing native jobs: actual installation/full
registration, both backends, external materials, manual Update, 34 selected
fixtures, normal Canal import, saved quality and fresh seam/width editing.
All runtime modules originate in the installed archive; its payloads remain
exact. Peak sampled private memory is 701,452,288 bytes under 0.75 GiB.
The initial combined check hit that cap and remains recorded as incomplete.
Separate fixture-family processes preserve all assertions and pass without
raising the limit. No user preferences are saved; the archive is not published.
Its immutable build manifest records pre-installation status; the linked audit
is the subsequent installation evidence. The complete records below are
explicitly **dev67** evidence.

The [current independent saved-gradient audit](../../experiments/streets/python_v1/slopes_gallery_dev68_01/slope_checks.json)
passes all 103 networks / 634 mesh parts / 688,258 quads. It checks four corners
and the complete bilinear bound with explicit float32 uncertainty. The test
compares uncertainty-adjusted gradients with the stored policy and physical
maximum; it does not prove a strict measured 10% bound without that allowance.

**Preceding implementation:** dev67 integrates the previously external merged
entrance/run construction into normal preparation, physical ownership and
Bundle output tasks. Explicit source-mouth aliases retain both endpoints and
stations, and a physical rim view supports E0 without dividing by exaggeration.
The [normal-backend proof](../../experiments/streets/python_v1/merged_backend_dev67_01/preparation_checks.json)
passes both modes, complete topology/interfaces, Update-equivalent edits,
local/full/furniture equality, restoration, unrelated identity and E0. The
[native proof](../../experiments/streets/python_v1/merged_native_dev67_01/native_merged_entrance_checks.json)
passes the actual Update operator, combined control/width edits, restoration
and cold reconstruction. Its scene has 102 editable fixtures; all original
100 source/mesh/furniture sets remain exact. These receipts belong to source
`58d5ea01db56b18d3c3215974b30e6781576e36306301c38a7f1f38a5f16f2d8`.
The native worker peaks at 498,216,960 sampled private bytes under 0.75 GiB,
one numerical thread and the unchanged 1 GiB reserve.

**Final source for this checkpoint:**
`30fac1a321c1defbe838d5889083cc7d31b1fa66ccef59005c02b8ebc0adfdd6`.
It includes disjointness, port-order, catalogue and consumed-centreline support
guards. The [normal-dispatch matrix](../../experiments/streets/python_v1/normal_merged_matrix_dev67_01/matrix_checks.json)
passes 21/22 with both baseline restorations exact; the retained failure is the
existing Exact translated-datum outer junction 210 quad grading. R2/R4, reversed
storage, unequal widths/crowns and T/skew/five-arm cases pass. E1 Bundle counts
are 3,625 / 4,738; whole-network counts are 7,720 / 8,801 including the separate
1,600-quad control road. The terminal failed worker peaks at 147,873,792 sampled
private bytes under 0.5 GiB, with the same reserve and one numerical thread.
The [fresh native reopening](../../experiments/streets/python_v1/merged_reopen_dev67_01/native_merged_entrance_checks.json)
passes on this final source for both new cases, with actual Update, edits,
local/full/furniture equality, restoration and cold identity. The 102-fixture
file is unchanged. It does not claim regeneration of the previous 100 cases.
The [37-test regression](../../experiments/streets/python_v1/merged_regression_dev67_01/preparation_checks.json)
belongs to the immediately preceding 8abf8a source, before the two-line consumed
support guard; its modules and source scope remain explicit in that receipt.
The [complete regression coverage audit](../../experiments/streets/python_v1/suite_dev67_resumed_audit_01.json)
verifies all 953 ordered test IDs in 109 modules on this final source: 545
observed passing statuses before a system-reserve interruption, plus 408
terminal passing suffix events. No evidence IDs are missing or duplicated;
failures, errors and skips are zero. It is not an uninterrupted single-process
run. Private peaks remain at or below 190,189,568 bytes, with one numerical
thread and the unchanged 1 GiB reserve. The suffix cap is 0.25 GiB.
The [authored-regeneration audit](../../experiments/streets/python_v1/gallery_dev67_regeneration_audit_01.json)
now verifies all 102 cases generated from original catalogue inputs using the
Update operator, without importing saved output meshes. Original source,
positions, quad connectivity and furniture remain exact against the independent
capture. The four-column gallery has 677,596 quads in total. Peak resident/private
usage is 560,672,768 / 636,633,088 bytes under 0.75 GiB, one numerical thread and
the unchanged 1 GiB reserve.

The [saved-file mesh audit](../../experiments/streets/python_v1/gallery_quality_dev67_01/import_quality_checks.json)
also passes all 102 networks: actual float32 buffers contain positive convex
quads, maximum edge ratio 7.99635, no cells above ratio 8, no overlaps and no
heights outside the global incoming-crown envelope. All reference positions
and topology are exact. This does not independently establish all physical
gradients. Broader gallery editing/reopening and installer qualification
remain open. The new [editable gallery](../../experiments/streets/python_v1/gallery_generated_dev67_01/python_streets_gallery.blend)
has SHA256 `76d46f773e677ab5495e580f8f4c01aa24c143d06d3be296a636e0cd1b02e32f`.

All [nine current real-input checks](../../experiments/streets/python_v1/dev67_real_import_checks_03.json)
are now terminal success. The [independent audit](../../experiments/streets/python_v1/dev67_real_receipt_audit_01.json)
reconciles script/input/source/resource hashes and preserves all prepared
Street/Section/junction/Bundle data against dev53, excluding fresh network UUIDs.
Berlin, Canal Street and Nashville contain 398,972 / 21,818 / 263,213 quads.
Each has two READY components, COMPLETE coverage and no GN modifiers. Saved
float32 quality and exact prior mesh comparisons pass, as do cold reconstruction,
seam/width edits, incremental/full equality, restoration and unrelated identity.
Observed local Updates take 1.51 / 7.44 / 0.90 s. Maximum sampled resident/private
usage is 655,794,176 / 707,268,608 bytes under 0.75 GiB, one numerical thread and
the unchanged 1 GiB reserve. These individual runs establish no controlled
performance improvement or universal network guarantee.
The [source/receipt audit](../../experiments/streets/python_v1/merged_integration_dev67_audit_01.json)
and [verified recovery](../../experiments/streets/python_v1/checkpoints/dev67_01/checkpoint_checks.json)
preserve this scope; neither is an installer or a completed-goal claim.

**Preserved gradient checkpoint:** `0.1.0.dev66`, fingerprint
`91b1fa8e0dca0de47318a4f94403898256ed6360237cde82e52a20fbd8868a39`.
The production physical-gradient API and angular crown derivative pass
[22 focused tests](../../experiments/streets/python_v1/profile_gradients_dev66_02/gradient_checks.json),
including eight new differential/API tests and fourteen retained regression
tests. The old ordinary-profile arithmetic and checked six-family buffers
remain exact. The separate
[22 combined-entrance surface cases](../../experiments/streets/python_v1/combined_dynamic_matrix_dev66_02/dynamic_matrix_checks.json)
consume those actual gradients through the maintained compound kernel; all
boundary-gradient discrepancies are zero. Every patch retains its 2 m corners,
zero outer curbs, complete profiles, crown envelope, direct positive quads and
publication slope/ratio guards. Both baseline restorations remain exact.

The [22 connector cases](../../experiments/streets/python_v1/combined_bridges_dev66_01/bridge_checks.json)
weld two full-width dynamic-profile strips to each patch with exact complete
mouth edges. Joined counts are 369–3,293 quads, including higher resolutions.
The connector's maximum checked angular longitudinal/entrance bounds are
0.348%/0.048%. All three workers are terminal success under 0.25 GiB, one
numerical thread and the unchanged 1 GiB reserve. Their peak sampled private
allocations are 75,669,504 / 98,410,496 / 65,961,984 bytes respectively.

These resolve the missing gradient API and establish short merged-profile
connectors. A subsequent
[complete-network matrix](../../experiments/streets/python_v1/merged_network_matrix_final_dev66_01/matrix_checks.json)
joins the full runs, dynamic separator islands and original rounded outer caps:
21/22 variants pass, with exact baseline restorations. The remaining Exact
translated-and-float32-rounded case exhausts the existing outer junction's
quad-grading budget. It is preserved as a construction defect. Baseline Bundle
counts are 3,625 / 4,738 quads; whole-network counts are 7,782 / 8,777 including
an independent 1,600-quad control road. Positive direct quads, whole interfaces,
topology, non-overlap, crown and slope guards pass for those baselines. The
terminal matrix has exit code 1 for the retained failure and peaks at
159,428,608 sampled private bytes under 0.5 GiB with the same reserve/threads.

The [native Blender review](../../experiments/streets/python_v1/merged_network_review_dev66_02/review_checks.json)
preserves all 100 editable fixtures exactly and appends two frozen geometry
proofs. It peaks at 482,508,800 sampled private bytes under 0.75 GiB. Those
proofs do not support addon Update; their source references are locked.
Logical-owner and source-end mappings, normal addon task dispatch and
Update/furniture/reopening for the new branch remain unqualified.
There is no full dev66 regression,
real-network or installer qualification yet. Older receipts below keep their
original source scope.

The [native dev66 regression](../../experiments/streets/python_v1/profile_gradient_native_dev66_01/generation_checks.json)
passes nine existing fixtures through normal Update: six core families, a
compound divided crossing and the Exact/Smooth internal-overlap examples.
Sources, published mesh buffers and furniture match the independent earlier
capture exactly. The terminal Blender worker peaks at 271,523,840 resident /
336,453,632 sampled private bytes under 0.75 GiB, with a 1 GiB reserve and one
numerical thread. This is selective native regression of existing addon
branches; it is not a full-gallery or new merged-entrance Update proof.

The [independent saved-buffer publication check](../../experiments/streets/python_v1/combined_joined_publication_dev66_01/publication_checks.json)
passes all 22 joined meshes, including float32 full bilinear gradients, positive
quad/ratio guards, global non-overlap and merged connector curb agreement.
The maximum measured Smooth connector error is 3.573 mm, within 5 mm. The
worker peaks at 46,845,952 resident / 67,334,144 sampled private bytes under
0.25 GiB with the same reserve/thread policy. This is an audit of saved meshes,
not normal addon dispatch. The
[source/scope audit](../../experiments/streets/python_v1/profile_gradients_dev66_audit_01.json)
also verifies the preceding four terminal jobs, the four changed/two added
Python files, separate branches and empty index. Later document/figure edits
do not change its source or numerical proof. The
[dev66 recovery archive](../../experiments/streets/python_v1/checkpoints/dev66_01/checkpoint_checks.json)
is for continuation, not installation or completed-goal qualification.

**Qualified dev65 geometry and regression:** source
`f8c603afe359f877d40d102de6d597b64a909c0b995a93ad1834061e7ca3b3fb`.
The maintained-gallery tooling checkpoint source is
`d4535b48695b0f4d594c2720229e7fc80175fe50bf91e47ca6769eca2a3a9546`.
It adds only `tests/_python_street_catalogue.py` and
`tests/generate_python_street_gallery.py`; all 1,222 previous Python files,
including runtime and regression definitions, remain exact. The
[promotion record](../../experiments/streets/python_v1/gallery_tooling_promotion_dev65_01.json)
records that distinction. The maintained generator and saved-panel reopening
checks pass separately; the exact-scope audit below verifies the unchanged
runtime rather than claiming a repeated numerical run on the added tools.
The internal angular pavement-union provider is now integrated through explicit
domain, profile-sampling and cap policies. The
[25 focused runtime tests](../../experiments/streets/python_v1/internal_unions_dev65_03/remaining_suite_checks.json)
pass both interpolation modes, original-ribbon union/area agreement, full
published bilinear slope bounds, exact separator interfaces, four grass islands,
connected filled topology and a replacement below 8,000 quads. Combined edits,
exact local/full/restoration, furniture, unrelated identity, reversed storage,
translated-datum E0 and existing internal/dependency cases pass as well. An
eligibility guard retains the old provider for free-ended runs.

The [100-fixture preview](../../experiments/streets/python_v1/internal_unions_native_dev65_01/python_streets_internal_unions.blend)
retains the previous 98 saved source/mesh/furniture sets exactly and appends the
two overlapping internal cases without spatial overlap. Their
[actual Update checks](../../experiments/streets/python_v1/internal_unions_native_dev65_01/native_internal_union_checks.json)
and [fresh selective reopening](../../experiments/streets/python_v1/internal_unions_reopen_dev65_01/native_internal_union_checks.json)
pass edits, local/full/furniture equality, restoration, remote identity and cold
reconstruction. The older separated internal fixtures also reconstruct exactly
on this source. Native private usage peaks at 505,831,424 bytes under a 0.75 GiB
cap; reopening uses 0.5 GiB. The reserve remains 1 GiB.

The [saved-buffer audit](../../experiments/streets/python_v1/internal_unions_published_dev65_01/published_mesh_checks.json)
passes the original-ribbon footprint/area reference, full physical slopes below
9.059%, crown bounds, four grass disks and one connected filled network without
cracks. All 504/574 grass rim edges match complete float32 paving edges. The
[longitudinal bound](../../experiments/streets/python_v1/internal_unions_longitudinal_dev65_01/longitudinal_checks.json)
is below 2.882%, retaining the 3% limit. The Bundle replacements contain
4,450/5,475 quads, including 850/981 grass quads. Full current-gallery,
real-network, installer, higher-resolution and broader regression qualification
remain incomplete. The complete 933-test result below belongs to
the preserved dev64 source.

The [dev65 source/receipt audit](../../experiments/streets/python_v1/dev65_internal_audit_01.json)
verifies all six scoped successes, native/reopened buffer equality, scripts,
resource guards, document links and the empty Git index. Relative to dev64,
1,212 Python files remain exact, six changed and four were added. Existing
angular body/curve certificates, profiles and intersection builders are exact.

The [dev65 recovery verification](../../experiments/streets/python_v1/checkpoints/dev65_01/checkpoint_checks.json)
passes every one of 1,309 payloads, all 1,222 Python files and the reconstructed
source fingerprint. Its 32,101,634-byte archive preserves the 100-fixture scene
and scoped evidence; it is not an installer. The current
[938-test/105-module run](../../experiments/streets/python_v1/suite_dev65_01/suite_checks.json)
has completed with 938 passes and no failures, errors or skips. Its
[coverage audit](../../experiments/streets/python_v1/suite_dev65_audit_01.json)
verifies exact ordered test coverage, all current modules and Python sources,
and terminal resource receipts. The single process took 3,271.73 s, with peaks
of 175,550,464 bytes resident and 191,934,464 sampled private bytes under the
0.25 GiB cap, 1 GiB reserve and one numerical thread. Native, gallery,
real-input and installer qualification remain separate.

The [reproducibility gap](../../experiments/streets/python_v1/catalogue_gap_dev65_01.json)
identified fourteen missing inputs in the original 86-fixture check. The new
repository catalogue includes all 100 and the canonical saved chain widths.
An [independent capture](../../experiments/streets/python_v1/catalogue_capture_dev65_01/catalogue_checks.json)
reads all saved source/mesh/furniture references without regeneration; the
[candidate's native run](../../experiments/streets/python_v1/gallery_generate_candidate_dev65_01/generation_checks.json)
regenerates all 100 exactly through Update, without importing an output mesh.
The [maintained repository entry point](../../experiments/streets/python_v1/gallery_generated_dev65_01/generation_checks.json)
also passes all 100 exact source/mesh/furniture comparisons on the current
tooling source, using 660,979 direct quads in total and a layout without overlaps.
It releases runtime caches between fixtures. Its
[saved-panel check](../../experiments/streets/python_v1/gallery_bootstrap_dev65_01/bootstrap_checks.json)
passes optional/idempotent registration, respect for an already registered
addon, actual cold Update and unchanged source/mesh/furniture/object identity
for the T and Smooth internal-union examples. Generation/reopening peak at
624,881,664/318,083,072 sampled private bytes under 0.75/0.5 GiB caps, retaining
the 1 GiB reserve and one numerical thread. The
[source/reproduction audit](../../experiments/streets/python_v1/gallery_reproduction_dev65_audit_01.json)
verifies both terminal jobs and all 1,222 unchanged qualified Python files.
These are reproduction and panel checks; complete gallery-quality qualification
remains separate. The
[qualification-input preflight](../../experiments/streets/python_v1/qualification_inputs_dev65_01.json)
verifies files, hashes and parser flags only; its 21 gallery and nine real-input
jobs have not run on this source.

### Latest standalone result: overlapping combined entrances

The [compound-kernel matrix](../../experiments/streets/python_v1/combined_compound_matrix_dev65_01/compound_matrix_checks.json)
passes 22/22 variants on the current unchanged runtime source: eleven inputs
in both Exact and Smooth, including T/skew/five-arm layouts, unequal widths
and crowns, edits, reversed storage, translated datum and resolutions 2/4.
Every surface uses direct quads, retains 2 m round corners and zero exposed
curbs, meets its incoming crown ceiling and complete mouth profile, and forms
one connected disk. Published float32 quads pass positive-corner/ratio-8 checks
and the complete bilinear 10% gradient bound. Entrance longitudinal bounds
checked on the construction mesh are at most 1.147%, below 3%. Both repeated
baseline meshes are exactly restored. Counts are 337–3,197 quads per patch,
including the higher-resolution cases.

The terminal matrix worker peaks at 80,498,688 bytes resident and 103,317,504
sampled private bytes under 0.25 GiB, one numerical thread and a 1 GiB reserve.
The source and six prototype-module fingerprints are recorded in the receipt.
This dev65 result remains an external prototype: normal addon integration and
dynamic profile-gradient forwarding are explicitly false/pending in that report.
The dev66 result above closes the gradient-forwarding gap on its separate source.
Bundle-body seams, original owner/endpoint mappings, actual Update/furniture,
reopening and installed-addon qualification are not established by this run.
See the [current completion tasks](STREETS_PYTHON.md#work-remaining-for-the-initial-python-deliverable)
and [geometric explanation](STREETS_PYTHON_GEOMETRY.md#combined-entrances-for-overlapping-members-prototype).

**Verified dev64 checkpoint:** source
`c228544813f2492e490cf6e2619cc18b4d2147b791ee4f3ad7eaa55a4131caef`.
The physical-port inventory fix retains a transverse Street's identity while
declaring only the Sections inside a Bundle. Source-endpoint restoration and
bounded quality continuation now enable the angular separator chart in the
addon. [Twenty-five focused tests](../../experiments/streets/python_v1/internal_angular_core_dev64_01/island_chain_checks.json)
pass both modes, including combined control/width edits, local/full/restoration,
whole float32 interfaces, furniture, remote locality, full profiles and E0.

The [98-fixture preview](../../experiments/streets/python_v1/internal_angular_native_dev64_02/python_streets_internal_crossings.blend)
retains the earlier 96 exactly and appends the two separated internal crossings
without overlap. [Actual Blender Update checks](../../experiments/streets/python_v1/internal_angular_native_dev64_02/native_internal_checks.json)
and [fresh selective reopening](../../experiments/streets/python_v1/internal_angular_reopen_dev64_01/native_internal_checks.json)
pass edits, full equality/restoration, furniture (including interior anchor
904), unrelated object identity and cold reconstruction. Native private usage
peaks at 459,718,656 bytes under a 0.75 GiB cap; reopening peaks at 330,985,472
under 0.5 GiB. System-reserve configuration remains 1 GiB. The first native
test failed because the harness did not reshape flat buffers before measuring
slopes. The corrected harness uses the same actual published buffers and no
addon change was needed.

The [independent physical/topology audit](../../experiments/streets/python_v1/internal_published_dev64_01/published_mesh_checks.json)
checks actual saved float32 buffers without mesh regeneration. Each paving
network is one component with Euler characteristic -1 (two holes); the two
grass components have Euler characteristic 2 and close those holes exactly.
The filled surface has Euler characteristic 1, a simple degree-two exterior
boundary and no cracks. All 556/652 grass boundary edges match paving edges
exactly. Full bilinear slopes are at most 8.995%; the crown ceiling and edge
ratio eight are retained. The separators contain 1,878/4,086 quads. Total
fixture counts 9,386/12,404 include their 1,600-quad control roads.
The [drawing](../../experiments/streets/python_v1/internal_drawing_dev64_02/drawing_checks.json)
uses these actual buffers; the formulas remain editable text in the geometry
guide.

The 933-test/104-module run at `suite_dev64_01` stopped at the system reserve
after 517 passes, with no test failures; its
[partial evidence](../../experiments/streets/python_v1/suite_dev64_01/partial_evidence.json)
preserves exact completed and pending IDs. The
[continuation](../../experiments/streets/python_v1/suite_dev64_resume_01/resume_checks.json)
runs only the remaining 416 in 16-test workers capped at 192 MiB, retaining
the 1 GiB reserve. All 26 workers have now completed and passed. The
[independent coverage audit](../../experiments/streets/python_v1/suite_dev64_audit_01.json)
verifies the exact discovered 933-test set, each passing once across the prefix
and batches, zero failures/errors/skips, unchanged dev64 source and all script
and resource receipts. Resumed coverage is separate from a complete
single-process run. The initial six-case failed probe remains historical
evidence. Dev65 qualifies the paired overlapping internal runs described above;
tighter requested-radius contacts and broader configurations still need work.
Full current gallery/real/installer/performance qualification is incomplete. The complete
926-test/96-fixture qualification below belongs to preserved dev63.
The [source/receipt audit](../../experiments/streets/python_v1/dev64_internal_audit_01.json)
verifies all five scoped successes, scripts, source, resource caps, branches and
the empty Git index. It confirms 1,210 exact Python files against dev63, six
changes and two additions, with existing angular pavement bridges, profile and
intersection builders byte-exact.
The [verified dev64 checkpoint](../../experiments/streets/python_v1/checkpoints/dev64_01/checkpoint_checks.json)
contains all 1,218 Python files, the 98-fixture preview and selected docs,
diagrams, evidence and resume tools. Every one of its 1,294 payloads and the
source fingerprint is verified. Archive size is 30,959,232 bytes, SHA-256
`7b6aae45c2a284d4ce317f95688c0e137d2f1c1a875c1263a9b7ec97ac1bb962`.
The archive is recovery data, not an installer. The live full regression and
remaining release gates stay separate; later document edits do not alter it.

**Verified dev63 checkpoint:** source
`d113fbe8d2456c4f5def95dee803a8d17419394fec4525e6858cfafaa5590aaa`.
Attached angular construction now handles one continuous separator or several
ordered islands, derives their component count, coordinates strip sampling
with actual rounded-cap boundary parity, and retains source curb corners and
the Smooth approximation budget. Carriageways and cap preparation share
`CarriagewayProvider.span_metadata`; cap planning requires no temporary road
mesh. Its actual port dependencies are projected to the relevant Bundle.

The [13-test focused result](../../experiments/streets/python_v1/angular_attached_core_dev63_02/island_chain_checks.json)
passes physical float32 interfaces, complete mouths, one/three/two-island edits,
exact local/full/restoration, reversed storage, numerical R2, translated datum
with E0, tighter curb fidelity, and both positive/negative dependency projections.
It repairs the separated-member anchor regression in the previous 923-test run
by consuming a separator digest only for plans that query the completed cap.
Peak focused resident usage is 75,010,048 bytes under a 0.25 GiB cap and 1 GiB reserve.

The [96-fixture preview](../../experiments/streets/python_v1/angular_states_native_dev63_02/python_streets_separator_states.blend)
contains six new Exact/Smooth examples: continuous separator, three islands,
and tighter 20 m inner-radius controls. Their Bundle paving/grass totals are
6,672/7,843, 4,465/5,396 and 2,292/3,580 quads respectively, excluding junctions,
crossroad Sections and control roads. The previous 90 saved source/mesh/furniture
sets remain exact; the new examples do not spatially overlap them.
[Actual native operations](../../experiments/streets/python_v1/angular_states_native_dev63_02/native_states_checks.json)
and [fresh selective reopening](../../experiments/streets/python_v1/angular_states_reopen_dev63_01/native_states_checks.json)
pass combined control/width edits, button-driven one/three/two-region updates,
full equality/restoration, furniture, incremental remote identity and cold
reconstruction. The initial native harness incorrectly required a remote
pointer to survive a deliberately requested force-full recreation. That failed
receipt remains; the successful harness checks immediate pre/post identity on
each incremental edit and retains all geometry/full/quality checks.
Native resident/private peaks are 394,932,224/447,352,832 bytes under a 0.75 GiB
cap and the unchanged 1 GiB reserve. Selective reopening uses a 0.5 GiB cap.

The [independent published-buffer audit](../../experiments/streets/python_v1/angular_states_physical_dev63_01/published_states_checks.json)
passes slopes, crown, whole float32 edges, and bidirectional physical contours,
area and connected-component/Euler checks against source-ribbon Boolean unions.
Exact rim discrepancy is at most 2.56 micrometres; Smooth discrepancy at most
3.83 mm. The numerical Exact test uses 0.1 micrometre; published checks add the
existing 20-micrometre float32 allowance. The Smooth physical-curve budget
remains 5 mm. The maximum longitudinal bound is 2.99684%, under 3%.
No output mesh defines the footprint reference and no mesh is regenerated by
this audit. The [six-panel drawing](../../experiments/streets/python_v1/angular_states_drawing_dev63_03/drawing_checks.json)
also reads actual saved buffers; [formulas](STREETS_PYTHON_GEOMETRY.md#source-preserving-entrance-samples-and-cap-parity)
remain editable LaTeX text.

The complete [926-test/103-module regression](../../experiments/streets/python_v1/suite_dev63_01/suite_checks.json)
now passes: zero failures, errors or skips, with the exact source fingerprint
unchanged before/after. Peak resident/private memory is
176,381,952/191,303,680 bytes under a 0.25 GiB process cap and the unchanged
1 GiB reserve. The separated-member anchor regression is covered by this full
run as well as the focused checks.
The [serial 96-fixture Blender qualification](../../experiments/streets/python_v1/gallery_dev63_qualification_01.json)
also passes all 20 jobs. Its [receipt/source audit](../../experiments/streets/python_v1/dev63_gallery_receipt_audit_01.json)
verifies all 96 current-source cold and forced-full reconstructions, exact saved
source/mesh/furniture data, published convexity, crown bounds, overlap and edge
ratio checks. All 22 closed angular pavement/grass interfaces use identical
whole float32 edges. The gallery has 620,608 quads in total across 96 networks;
this is not a per-Bundle count. Maximum edge ratio is 7.99635, within 8.
Peak resident/private usage is 284,299,264/341,659,648 bytes under the 0.5 GiB
process cap and unchanged 1 GiB reserve. Delivery/reference scene bytes remain
exact. This qualification reads the R1 examples; higher-resolution native
checks, current real-network/installer qualification and the remaining
construction/performance gates stay separate. The broader goal is incomplete.
The [increment receipt/source audit](../../experiments/streets/python_v1/dev63_increment_audit_01.json)
checks six terminal successes, exact scripts/source, caps/reserve,
document links, branches and an empty staging index. It confirms 1,206 exact
Python files against dev62, eight changes and two additions limited to the
backend/tests. The older angular surface/island builders are byte-exact.
The [dev63 recovery archive](../../experiments/streets/python_v1/checkpoints/dev63_03/streets_python_dev63_recovery.zip)
and [verification](../../experiments/streets/python_v1/checkpoints/dev63_03/checkpoint_checks.json)
preserve 1,342 payloads, including all 1,216 repository Python files, documents,
corrected drawings, the 96-fixture preview and complete regression/gallery
receipts. Every payload and the exact Python fingerprint is checked.
The archive is 29,326,260 bytes, SHA-256
`d2118b455ff4738e806a326c106e0815381edccf36ac7de9b5cf7390ac6d3111`.
It is not an installer or a whole-goal completion claim; historical evidence
and the current regression remain separate. Later document edits do not alter it.

**Previous dev62 development:** source
`2bf4a8fd1425ed0929177822e9a0c845e1d1fe70cb2718a8d6d9ec74cea49a84`.
The new `AttachedAngularPlan` fallback generates the previously rejected
curved, overlapping Bundle with four actual rounded T-junction connections.
Two single-Section members retain their original paths, widths, common crown
and complete entrance profiles. Two grass islands join the actual junction
and crossroad boundary samples; every mesh face is a directly constructed quad.
This is the focused attached-angular family, not all possible attached Bundles.

The [six permanent-test qualification](../../experiments/streets/python_v1/angular_attached_core_dev62_01/attached_core_checks.json)
passes Exact/Smooth float32 interfaces, complete mouth partitions, reversed
storage, resolution 2, local/full vertex and width Updates/restoration,
remote/furniture preservation and physical queries at zero exaggeration on a
translated datum. Peak resident usage is 68,730,880 bytes under a 0.25 GiB cap
and the unchanged 1 GiB reserve. Furniture now explicitly depends on its Bundle
separator, including actual cap height reconstruction.

The [native publication/edit check](../../experiments/streets/python_v1/angular_attached_native_dev62_01/native_attached_checks.json)
and [fresh-process reopening](../../experiments/streets/python_v1/angular_attached_reopen_dev62_01/native_attached_checks.json)
pass actual Update-button edits, exact full rebuilding/restoration, remote
identity and cold reconstruction. The
[combined preview](../../experiments/streets/python_v1/angular_attached_native_dev62_01/python_streets_attached_angular.blend)
contains **90 fixtures**; the previous 88 retain exact saved source, geometry
and furniture and the new two do not spatially overlap them. New Bundle
pavement/grass totals are **3,884 Exact / 4,867 Smooth quads**, excluding the
four junction meshes, crossroad Sections and control road. Native resident/
private peaks are 374,603,776/450,011,136 bytes under a 0.5 GiB cap and 1 GiB
reserve. Initial edits take 0.901/1.279 s; reopened edits take 0.931/1.177 s.
These timings describe those edits, not universal responsiveness.

The [saved float32 audit](../../experiments/streets/python_v1/angular_attached_slopes_dev62_01/published_slope_checks.json)
passes pavement slopes 4.113%/3.972%, longitudinal bounds 2.946%/2.578%, the
incoming crown ceiling and all 328/380 whole separator boundary edges. Grass
slopes are 7.966%/6.231%, within the same 10% physical limit. The audit never
regenerates meshes. The [illustration record](../../experiments/streets/python_v1/angular_attached_drawing_dev62_01/drawing_checks.json)
likewise uses actual saved buffers; see [construction/formulas](STREETS_PYTHON_GEOMETRY.md#curved-bundles-attached-to-real-junctions).
The [increment audit](../../experiments/streets/python_v1/dev62_increment_audit_01.json)
checks the five terminal jobs and scripts, exact source, memory caps/reserve,
document links, branch separation and empty staging index. It confirms 1,205
unchanged Python files against dev61, with four changes and five additions
limited to the backend/tests; the old angular surface/island builders are exact.

The complete [923-test, 103-module regression](../../experiments/streets/python_v1/suite_dev62_01/suite_checks.json)
is terminal with one anchor-dependency failure and no errors, on unchanged
source. Dev63's focused checks above resolve that failure. Dev62 therefore did
not pass the full release gate; it also did not establish that all 88 older
fixtures were regenerated on this source or that real-import/installer gates
were complete. Internal
junctions, more general attached topology/profile changes, mixed bent joins
and expensive Updates remained open at that checkpoint. Its regression worker
is terminal; the current dev63 run is recorded separately above.

The [dev62 recovery archive](../../experiments/streets/python_v1/checkpoints/dev62_01/streets_python_dev62_recovery.zip)
and its [verification](../../experiments/streets/python_v1/checkpoints/dev62_01/checkpoint_checks.json)
preserve 1,286 payloads including all 1,214 repository Python files, current
documents/drawings, the combined 90-fixture preview and selected resume/evidence
files. Every payload and the reconstructed source fingerprint is verified.
The archive is 27,031,945 bytes, SHA-256
`189db5ff4b22ba2975d7e78b3df20647811dd13a0c922243fe7d1491d444ebf3`.
It is not an install-qualified or published release. Full historical evidence
and the running regression remain outside it. Later documentation edits do
not change the frozen archive.

**Previous dev61 development:** source
`5341dcbd5018646a68a50908806cbccab232fe7e2a3726170df3d40bc579244e`.
The [19-test focused qualification](../../experiments/streets/python_v1/angular_island_core_dev61_01/island_chain_checks.json)
passes Exact/Smooth three-island float32 construction, whole interfaces,
vertex/width local-full equality, exact restoration, remote furniture and mesh
preservation, reversed input storage, two/four-island examples and the existing
union/spatial-index tests. Peak resident usage is 84,320,256 bytes under the
same 1 GiB process cap/reserve. The
[four-island Smooth result](../../experiments/streets/python_v1/angular_four_islands_dev61_02/multi_island_checks.json)
has 7,907 pavement and 476 grass quads; a bounded local layered chart avoids the
previous 9,029-quad whole-body fallback on one half. The original source layout,
curve tolerance, shape/crown/slope limits and overlap checks are retained.
The initial development failures and corrected test-harness error remain in
their original receipts. These focused results alone do not establish native
Blender, the complete current suite, real-network or installed-copy qualification.

The separate [native qualification](../../experiments/streets/python_v1/angular_island_native_dev61_03/native_island_checks.json)
and [fresh-process reopening](../../experiments/streets/python_v1/angular_island_reopen_dev61_01/native_island_checks.json)
now pass both three-island examples: actual Update-button vertex/width edits,
exact full rebuilding and restoration, furniture dependency updates, unchanged
remote meshes/poses and cold reconstruction. The
[combined scene](../../experiments/streets/python_v1/angular_island_native_dev61_03/python_streets_island_chains.blend)
contains 88 fixtures, with exact saved source/meshes/furniture for all 86 prior
fixtures and no spatial overlap of the new examples with them. New Exact/Smooth
Bundle counts are 7,886/8,386 quads, excluding their control roads. Peak native
resident/private usage is 359,743,488/411,897,856 bytes under the unchanged
1 GiB cap/reserve. Native R1 quality uses the original published-buffer auditor;
only its SAT batch count is excluded from equality after a force-full changes
object iteration order. All geometry comparisons and acceptance limits remain.

The [independent saved-buffer slope audit](../../experiments/streets/python_v1/angular_island_published_slopes_dev61_01/published_slope_checks.json)
passes 2.983%/2.754% longitudinal bounds, physical slopes under 10%, the
incoming crown ceiling and exact whole float32 grass/pavement edges, without
regeneration. The [density boundary audit](../../experiments/streets/python_v1/angular_island_density_boundaries_dev61_01/density_boundary_checks.json)
also proves complete physical polygon coverage in both directions and identical
grass arrays for the four-island comparison. A complete current regression,
fresh real imports, installer, and higher resolutions of the new long examples
are still separate gates. Broader
attached/internal junctions, profile changes, mixed bent interpolation and
expensive Updates remain construction/performance work before the final release
qualification. The dev59 gallery remains the last complete 86-fixture audit.
The [increment receipt/source audit](../../experiments/streets/python_v1/dev61_increment_audit_01.json)
verifies all six terminal successful jobs and their actual scripts, limits and
source hashes. Against the dev58 recovery source, 1,199 Python files are exact;
the five changes and five additions are confined to the numerical street
backend/tests. The old angular surface builder is byte-exact. This is scoped
preservation evidence, not a substitute for the complete current regression.

The [dev61 recovery archive](../../experiments/streets/python_v1/checkpoints/dev61_01/streets_python_dev61_recovery.zip)
preserves 1,288 payloads, including all 1,209 repository Python files, the
88-fixture preview, resume helpers and selected receipts. Its
[check](../../experiments/streets/python_v1/checkpoints/dev61_01/checkpoint_checks.json)
verifies every payload and reconstructs the exact source fingerprint. It is
25,690,410 bytes, SHA-256
`660d4701ba8eea0ee39c680aba1468b5a99d5e166d214b093237e1a07acbcb04`.
This is a local recovery checkpoint, not an installer or a completed release.
The full `suite_dev61_01` run stopped at the system-memory reserve after
**687 of 917 tests** completed, with no recorded test failures. Its
[partial evidence and remaining test IDs](../../experiments/streets/python_v1/suite_dev61_01/partial_evidence.json)
preserve those observations against the unchanged source above; **230 tests
remain**, including the interrupted test. Peak worker resident/private usage was
172,810,240/182,247,424 bytes. The stop was a system-wide reserve event, not a
worker-cap overrun or a passing full-suite result. The first remaining-tests
attempt did not start because preflight memory was insufficient. Use
`experiments/streets/python_v1/run_remaining_suite_dev61.py` and the preserved
plan for a bounded continuation; it records each completed test separately.
No worker from those attempts remains active. Later documentation edits do not
alter the recovery archive.

The next attached-angular input is authored in
`experiments/streets/python_v1/attached_angular_fixture.py` outside the addon.
Its [source-only check](../../experiments/streets/python_v1/attached_angular_inputs_dev61_01.json)
passes original curve regularity, continuous Bundle-centerline separation and
four coincident three-arm junction incidences in both modes. The
[actual mouth preparation](../../experiments/streets/python_v1/attached_angular_ports_dev61_01/attached_port_checks.json)
also passes without building surfaces. A
[baseline diagnostic](../../experiments/streets/python_v1/attached_angular_baseline_dev61_01/attached_angular_checks.json)
confirms that the current production construction rejects this valid input;
successful diagnostic execution is not successful network generation.

Two external prototypes now pass focused numerical checks in both modes:
[the four full-profile entrance connectors](../../experiments/streets/python_v1/angular_bridge_candidate_dev61_01/angular_bridge_checks.json)
retain exact original mouth positions, direct quads and the existing slope
limits; [the two open-ended grass body pieces](../../experiments/streets/python_v1/angular_open_candidate_dev61_01/open_island_checks.json)
retain the physical source curbs within the original tolerance and pass shape
and overlap checks. The latter check peaks at 42,536,960 resident bytes under
a 0.25 GiB process cap with the unchanged 1 GiB system reserve. These prototypes
are outside the addon. Joining the pavement body, grass pieces and actual
rounded junction caps, then integrating Update dependencies and qualifying
complete native output, remains undone. They do not qualify an attached network.

**Previous dev60 boundary checkpoint:** source
`f92ac0f8c040c0de571fe68dc9dd4cb41c7ccc1c033ac9ae48c613409ab639e8`.
The [14-test boundary qualification](../../experiments/streets/python_v1/region_pair_core_dev60_02/region_pair_checks.json)
passes at 40,710,144 peak resident bytes. Conservative spatial pruning handles
arrangements that previously exceeded the exhaustive pair estimate, preserving
the exact samples/intersection formulas and existing work/segment limits.
The [actual Smooth Bundle probe](../../experiments/streets/python_v1/angular_multi_spatial_dev60_02/multi_island_checks.json)
now obtains one outer boundary and three separator holes from 4,128 samples;
it still fails at missing multi-island mesh decomposition. Its successful
diagnostic process exit does not mean generation succeeds. No dev60 native
gallery, full suite, fresh real-network or installer gate is claimed.

**What remains to reach the goal:**

1. Complete broader valid Bundle construction: attached Street connections and
   internal junctions, general interior islands and topologies beyond disjoint
   ordered angular island intervals. Dev62 now qualifies the two-member,
   single-Section attached-angular family with four actual mouths. Dev63 extends
   it to one continuous separator, three ordered islands and tighter entrances
   retaining source curbs. Other attached arrangements still need construction.
   The [internal-crossing inventory probe](../../experiments/streets/python_v1/internal_inventory_dev63_01/inventory_checks.json)
   captures an incidence-only rejection: transverse Street 800 has exterior
   Sections 300/302 and interior Section 301 between physical junctions 500/501.
   The interior inventory correctly declares only 301, but the dev63 validator
   requires the complete global Street inventory. Changing only 301's Street ID
   makes the same incidence acceptable. Dev64 now fixes that inventory rejection
   while preserving the original identity. The separate six-case surface probe
   still fails; inventory acceptance is not proof of surface-generation support.
   Dev61 above constructs the focused
   multiple-island family. Changing-width angular examples are qualified in
   dev59 below; broader constructions still require their own evidence.
2. Complete outstanding profile-family changes at junctions, multiply claimed
   dimensional seams and bent mixed Exact/Smooth joins. Preserve valid-path
   and curvature checks for other sharp Smooth cases.
3. Reduce expensive affected-region Updates with measured gains on identical
   inputs and preserved mesh quality/dependency behavior.
4. Qualify one final source across the complete regression suite, combined
   Blender scene, independent real OSM inputs, Update/full equivalence,
   furniture dependencies, save/reopen, slopes/crowns, density and the 1 GiB
   process cap. Build and check the installer/documentation from that source.

The focused attached-angular Bundle is now implemented and qualified as dev62
above. Broader construction and complete regression of the combined
implementation remain separate gates. The previous gallery and real-network
receipts below retain their original source scope.

**Qualified dev59 gallery checkpoint:** source
`c04f837fc5a81a841a6d421cdcd66cdcc9bd1f08fc3f4e51bde510a61b162f9c`
adds fixed-boundary grading after separator refinement. The prior queue is
terminal failure at [Exact180 R4](../../experiments/streets/python_v1/angular_selective_resolution_01_dev58_01/failure.txt):
four grass quads exceed ratio eight after subdivision. The
[scoped native experiment](../../experiments/streets/python_v1/refined_separator_native_01/refinement_checks.json)
repairs 9.53067 to 7.58178 in float32 without moving the pavement/grass boundary
or changing connectivity/heights. Its peak resident usage is 315,535,360 bytes.
The [12-test focused check](../../experiments/streets/python_v1/separator_refinement_core_dev59_01/refinement_core_checks.json)
now passes on unchanged source, zero failures/errors, with 141,787,136 peak
resident bytes. The [actual Blender retry](../../experiments/streets/python_v1/angular_selective_resolution_01_dev59_01/case_checks.json)
now passes R2/E0 and R4/E10, exact full rebuilding and restoration of source,
meshes and furniture, and the unchanged remote network. Peak resident/private
usage is 409,124,864/474,464,256 bytes. The
[Smooth225 native R2/R4 sequence](../../experiments/streets/python_v1/angular_selective_resolution_02_dev59_01/case_checks.json)
also passes, peaking at 430,051,328/483,631,104 resident/private bytes.
The [current gallery qualification](../../experiments/streets/python_v1/gallery_selective_dev59_qualification_02.json)
now passes all 86 fixtures and 535,485 published quads, including control roads.
All 74 unchanged fixtures retain exact reference source/mesh/furniture; twelve
inputs differ only by the intended width changes. All current saved buffers and
furniture survive upgrade, cold reconstruction and forced full rebuilding.
Twelve closed angular interfaces match whole float32 edges exactly. Maximum
edge ratio is 7.99635; gallery peak resident/private usage is
281,169,920/341,639,168 bytes. The
[receipt audit](../../experiments/streets/python_v1/dev59_gallery_receipt_audit_01.json)
verifies all eighteen workers and both selected resolution sequences. Current
full-suite, real-import, installer and broader construction/performance gates
remain incomplete;
the passing 904-test result
below belongs to the previous source, not retroactively to dev59.

**Previous dev58 development:** source
`ac25fe4cf0007fc847715d339028d354f3e7b2b8b2888e0274f7fd4b5e6f6391`
adds a physical-slope-triggered interior quad repair. The
[captured native failure](../../experiments/streets/python_v1/angular_selective_edit_10_dev57_01/failure.txt)
is a +0.1 m move in the Exact 270-degree changing-width island. A nearly
collinear quad amplifies the original profile to a measured 12.2625% slope;
the [fresh Blender diagnostic](../../experiments/streets/python_v1/angular_edit_slopes_native_01/diagnostic_checks.json)
reproduces it without incremental caches. The original field sampled in that
cell has about 2.8% slope. The
[scoped native repair](../../experiments/streets/python_v1/angular_edit_shape_native_01/diagnostic_checks.json)
reduces the measured slope to 4.5688%, keeping exact boundary XYZ and all
connectivity. This experiment is evidence for the repair, not completed addon
qualification. The
[integrated Blender edit](../../experiments/streets/python_v1/angular_selective_edit_10_dev58_02/case_checks.json)
now passes actual vertex/width Update operations, exact local/full meshes and
furniture, restoration, cold reconstruction and preservation of unrelated
outputs. Both appended networks retain their saved source/mesh/furniture after
the version upgrade. Peak resident/private usage is 292,651,008/361,639,936
bytes. The [Smooth270 edit sequence](../../experiments/streets/python_v1/angular_selective_edit_11_dev58_01/case_checks.json)
also passes, peaking at 295,870,464/364,388,352 resident/private bytes. The
[serial dev58 queue](../../experiments/streets/python_v1/profile_shape_dev58_qualification_01.json)
now passes the [complete 904-test regression](../../experiments/streets/python_v1/suite_dev58_01/suite_checks.json)
in 100 modules, with zero failures/errors and matching before/after source
fingerprints. Peak resident/private usage is 165,425,152/187,260,928 bytes.
The first native R2/R4 sequence passes; the second stops on the refined
separator shape above. The third resolution sequence and ten other native
edits did not run. Combined saved-quality/reopening,
real-import and installer qualification remain separate later gates.
The [dev58 source recovery record](../../experiments/streets/python_v1/checkpoints/dev58_01/checkpoint_checks.json)
also verifies the archived code/assets and nine resume helpers. Recovery
verification does not substitute for any pending qualification gate.

**Previous dev57 development:** local physical quad tiles are integrated on source
`444a153c8d4eeff1cf3ae5e62f7db711a71e6aabeca6a6164ac2f57b45c2a569`.
The [40-test focused run](../../experiments/streets/python_v1/angular_tiles_core_dev57_01/tile_checks.json)
now passes all new Exact float64/float32 inputs, the captured width-edit/R4
regressions, independent physical shared-edge/curve checks and old angular/axial
contracts. Peak resident/private usage is 143,855,616/143,355,904 bytes. The
[new 86-fixture generation preview](../../experiments/streets/python_v1/angular_width_gallery_dev57_01/python_streets_width_changes.blend)
has all twelve width-change cases saved after actual Update operations. Its
[generation record](../../experiments/streets/python_v1/angular_width_gallery_dev57_01/generation_checks.json)
passes; the [original serial qualification](../../experiments/streets/python_v1/angular_tiles_dev57_qualification_01.json)
stopped at the system-memory reserve after seven completed native edit cases.
The [selective continuation](../../experiments/streets/python_v1/angular_native_resume_dev57_01.json)
completes three more, then captures the 270-degree Exact failure above. The
smaller workers append only the tested fixture and an unrelated reference
from the unchanged combined scene. This lowers the observed working set:
the first complete selective edit peaks at 273,358,848 resident bytes.
The original full regression did not start. Initial generation is not complete editing or release
qualification. The new classes and continuation points are described in the
geometry document. A separate local-first work diagnostic is external only;
its concurrent timings are not a qualified native performance claim. Its expanded
16-input matrix generates every input but changes several mesh buffers, so
that dispatch shortcut has not been adopted.

**Previous dev56 development:** source-anchored sampling passes the
[37-test focused regression](../../experiments/streets/python_v1/angular_sampling_core_dev56_01/sampling_checks.json)
on `59f32c84ff695bf173f6f612181d8843d5b85da665b42c0c5a9b2f49b2abb622`.
All three changing-width Smooth islands generate from float32 coordinates;
the six complete-overlap cases also retain their contracts. Source sampling
phase stays unchanged away from a trim and allocation budgets reject excessive
sampling before allocation. This is nine initial constructions, not nine
fully qualified editing configurations. The
[86-fixture inspection scene](../../experiments/streets/python_v1/angular_width_gallery_dev56_01/python_streets_width_changes.blend)
saves all nine initial Update results. The following native edit job fails the
180-degree Smooth island's +0.1 m first-member width edit after six complete-
overlap edit sequences pass. The R4 Smooth island and three Exact islands also
remain open. No new installer qualification or full dev56 regression is claimed.
The [independent saved-mesh audit](../../experiments/streets/python_v1/angular_width_gallery_quality_dev56_01/import_quality_checks.json)
now passes all 86 networks and 529,023 quads: maximum ratio 7.9963482749, zero
cells above ratio 8, with 393,900,032 peak resident bytes. This qualifies the
saved R1 geometry, not the failing width-edit or R4 inputs. The broader goal
remains active.

**Previous dev55 development:** finite changing-width curbs now participate in the
angular Bundle provider with original Section/Connection IDs, certified
polynomial boundaries and composite transition ownership. The
[before diagnostic](../../experiments/streets/python_v1/angular_dimensions_before_01/dimension_probe.json)
records all twelve authored inputs rejected by dev54. Dev55's
[current construction check](../../experiments/streets/python_v1/angular_dimensions_coverage_01/dimension_checks.json)
is still failing as a whole: seven cases generate, including all six complete
overlaps and the 225-degree Smooth island; five island layouts remain
unsupported. These tests are development
evidence, not passing acceptance. The complete twelve-case construction,
independent curve/corner audits, Update/local-full/furniture checks and broad
native regression are required before qualifying this increment. The
[33-test focused check](../../experiments/streets/python_v1/angular_dimensions_core_01/dimension_core_checks.json)
passes new overlap/ownership/Update/query contracts and older angular/axial
regressions on `b2e91ad20ccd...`. A separate
[float32 before regression](../../experiments/streets/python_v1/angular_dimensional_float32_before_01/float32_checks.json)
fails two Exact seams; its
[repair check](../../experiments/streets/python_v1/angular_dimensional_float32_after_01/float32_checks.json)
passes all six inputs on `668f7f803169...`. The existing ten-micrometre limit
is reused and the actual boundary discrepancy is reserved, not hidden.
The [current nine new-contract tests](../../experiments/streets/python_v1/angular_dimensions_new_01/dimension_core_checks.json)
pass on `a92f0427585a...`, including all six float32 source cases, width/vertex
edits, returning to compatible widths, local/full meshes/furniture, restoration
and released attachment queries. Peak resident/private usage is
65,507,328/83,791,872 bytes. The remaining
[partial-core qualification](../../experiments/streets/python_v1/angular_dimensions_dev55_partial_qualification_01.json)
has completed on that frozen source. It distinguishes passing core checks from the separate
twelve-case construction gate; an unfinished or failing island is not support.
The [full regression](../../experiments/streets/python_v1/suite_dev55_01/suite_checks.json)
now passes 894 tests in 98 modules, zero failures/errors, on unchanged
`a92f0427585a...`. Its peak resident/private usage is
162,877,440/188,792,832 bytes. The full new construction and native gates remain
separate from this passing regression. All three existing-gallery native gates
now pass: [upgrade](../../experiments/streets/python_v1/angular_dimensions_old_gallery_dev55_01/upgrade_checks.json),
[saved quality](../../experiments/streets/python_v1/angular_dimensions_old_quality_dev55_01/import_quality_checks.json)
and [fresh reopening/editing](../../experiments/streets/python_v1/angular_dimensions_old_reopen_dev55_01/reopen_checks.json).
All 86 fixtures retain their exact source/mesh/furniture, with 525,028 quads.
Maximum native resident/private usage is 519,278,592/592,629,760 bytes.
The separate [seven new native R1 edit sequences](../../experiments/streets/python_v1/angular_dimensions_native_demo_02/cases_partial.json)
pass; that broader job then fails its R4 generation guard. The failed receipt
is preserved and no release qualification is inferred from those seven rows.
All earlier
native scenes and installers remain unchanged. The broader goal stays active.

**Qualified dev54 native increment:** the worktree adds exact invalid-corner rejection to the
call-local quad energy workspace, without changing finite energies, gradients
or acceptance guards. The
[focused integrated check](../../experiments/streets/python_v1/integrated_quad_witness_dev54_01/integrated_witness_checks.json)
passes 600 exact energy comparisons, six captured actual descents and nine
workspace tests on source `5e22994068286de9d4b5e76b263499d5534e134832dae783dcca64165921e61c`.
The [test-only native comparison](../../experiments/streets/python_v1/canal_quad_witness_dev53_01/witness_update_comparison.json)
passes eight interleaved same-process measurements: median whole Update time
8.5096 s reference / 7.3316 s candidate (13.84% lower). Both methods retain exact
edited/local-full/restored buffers, source records and the other component;
14/20 output objects remain. Candidate peak resident usage is 279,859,200 bytes.
The [integrated full regression](../../experiments/streets/python_v1/suite_dev54_01/suite_checks.json)
now passes **885 tests in 96 modules**, zero failures/errors, on unchanged
source. Its peak resident usage is 163,074,048 bytes. The separate
[integrated runtime comparison](../../experiments/streets/python_v1/canal_quad_witness_dev54_01/witness_update_comparison.json)
now passes eight interleaved measurements: 8.5162 s reference / 7.4205 s
production median wall time, **12.87% lower**, and 12.70% lower median CPU time.
Edited/local-full/restored buffers and original source records are exact;
14/20 outputs plus the other component remain. Its peak resident/private usage
is 287,305,728/355,299,328 bytes. The
[native scene upgrade](../../experiments/streets/python_v1/quad_witness_native_dev54_01/upgrade_checks.json)
passes all 86 fixtures with exact source/mesh/key/furniture preservation. Its
[independent quality/reference check](../../experiments/streets/python_v1/quad_witness_quality_dev54_01/import_quality_checks.json)
passes all 525,028 direct quads, no overlaps, ratio/crown guards and exact dev53
buffers. Native peak resident/private usage is 468,660,224/530,558,976 bytes.
The [fresh-process reopening](../../experiments/streets/python_v1/quad_witness_reopen_dev54_01/reopen_checks.json)
now also passes all 86 controllers, both twelve-case angular families, actual
vertex/width Update operations, exact local/full meshes/furniture, restoration
and R2/R4 publication. All four jobs in the
[serial qualification](../../experiments/streets/python_v1/quad_witness_dev54_qualification_01.json)
are terminal success. The
[independent audit](../../experiments/streets/python_v1/dev54_witness_receipt_audit_01.json)
verifies these plus the focused/full regression receipts, scripts and frozen
source. Only the energy workspace, its test module and version differ from
dev53's 1,196 Python files. Maximum resident/private usage is
521,256,960/585,773,056 bytes under the unchanged 1 GiB cap/reserve. These checks
qualify the optimization and gallery, not a new installer or three fresh OSM
pipelines on dev54. Dev53 remains the fully qualified
installable preview below. The broader goal remains active.

**Qualified dev53 preview:** this increment adds compatible multi-Section angular
routes without replacing the original paths. The
[focused native run](../../experiments/streets/python_v1/angular_chains_focused_native_02/angular_chain_native_checks.json)
now passes all twelve new Update/restore/cold-recovery sequences and actual
R2/E0, R4/E10 publication on source `841424fe0c2f...`. Its maximum resident usage
is 364,007,424 bytes. The
[standard combined native check](../../experiments/streets/python_v1/angular_chains_native_04/scene_checks.json)
now passes all 86 fixtures and the twelve new button-edit sequences. Its
[independent saved quality/reference check](../../experiments/streets/python_v1/angular_chains_quality_04/import_quality_checks.json)
passes all 525,028 direct quads, no overlaps and exact unchanged positions and
connectivity for the previous 74 fixtures. Native peak resident/private usage
is 556,879,872/616,906,752 bytes.
[Fresh reopening](../../experiments/streets/python_v1/angular_chains_reopen_04/reopen_checks.json)
also passes all 86 controllers, both angular fixture families, button edits,
restoration and R2/R4 publication. The
[full regression](../../experiments/streets/python_v1/angular_chains_suite_04/suite_checks.json)
passes **882 tests in 96 modules**, zero failures/errors, on the identical
`841424fe0c2fcfe93ca9f7da80f345ddbdba201b30e8ccf65b1599e2251df856`
source. The [detached terminal audit](../../experiments/streets/python_v1/dev53_angular_receipt_audit_01.json)
verifies all four actual jobs, hashes and unchanged cap/reserve. The
first native attempt caught an internal-junction dispatch regression; that
incidence now retains its established provider. Dev53's
[nine fresh real-network jobs](../../experiments/streets/python_v1/dev53_real_import_checks_01.json)
now also pass: Berlin/Canal/Nashville have COMPLETE coverage, two READY
components each, identical input hashes and exactly unchanged dev52 mesh buffers.
Fresh seam edits, local/full equality, restoration and unrelated-output identity
pass. The [detached real audit](../../experiments/streets/python_v1/dev53_real_receipt_audit_01.json)
verifies all nine terminal receipts and their hashes. Maximum resident/private
usage is 656,871,424/709,816,320 bytes, below the same cap/reserve. Local edits
take 1.53/8.40/0.87 seconds; these are individual timings, not a controlled
speedup comparison.

The [dev53 addon ZIP](../../experiments/streets/python_v1/packages/dev53_v1/blosm_streets_python_dev53.zip)
now passes [all four installed-copy checks](../../experiments/streets/python_v1/package_dev53_checks_01.json).
The actual ZIP installs and fully registers/unregisters in an isolated profile;
both backend choices/external materials, 31 selected current fixtures, actual
manual Update/local-full/restoration, normal Canal import and fresh seam edits
pass using installed code only. Its
[detached qualification](../../experiments/streets/python_v1/package_dev53_qualification_01.json)
checks all 1,246 archive/installed payloads, the 1,196 Python files, source
fingerprint, actual exits and unchanged thread/memory settings. Maximum
resident/private usage is 833,974,272/899,420,160 bytes under 1 GiB. The archive
is 18,715,490 bytes, SHA-256
`e6043689ff8bbafc5305519b4e77035caa81428cfaeb42483581a0d874d1b8d1`.
The installed Canal edit takes 8.45 s and retains 14/20 outputs plus the other
component. Original preferences are not saved; the ZIP has not been published.
Its internal creation receipt stays at the immutable build-time status;
the detached installed-copy qualification supplies the completed result.
The broader geometry/performance goal remains active.

**Previous qualified checkpoint — dev52.** Dev52 passes the standard
74-fixture native/quality/reopening checks and all nine fresh real-network jobs.
The current d8d9 tooling source has exact runtime/unaffected-test equivalence
with the f34 874-test checkpoint and passes all seven affected numerical tests.
Its [preview ZIP](../../experiments/streets/python_v1/packages/dev52_v1/blosm_streets_python_dev52.zip)
now also passes all four installed-copy checks and the
[detached archive/payload audit](../../experiments/streets/python_v1/package_dev52_qualification_01.json).
Nineteen current fixtures, both backends/external materials, actual Update,
normal Canal import, saved quality and fresh seam edits pass using installed
code only. Peak resident/private usage is 799,191,040/864,784,384 bytes under
the unchanged 1 GiB cap and system reserve. The ZIP is 18,633,993 bytes,
SHA-256 `b6b1644f42246072ce8c6ab597d3cf5900796d894a50487f9ca0ed278f3b2345`.
It is a qualified local preview, not a published release or completion of the
broader construction/performance work.
The preceding complete core/real checkpoint was dev51, whose
shared-road-sampling integration passes all nine core and real-network jobs,
including 863 tests. Its corresponding installed-copy-qualified package was
dev49 below. The preceding dev50 core passes all six jobs, including 862
tests in 93 modules, 62-fixture saved quality, fresh edits/reopening and independent
compound/spacing audits. The installed-copy-qualified dev49 geometry and
regression checkpoint passes all 855 tests in 92 modules, the 61-fixture
gallery, fresh Update/reopen checks and all three real-network checks. Its
source is `d4ed8df50b4d699aa3a1147ef4391efae799d91be3ca7cafafe026f2d7c3729a`.
The [dev49 archive](../../experiments/streets/python_v1/packages/dev49_v1/blosm_streets_python_dev49.zip)
also passes [all four installed-copy checks](../../experiments/streets/python_v1/package_dev49_checks_01.json),
including six close-junction fixtures, actual installation/full registration,
both backend choices, external materials, normal OSM import and fresh Update.
The gallery has 362,283 convex direct quads and zero GN modifiers. Every
installed runtime module originates in that package; its source and all payloads
remain unchanged. The detached package audit verifies script hashes, terminal
exits and the 1 GiB resident/private cap with the same system reserve. Its
largest installed-copy worker uses 710,852,608 resident / 768,606,208 private bytes.
Broader construction work below remains unfinished. Earlier
checkpoints and versioned failures are retained as historical evidence.

### Dev52 integrated angular branch: focused and native checks pass

The first integrated angular milestone is `0.1.0.dev52`, source
`f34bbc9972b356adcd04551974a231dd3258ce53227aaa5d9714e3d093ec3f40`.
The preceding dev51 summary describes the last complete core/real qualification,
not the newer in-progress source. Dev52 replaces the external-only status for
the initial free-end, two-single-Section angular construction. Runtime classes
read dimensions, traversal and identities from Blosm and retain the existing
provider/dependency/Update/attachment contracts.

[Eleven focused tests](../../experiments/streets/python_v1/angular_runtime_core_07/angular_runtime_checks.json)
pass all 12 reference inputs and additional width/offset/rise/datum/orientation,
resolution 2/4, zero exaggeration, continuous-contact rejection, local/full,
restoration and released-query checks. The
[native Blender check](../../experiments/streets/python_v1/angular_native_01/angular_native_checks.json)
passes all new source edits and preserves the 62 previous fixtures exactly in
one 74-controller scene, using at most 476,282,880 resident bytes. The
[fresh-file operator/reopen check](../../experiments/streets/python_v1/angular_native_reopen_01/angular_reopen_checks.json)
passes all 12 new fixtures. The
[all-gallery quality audit](../../experiments/streets/python_v1/angular_gallery_quality_01/import_quality_checks.json)
passes all 74 networks and 447,444 quads, with exact old reference geometry,
no overlaps and maximum published edge ratio 7.97745. The
[full regression](../../experiments/streets/python_v1/angular_suite_dev52_01/suite_checks.json)
passes **874 tests in 95 modules**, zero failures/errors, on the same frozen
source. The [detached five-job audit](../../experiments/streets/python_v1/dev52_angular_receipt_audit.json)
verifies actual terminal exits, code/script hashes, threads1 and the unchanged
1 GiB cap/reserve. Highest resident/private usage is
476,282,880/540,127,232 bytes. The subsequent standard-tool and fresh real-import
qualification below now pass. Installed-copy qualification also passes for
the exact d8d9 package. The broader completion gates remain pending.

The [dev52 source recovery archive](../../experiments/streets/python_v1/checkpoints/dev52_01/streets_python_dev52_source_checkpoint.zip)
preserves this qualified core/native milestone. Its
[checkpoint record](../../experiments/streets/python_v1/checkpoints/dev52_01/checkpoint_checks.json)
verifies all 1,256 payloads, including 1,193 Python files, and reconstructs the
exact source fingerprint above. The archive is 23,390,714 bytes, SHA-256
`b0f51b2dc488c22b11674325a59481a414a6922d6c4a9980dba302004658b805`.
It is recovery data, not an installed-copy-qualified release.

## Remaining work before goal completion

The completed reproducibility increment moves angular fixture authoring into
`tests/_python_angular_fixtures.py`, shared by the numerical tests and ordinary
`check_python_street_scene.py` generator. The standard scene/reopen checks now
include the twelve angular inputs, real Update-operator edits, furniture poses
and R2/R4 float32 publication. This changes test tooling, not runtime geometry.
Qualification of those entry points now passes; the receipts below distinguish
the exact tooling source from the preserved runtime checkpoint.
The [tooling source audit](../../experiments/streets/python_v1/dev52_tooling_source_audit_02.json)
verifies all 1,052 runtime modules and every unaffected test file are exact
against that checkpoint. Only the shared fixture helper and three related test
entry points differ. The [focused check](../../experiments/streets/python_v1/angular_fixture_tools_01/fixture_tool_checks.json)
passes all seven affected tests and compares the twelve reference inputs and
four additional authored variants exactly with the preserved fixture factory.
The first native worker passed all angular edits and R2/R4, but its global
cache clearing invalidated an older test's warm-baseline assumption. The
corrected generator completes the older dependency checks before the angular
cold-cache checks. The [corrected native worker](../../experiments/streets/python_v1/angular_standard_native_02/scene_checks.json)
now passes all 74 fixtures, all twelve button-edit sequences and actual R2/E0
and R4/E10 publication. Its peak resident/private usage is
523,468,800/594,665,472 bytes. The
[independent saved quality check](../../experiments/streets/python_v1/angular_standard_quality_02/import_quality_checks.json)
passes all 447,444 quads and exact reference positions/connectivity. The
[standard fresh reopening](../../experiments/streets/python_v1/angular_standard_reopen_02/reopen_checks.json)
also passes all 74 controllers, twelve angular button/pose sequences and R2/R4.
All [three actual jobs](../../experiments/streets/python_v1/angular_fixture_tool_qualification_02.json)
are terminal success on source d8d9f2c4...;
the failed worker is retained and is not counted as a passing scene gate.
All [nine fresh real-input jobs](../../experiments/streets/python_v1/dev52_real_import_checks_01.json)
now pass on the same source: complete imports of Berlin (398,972 quads), Canal
(21,818) and Nashville (263,213), independent exact saved-buffer comparisons
against dev51, and fresh seam-edit/local-full/restoration checks. Each import
has both components READY and zero GN modifiers. The largest worker uses
657,723,392 resident / 703,758,336 private bytes, below the same 1 GiB cap with
the same reserve. The measured seam Updates take 1.62/29.04/3.23 seconds;
Canal remains expensive. These are individual timings, not a controlled
performance comparison. All four installed-copy jobs and their detached audit
also pass as recorded at the beginning of this document. The installed Canal
edit takes 8.64 seconds; variation from the separate 29.04 s source run means
these individual measurements establish no new speedup claim.

The historical [12-case curved Bundle diagnostic](../../experiments/streets/python_v1/nonmonotone_bundle_probe_01/nonmonotone_probe.json)
makes the nonmonotone-union gap concrete. All individual Exact/Smooth offset
paths and continuous centerline separation proofs pass. Constant 4 m gaps
between two 6 m carriageways require a pavement union with no separator;
4-to-8-to-4 m gaps require a separator island within that union. The dev51
180/225/270-degree inputs all fall through to separated construction and fail
there. A successful source validity check is not support for the resulting
surface. The diagnostic changes no production source and peaks at 40,615,936
resident bytes. Dev52 now implements and qualifies these 12 free-end examples;
their failure is a dev51 record, not the current runtime status.

The separate [angular preparation](../../experiments/streets/python_v1/polar_bundle_frame_01/polar_frame_checks.json)
now passes all 12 source inputs, with correct union/hole topology and original
end cuts retained. Its worker peaks at 54,374,400 resident bytes. The
[six constant-overlap planar quad fills](../../experiments/streets/python_v1/polar_bundle_quads_03/polar_quad_checks.json)
also pass positive convex quads, no world overlaps and ratio below 6, using
863–2,098 quads and 44,126,208 peak resident bytes. These are preparation/planar
experiments and do not count as runtime support. Subsequent
[finite-cap profiled construction](../../experiments/streets/python_v1/polar_crown_patches_06/polar_crown_patch_checks.json)
passes all 12 examples with 1,985–3,772 direct quads per Bundle, including
180–269 grass quads for each of the six islands. World ratio, physical slopes,
an angular longitudinal-direction bound, ceiling, zero outer curbs, drainage
and complete shared pavement/grass edge partitions pass. The worker peaks at
62,840,832 resident / 81,473,536 private bytes under the same 1 GiB cap/reserve.
The independent
[island boundary audit](../../experiments/streets/python_v1/polar_island_boundary_audit_01/island_boundary_audit.json)
and [surface audit](../../experiments/streets/python_v1/polar_crown_surface_audit_02/crown_surface_audit.json)
pass whole-curve accuracy, cap-corner retention and simulated float32 publication
buffers at exaggerations 1/10. Those earlier files remain external prototypes.
The dev52 runtime and native checks above now qualify provider integration,
additional authored dimensions/orientations, Blender publication, Update and
furniture dependencies for the initial branch. Dev53 now adds compatible
multi-Section angular routes as qualified above. Angular dimensional seams,
actual attached Street ports and more general island topology remain open.
The [drawing and formulas](STREETS_PYTHON_GEOMETRY.md#angular-charts-for-curved-bundle-unions-experimental)
distinguish this evidence from addon support.

1. **Reduce costly Updates.** Dev49 reduces the paired Canal median by about
   20%, with exact output geometry and the same quality checks. The run has
   substantial absolute timing variation; it does not establish stable latency.
   Preserving unaffected outputs already works, but reconstructing an affected
   complex intersection is still expensive. Fixed-rim quad quality descent is
   a remaining profiling target. Further changes must preserve the geometry
   and demonstrate a gain on identical inputs.
   Subsequent native profiling confirms fixed-rim quad quality descent dominates
   the expensive Canal edit. The
   [first balanced candidate comparison](../../experiments/streets/python_v1/canal_quad_reuse_compare_02/quad_reuse_comparison.json)
   preserves exact output and shows only a modest experimental gain. A
   [second comparison](../../experiments/streets/python_v1/canal_edge_scalar_compare_01/quad_reuse_comparison.json)
   also preserves edited/local-full/restored buffers, but its overall candidate
   median is worse (9.83 s versus 8.36 s), with substantial timing drift during
   the run. Neither candidate has been adopted into the addon; these tests do
   not establish that responsiveness is solved.
   Dev54's separately qualified integrated invalid-corner shortcut now lowers
   the paired whole-Update median from 8.5162 s to 7.4205 s (12.87%), preserving
   exact output and the original guards. All 885 tests and the 86-fixture
   native/quality checks pass. This reduces redundant rejected-trial work;
   complex Updates still take seconds and there is no general latency bound.
2. **Complete the outstanding valid-input construction branches.** These include
   nonmonotone overlapping Bundles beyond dev53's compatible multi-Section
   free-end angular routes: dimensional changes, bent Exact Section joins,
   external Street ports and more general separator islands/internal
   junctions; further profile-family changes absorbed into junctions (dev45
   implements the same-family crown-rise case); dimensional seams
   claimed by multiple junctions beyond dev46's complete full-profile exterior
   construction (dev48 also covers complete compatible internal multi-Section
   routes). The 60 m, 150-degree Smooth seam is already qualified by the
   dev50/dev51 work described below. Other sharp Smooth configurations and
   bent mixed Exact/Smooth joins remain open; `SectionTransition.solve`
   explicitly rejects a bent mixed-mode join before construction.
   The qualified examples cover specific constructions, not every configuration
   in these families. Clear rejection messages do not count as support for valid
   configurations.
3. **Qualify and package the final implementation.** After those code changes,
   verify the combined scene, independent real OSM imports, incremental/full
   agreement, reopening, mesh/profile/slope quality and the existing 1 GiB
   process-memory limit. Refresh the installable package and documentation from
   that exact source. The qualified dev53 preview remains available while this
   work continues; older immutable packages are kept for comparison.

### Dev51 common boundary sampling: qualified core and real networks

The [before audit](../../experiments/streets/python_v1/smooth_rim_contract_before_01/rim_contract_checks.json)
finds that the Bundle helper reconstructs the new roads with the legacy global
sampling policy: 215/377 samples versus the actual 149/207. All four endpoints
agree but all four sampled edge chains differ. The existing core suite did not
cover this reconstruction contract.

Dev51 introduces immutable `RoadSampling` for the retained source interval and
its local-curvature policy. Pavement and Bundle boundaries use the same class;
the policy flag survives in road metadata after numerical mesh buffers are
released. Existing nonfallback roads retain their previous metadata and rows.
The [after audit](../../experiments/streets/python_v1/smooth_rim_contract_after_01/rim_contract_checks.json)
matches all four chains exactly, and the new regression covers both authoring
directions and geometry-buffer release. The
[focused test run](../../experiments/streets/python_v1/smooth_rim_core_after_01/smooth_core_checks.json)
passes all 30 tests on source
`8ca473ba7dd1453b01c8043849475f2f1201f42a9c76ad12da022c1c810d2615`.
The [62-fixture upgrade](../../experiments/streets/python_v1/gallery_upgrade_dev51_01/upgrade_checks.json)
preserves all source records, mesh buffers and semantic keys. Its
[saved quality audit](../../experiments/streets/python_v1/gallery_quality_dev51_01/import_quality_checks.json)
passes all 369,835 quads and the unchanged ratio-8/crown/overlap guards.
[Fresh editing/reopening](../../experiments/streets/python_v1/reopen_56/reopen_checks.json)
also passes all 62 controllers. The
[complete regression](../../experiments/streets/python_v1/suite_57/suite_checks.json)
passes **863 tests in 93 modules**, zero failures/errors, on the identical source.
Its worker uses 159,162,368 resident / 180,002,816 private bytes. All three
real-network upgrades now pass with exact source/mesh/semantic keys: Berlin
398,972 quads, Canal 21,818 and Nashville 263,213. The
[terminal receipt audit](../../experiments/streets/python_v1/dev51_receipt_audit.json)
checks all nine actual script hashes, source fingerprint, terminal exits and
memory caps. The largest worker uses 609,128,448 resident / 674,975,744 private
bytes, below 1 GiB. These checks do not prove a complete new sharp Bundle family.

The [dev51 source recovery archive](../../experiments/streets/python_v1/checkpoints/dev51_01/streets_python_dev51_source_checkpoint.zip)
preserves all 1,186 Python files, required libraries, current documents/drawings,
the qualified 62-fixture scene and selected evidence, including the explicitly
experimental angular charts. The [checkpoint record](../../experiments/streets/python_v1/checkpoints/dev51_01/checkpoint_checks.json)
verifies every payload and reconstructs the exact 8ca473ba... source fingerprint.
The archive is 22,026,952 bytes, SHA-256
`ffc7d2abebc1bb8dadbcb85cc3c3081d3abe6f7efa47e5011d5ce3a2f0705b31`.
It is recovery data, not an installed-copy-qualified or published release.

The new [drawing and formulas](STREETS_PYTHON_GEOMETRY.md#local-curvature-certificate-for-smooth-seams-dev50)
explain the unchanged curve, tighter certificate and identical shared rows.
The drawing export reads the actual dev51 saved Blender buffers and verifies
they stay unchanged after reconstruction; it is separate from the independent
quality audit. The corresponding installed-copy checkpoint was dev49; the
current dev52 installer is qualified separately above.

### Dev50 local curvature proof: qualified geometry core

The new `sharp_smooth_seam_60` fixture is the 60 m, 150-degree dimensional
seam in Smooth mode. The legacy global curvature estimate rejects a valid
quintic. A bounded dyadic Bernstein proof certifies its unchanged polynomial
against the same 2 m minimum radius. This retry runs only after all primary
Smooth candidates fail, preserving prior accepted candidate selection.

The [before check](../../experiments/streets/python_v1/smooth_core_before_01/smooth_core_checks.json)
records eight failing new subcases with the old production implementation.
The [focused production check](../../experiments/streets/python_v1/smooth_core_after_04/smooth_core_checks.json)
passes all 29 tests across the continuation, Smooth seam, Exact seam and fitted
Exact modules. The first after-run is retained: one new transform assertion
used less than four absolute-coordinate ulps at kilometre coordinates; another
test tried to compare deliberately released remote numerical mesh buffers.
The corrected tests budget coordinate roundoff and combine the verified retained
buffer with rebuilt parts, without changing geometric guards or production code.

The initial 62-fixture generation and Update checks pass, but the
[saved quality audit](../../experiments/streets/python_v1/gallery_quality_55/import_quality_checks.json)
rejects one new `sharp_smooth_seam_60 / road:2` cell at ratio 8.05471. All 342
preceding mesh parts stay exactly unchanged; all quads are convex, have no
published overlaps and stay within their crown envelope. This is a real new
grading failure, not an accepted tolerance exception.

The [visible-interval experiment](../../experiments/streets/python_v1/visible_curvature_probe_01/visible_curvature_probe.json)
shows the source-wide curvature proof includes the portion already consumed
by the seam. Restricting the Bernstein proof to the actual visible source span
reduces its bound from 0.148174 to 0.105698 per metre. The same 5 mm sampling
budget then needs fewer road rows and improves the new road's edge ratio to
4.39541. Production applies this only beside a new local-proof Smooth fallback;
existing qualified roads keep their row layout. The additional focused test
checks both sampled curbs against 12,001 dense source queries per road.

The corrected [62-fixture scene](../../experiments/streets/python_v1/native_57/python_streets_checkpoint.blend)
passes native generation, edits and restoration with **369,835 direct quads**,
zero GN modifiers. Its worker peaks at 424,243,200 resident / 486,817,792 private
bytes. The new complete fixture uses 7,552 quads, including its control road;
the transition itself still uses 960. The
[independent saved quality comparison](../../experiments/streets/python_v1/gallery_quality_56/import_quality_checks.json)
passes all 62 networks, preserves all 342 previous mesh parts exactly, and has
zero cells above ratio 8 (gallery maximum 7.97745), no overlaps and bounded crowns.

[Fresh reopening](../../experiments/streets/python_v1/reopen_55/reopen_checks.json)
passes all 62 controllers, actual edits, local/full agreement, restoration,
attributes, separators, furniture and reconstruction. Its worker peaks at
425,246,720 resident / 487,698,432 private bytes. The
[independent saved-buffer audit](../../experiments/streets/python_v1/gallery_saved_dev50_01/compound_saved_mesh_checks.json)
passes all six close-junction surfaces, whole entrance slopes, zero curbs,
crown limits and graph-depression checks. All 62 fixture world boxes stay
disjoint by at least 10 m. No addon geometry builder is imported by that audit.

All six jobs in the [dev50 core qualification](../../experiments/streets/python_v1/dev50_checks_v3.json)
pass on the same corrected source. The
[complete suite](../../experiments/streets/python_v1/suite_56/suite_checks.json)
passes **862 tests in 93 modules**, no failures/errors, with 160,006,144 resident
and 176,361,472 private bytes. The subsequent shared-rim correction above has
its own source fingerprint and pending whole-release checks. The 30 m Smooth source still fails its offset-curvature proof;
mixed-mode and broader sharp configurations remain open construction work.

### Dev49 measured quad optimizations: qualified core and real networks

The two candidates described below are now integrated. The
[focused production tests](../../experiments/streets/python_v1/quad_hotspot_production_01/production_checks.json)
pass 33 tests on source
`d4ed8df50b4d699aa3a1147ef4391efae799d91be3ca7cafafe026f2d7c3729a`.
New differential tests compare the complete sorted corner options against the
scalar dev48 oracle for different work chunks and transformed rims, and check
binary64 gradient accumulation against the legacy ordered scatter. The existing
independent derivative and complete-descent tests continue to pass. A direct
quantization-sensitive test proves that disabled height projection stays disabled.

The [serial production qualification](../../experiments/streets/python_v1/dev49_checks.json)
has passed exact all-fixture/Canal upgrades and the
[paired whole-Update production comparison](../../experiments/streets/python_v1/canal_quad_hotspot_production_01/comparison_checks.json).
All 61 fixture source records, mesh buffers and semantic keys remain exact;
both Canal components retain all 21,818 quads exactly. The paired production test loads
the original hotspot functions from the
[verified dev48 source checkpoint](../../experiments/streets/python_v1/checkpoints/dev48_01/checkpoint_checks.json),
whose ZIP hash is fixed in the test. It compares those methods with the actual
dev49 production classes; the experiment does not substitute another candidate
implementation for the new production code.

The dev48 checkpoint preserves all 1,184 Python files, street documentation and
required street libraries, with every payload hash checked. It is a source
checkpoint, not an installed-copy-qualified addon release. All five jobs in
the dev49 core qualification now pass. Remaining geometry families and
the final whole-goal release stay open; the real-network qualification below
and this checkpoint's installed-copy qualification have completed.

The actual production methods reduce the paired ABBA median from **22.2428 s
to 17.7977 s**, approximately **20.0%**. Individual reference times are
33.9181/10.5675 s; production times are 26.9353/8.6602 s. The two timing regimes
vary substantially, so these figures do not establish a stable absolute latency.
Both slow and fast pairs show an improvement. Every edited buffer matches the
archived methods exactly; candidate incremental/full results agree, source and
meshes restore exactly, 14/20 objects remain and the other component is untouched.
Peak resident/private memory is 292,962,304 / 357,208,064 bytes.
[Fresh reopening](../../experiments/streets/python_v1/reopen_54/reopen_checks.json)
also passes all 61 controllers. The
[full numerical regression](../../experiments/streets/python_v1/suite_55/suite_checks.json)
passes **855 tests in 92 modules**, with no failures or errors.

#### Real-network publication precision and saved-scene qualification

All seven jobs in the
[real-network and saved-gallery qualification](../../experiments/streets/python_v1/dev49_real_checks_v2.json)
pass on the identical dev49 source. Berlin and Nashville pass whole-scene
quality audits, fresh reopening, actual source edits, incremental/full equality
and restoration. The largest worker peak is 659,161,088 bytes resident and
705,019,904 bytes private, below the unchanged 1 GiB cap. The independent saved
gallery audit passes all six close-junction cases and all 61 fixture boxes.
Canal's exact upgrade, actual Update and reopening are covered by the core queue.

The first strict upgrades correctly rejected changes in Berlin and Nashville;
their failed receipts remain. These are dev48 publication-precision repairs,
not different XY geometry from the dev49 optimization. Disabling only that
guard in a test worker reproduces every old real-network mesh exactly;
restoring it reproduces the repaired buffers exactly. Berlin changes 21
interior height buffers by at most **0.147 mm**, Nashville 24 by at most
**0.101 mm**. Every XY coordinate, quad, complete XYZ boundary, source record,
semantic ID and object transform stays exact. Independent audits require the
old saved buffers to exceed their unchanged slope limits and the repairs to
satisfy them; crown/datums and graph depression checks also pass.

The [Berlin audit](../../experiments/streets/python_v1/berlin_publication_audit_dev49_02/publication_audit_checks.json)
and [Nashville audit](../../experiments/streets/python_v1/nashville_publication_audit_dev49_01/publication_audit_checks.json)
use declared source entrance rows and actual offset velocities. Their maximum
entrance longitudinal slopes are 1.6695% and 1.0611%, below the same 3% limit.
An earlier checker mistakenly classified raised inner separator curbs as mouths;
its failure is retained. The corrected checker uses actual source ports,
without relaxing the entrance limit or treating all positive boundary heights
as entrances. The qualified scenes retain 398,972 and 263,213 direct quads.

#### Dev49 installed development deliverable

The [package record](../../experiments/streets/python_v1/packages/dev49_v1/package_checks.json)
verifies 17,371,497 bytes, all 1,231 payload files and 1,185 Python files.
Archive SHA-256 is
`05525f86466c79ae498166b802b7827f96287b54dfe909b9224433231d2e8bad`.
The immutable ZIP contains its build-time documentation; this detached record
adds completed installation evidence without rewriting the tested archive.
The installed profile is inside the workspace; original Blender preferences
were not saved or changed. No runtime imports come from either development checkout.

The installed Canal scene has 21,818 quads across two components, exactly
matching the development checkpoint. Its actual seam/width edit takes 9.0226 s,
forced full regeneration 9.5514 s and restoration 8.1377 s on this run. All
buffers agree and restore exactly; 14/20 outputs and the other component are
retained. These separate installed timings are not a controlled speed comparison.
The core ABBA experiment above supplies the measured relative optimization.

### Dev48 internal routes: qualified geometry core

Dev48 implements full original-interval coverage for internal
multi-Section routes, conditional compound-only height reserves, and coincident
computational cuts that retain all original source IDs. The
[focused production regression](../../experiments/streets/python_v1/compound_routes_production_03/focused_checks.json)
passes **26 tests** on source
`a0543f7d8a81433136a7223a2dc668b06efb91092cfe214b178cc075c01122e1`.
It covers both interpolation modes, width/rise combinations, authoring
directions, common datum transforms, free terminals, incompatible seams,
transactional failed Updates, and released-buffer local/full/restoration checks.

The initial native scene (`native_54`) passes all 61-fixture button/edit/restore
checks. Its independent saved-buffer audit then finds a float32 slope excess
in `compound_short_prefix`: 0.10000038, against the unchanged 0.10 limit.
The current source adds a conditional publication-rounding correction, keeping
the original rim, XY and direct quads fixed. Its new regression covers
exaggerations 0/1/10 and datums 0/37 m. The initial scene is retained as evidence
of that failure, not advertised as the final corrected checkpoint.

The corrected [61-fixture scene](../../experiments/streets/python_v1/native_55/python_streets_checkpoint.blend)
passes generation, actual button edits, local/full agreement and restoration
checks: **362,283 convex direct quads**, zero GN modifiers. Its native worker
peaks at 414,564,352 resident / 483,827,712 private bytes. The
[independent saved-buffer audit](../../experiments/streets/python_v1/compound_routes_saved_dev48_03/compound_saved_mesh_checks.json)
passes all six close-junction surfaces. All 64 whole entrance edges match;
outer curbs stay zero, the exposed ceiling is 8 cm, every curb chord is within
4.937 mm of the original 2 m circle/tangent construction, and graph depression
fill is zero. All 61 fixture bounding boxes are disjoint by at least 10 m.
The new short-prefix surface has 554 quads and a maximum saved bilinear slope
of 0.09999884, below the unchanged 0.10 limit.

The exact baseline comparison initially rejects a changed `junction_mesh:1`.
The [independent delta audit](../../experiments/streets/python_v1/publication_delta_dev48_01/publication_delta_checks.json)
then establishes that all **323 other previous mesh parts are exact**, all
58 previous source records and furniture poses are unchanged, and the only
repair belongs to `profiled_corridor_crossing`. Its maximum physical Z change
is 0.5484 micrometres; XY, quads and the entire XYZ boundary remain exact. The
old saved slope 0.14080039 exceeded its 0.14080000 limit; the corrected value
is 0.14079785. Only that proved change is accepted in the resumed comparison.

The [continuing serial qualification](../../experiments/streets/python_v1/dev48_checks_v2.json)
reuses the two already completed native/audit jobs on the identical source.
[Broad gallery quality](../../experiments/streets/python_v1/gallery_quality_53/import_quality_checks.json)
passes all 61 networks: positive convex quads, no published overlaps, bounded
edge ratios and the common crown envelope. The paired exact delta audit above
limits the permitted historical difference to that one named mesh and repair.
[Fresh reopening](../../experiments/streets/python_v1/reopen_53/reopen_checks.json)
also passes all 61 controllers, including attribute/centerline edits, exact
restoration, dependencies, separators, furniture and zero GN modifiers.
Its worker peaks at 419,790,848 resident / 486,244,352 private bytes.
The [full numerical suite](../../experiments/streets/python_v1/suite_54/suite_checks.json)
passes **852 tests in 92 modules**, zero failures/errors, on the same frozen
source. All five jobs in the serial core qualification are terminal and pass.
The [Canal version upgrade](../../experiments/streets/python_v1/canal_upgrade_dev48_01/upgrade_checks.json)
also preserves both components' source records, semantic keys and all 21,818
published quads exactly.
Real-network regeneration, resource/performance comparison and installer
qualification remain separate pending gates.
Dev47 remains the latest fully qualified installer while these finish.
The [geometry section](STREETS_PYTHON_GEOMETRY.md#internal-section-routes-in-close-junctions-dev48)
documents the ownership, reserve and partition formulas.

#### Measured optimization candidates, not yet integrated

The [test-only kernel comparison](../../experiments/streets/python_v1/quad_hotspot_candidates_02/hotspot_checks.json)
checks ordered `bincount` gradient accumulation and bounded batched corner-set
enumeration in Blender 5.2 / NumPy 2.3.4. Six saved close-junction meshes and a
4,096-quad regular grid retain exact gradient values, full bounded-descent
outputs, and every sorted corner candidate/score across chunk sizes.

The [controlled whole-Update comparison](../../experiments/streets/python_v1/canal_quad_hotspots_01/comparison_checks.json)
then uses the same Connection 12 move and width edit in Canal Street. Interleaved
reference/candidate/candidate/reference measurements give median wall times
**11.6807 s / 9.0678 s**, a 22.4% reduction on this run. CPU medians are
11.5547 s / 8.8516 s. Every edited buffer is exact across methods; candidate
incremental/full builds agree, all original inputs and buffers restore exactly,
14/20 output objects remain, and the other component stays untouched. Peak
resident/private memory is 289,636,352 / 347,660,288 bytes. The candidate profile
still identifies fixed-rim quality descent as the main remaining cost.

These paired timings support the dev49 integration above. They do not establish
a speedup from the earlier 33.87 s installed qualification, measured under
different conditions. Production source remains unchanged during the experiment;
the production integration has separate focused tests and ongoing qualification.

### Dev47 exposed-profile limits: qualified installable checkpoint

The [nine focused boundary tests](../../experiments/streets/python_v1/boundary_scope_after_02/scope_checks.json)
verify isolated layout records, exposed crown limits and the existing related-
source common-minimum invariant. The previous implementation fails the new
high-internal-crown and unrelated-record regressions. The
[complete suite](../../experiments/streets/python_v1/suite_53/suite_checks.json)
passes all 843 tests in 91 modules. The
[core receipt audit](../../experiments/streets/python_v1/dev47_audit_checks.json)
verifies every terminal exit, script hash, source fingerprint and memory limit.
The [current scene](../../experiments/streets/python_v1/native_53/python_streets_checkpoint.blend)
and [actual-mesh drawing](STREETS_PYTHON_GEOMETRY.md#exposed-profiles-determine-the-height-ceiling-dev47)
show the internal 18 cm source profile beneath an exposed 8 cm physical ceiling.

The [package record](../../experiments/streets/python_v1/packages/dev47_v1/package_checks.json)
verifies 17,228,589 bytes, 1,226 payload files and all 1,181 Python files.
Archive SHA-256:
`7617b685971bbd2f417ad344e5e0b9c211d4cdc03fb3b92be77688bddbe8e8e2`.
The actual Install Add-on operator, complete registration/unregistration, GN
and Python choices, external material library, manual Update, all three new
compound fixtures and exact restoration pass using exclusively installed
runtime modules. Registration peaks at 693,796,864 bytes resident. Normal Canal
import generates both components and 21,818 quads; its 21 saved mesh parts are
exact against the development result. Fresh Connection-12 editing preserves
14/20 output objects and the other component, agrees with a forced full rebuild,
and restores exactly. Local/full/restore timings are 33.87/35.17/31.51 s in
that separate qualification run. Every worker stays below 1 GiB resident/private.
The archive bytes remain unchanged after qualification; detached receipts record
the final result, while its packaged documentation retains build-time status.

The immediate construction experiment addresses a dimensional seam claimed by
two nearby junctions. An external
[complete-port prototype](../../experiments/streets/python_v1/compound_exterior_ports_03/compound_exterior_ports_checks.json)
passes all eight previously failing synthetic cases, seven source-reversal or
transformed-resolution cases, and a released-cache incremental/full/restoration
check. The worker peaks at 73,129,984 bytes resident; production source remains
unchanged. This was evidence for physical entrance planning, not shipped
support. The production follow-up below extends its uniform-radius scope and
qualifies local dependencies and Blender publication. General compound/Bundle
policies remain unfinished. See the
[geometric construction review](STREETS_PYTHON_GEOMETRY.md#next-review-complete-ports-for-a-compound-dimensional-seam).

The production `CompoundPortPlanner` now replaces those test-only patches in
dev46, retains original source/owner identities, and includes bounded whole-
segment centerline support checks and local dependency signatures. Its
[focused 20-test run](../../experiments/streets/python_v1/compound_ports_unit_02/unit_checks.json)
passes heterogeneous strict radii, transformed datums/resolutions, complete
profile edges, zero curbs, furniture, cache release/local-full restoration,
unrelated-compound locality and negative support cases. The
[57-fixture native scene](../../experiments/streets/python_v1/native_51/python_streets_checkpoint.blend)
passes actual Update, parameter and source edits. Its
[saved-buffer comparison](../../experiments/streets/python_v1/gallery_quality_50/import_quality_checks.json)
preserves all 306 previous parts exactly. Fresh reopening also passes for
all 57. The [current-source suite](../../experiments/streets/python_v1/suite_52/suite_checks.json)
passes 840 tests in 91 modules, zero failures/errors. The
[three-job queue](../../experiments/streets/python_v1/dev46_checks.json)
verifies terminal success, entry-script hashes and unchanged source
`3e3485bd61dded8fd92b3d1ae9db55325be7722f5c331cbd43827f78f443c66d`.
The [independent saved-mesh audit](../../experiments/streets/python_v1/compound_ports_saved_03/compound_saved_mesh_checks.json)
imports no runtime geometry. Each new 721-quad intersection has 64 complete
shared profile edges, zero curbs, a maximum physical height below 0.08 m,
complete bilinear-quad slopes below 9.362%, entrance longitudinal slopes below
0.169%, and zero graph depression fill. All complete curb chords are within
4.937 mm of independent analytical circles/tangents under the 5 mm limit.
The first reference used 1,025 samples per chord and its conservative margin
gave a 5.039 mm upper bound; increasing only the test reference to 4,097 samples
tightens that bound to 4.937 mm. Runtime meshes, sampling and tolerance are
unchanged. The failed earlier audit receipts remain historical evidence.
Generation/reopening peaks are 400,310,272/406,745,088 bytes resident, with
sampled private memory below 1 GiB in every worker. This extends a construction
branch; it does not qualify general
half-profile, separator-hole or arbitrary compound incidence policies.

The current [three-repeat optimizer comparison](../../experiments/streets/python_v1/quad_workspace_compare_02/workspace_checks.json)
replays six actual Canal quad-quality calls with alternating implementation
order in one Blender process. The prototype produces bit-identical positions,
topology, iteration counts and energies on every replay. The sum of per-case
median times falls from 26.714 s to 22.149 s, a 17.1% reduction for those
optimizer calls. This is not a measured full-Update speedup. The prototype
was outside the runtime for this experiment; its integration is qualified below.

### Dev44 optimizer integration: qualified development checkpoint

Development Python-source SHA-256:
`1ace9e03256f06763f5adb7cc243ac2d48d06cace1a870c1564fa379f6a097e2`.
The [five maintained focused tests](../../experiments/streets/python_v1/quad_workspace_unit_02/unit_checks.json)
pass independent derivative, trial-buffer isolation, exterior/hole-rim,
height/topology/material and descent-budget checks. The
[native integrated replay](../../experiments/streets/python_v1/quad_workspace_integrated_01/integrated_workspace_checks.json)
matches all six captured calls exactly, including positions, topology, accepted
steps and energies. Neither check changes the pinned source.

The [same-process Canal comparison](../../experiments/streets/python_v1/canal_workspace_update_01/update_comparison_checks.json)
passes two alternating old/new local-Update pairs on Connection 12. Median time
is 11.860 s with the previous optimizer and 10.583 s with dev44: 10.8% less time.
All edited mesh buffers are exact across implementations and repetitions.
Incremental/full equality, exact source/mesh restoration, 14/20 retained output
objects and an unchanged other component pass. The worker peaks at 288,763,904
bytes resident under the unchanged 1 GiB process cap/reserve. The legacy method
is selected only inside this test process, never by shipped runtime code.
These timings do not establish a reduction from the earlier installed 37.51 s
run, whose execution conditions differed. This measured gain is useful but does
not resolve the remaining Update cost.

The [full numerical regression](../../experiments/streets/python_v1/suite_49/suite_checks.json)
passes 832 tests in 89 modules, zero failures/errors and identical source
fingerprints. The [combined scene](../../experiments/streets/python_v1/native_49/python_streets_checkpoint.blend)
passes [generation and Update](../../experiments/streets/python_v1/native_49/scene_checks.json)
and [fresh reopening/editing](../../experiments/streets/python_v1/reopen_49/reopen_checks.json)
for all 53 fixtures. Its [independent saved-buffer audit](../../experiments/streets/python_v1/gallery_quality_48/import_quality_checks.json)
proves all 296 previous mesh parts exact in positions and topology: 326,458
convex quads, no overlaps and maximum edge ratio 7.97745. The
[complete seven-job queue](../../experiments/streets/python_v1/dev44_checks.json)
also passes exact source, mesh-buffer and semantic-key upgrades of Berlin
(398,972 quads), Canal Street (21,818) and Nashville (263,213), with both
components retained in each. The
[receipt audit](../../experiments/streets/python_v1/dev44_audit_checks.json)
verifies every terminal exit status and entry-script hash against the unchanged
1,178-file Python source. The largest worker peaks at 604,377,088 bytes resident
and 645,521,408 bytes sampled private memory, below the 1 GiB process cap.
[Actual saved bootstrap/demo execution](../../experiments/streets/python_v1/scene_demo_dev44_01/demo_checks.json)
also passes both interpolation choices and exact restoration. The
[dev44 ZIP](../../experiments/streets/python_v1/packages/dev44_v1/blosm_streets_python_dev44.zip)
passes [all four installed-copy jobs](../../experiments/streets/python_v1/package_dev44_checks.json).
The actual Install Add-on operator, full registration/unregistration, both
backends, external material library, Update button, normal Canal import,
published-buffer audit and fresh Connection-12 Update all pass. Runtime origins
are exclusively in the isolated installed copy; its source remains unchanged.
The [package record](../../experiments/streets/python_v1/packages/dev44_v1/package_checks.json)
verifies the unchanged 17,016,747-byte archive, SHA-256
`3643367c12b60aba1809761e57810a9f1afe959a562dbd76716a7c1809770bc0`.
Its 1,220 payload files contain all 1,178 qualified Python files. Installed
registration peaks at 660,336,640 bytes resident; all jobs stay below 1 GiB
resident/private memory. The installed Canal local Update takes 34.01 s in its
separate qualification run; use the controlled comparison above for the measured
optimization gain. Dev45 development began after freezing this ZIP; its tests
still use only the unchanged installed dev44 source.
An initial queue-wrapper failure after successful generation is
preserved in `native_49/qualify_dev44_failed.py`: it expected `passed` instead
of the generator's `completed` success field. The corrected queue reused only
terminal successful receipts with matching entry-script hashes; it did not
rerun the completed suite or generation.

### Dev45 crown-rise ownership: qualified development checkpoint

Source fingerprint `4cb19348505fd1f104e897a6cb65484ce0d6d6b59790dc4e8355c3f195a14e93`.
[Three maintained focused tests](../../experiments/streets/python_v1/absorbed_profile_unit_02/unit_checks.json)
and [five previous ownership tests](../../experiments/streets/python_v1/absorbed_legacy_unit_dev45_01/unit_checks.json)
pass. New coverage includes 64 mode/width/rise/reversal combinations, transformed
datum/resolution cases, complete mouth profiles, zero curbs, 2 m circles, crown
ceilings, furniture, released buffers, locality and incremental/full restoration.
The [two Blender fixtures](../../experiments/streets/python_v1/absorbed_profile_native_01/profile_native_checks.json)
pass actual rise Update, ownership reversal, anchor movement, remote-object
retention, exact restoration and saved-scene cold reconstruction. Each complete
fixture has 4,832 quads, including 576 junction quads. Runtime height, slope,
shape and overlap guards are unchanged. The
[combined 55-fixture scene](../../experiments/streets/python_v1/native_50/python_streets_checkpoint.blend)
passes [generation/Update](../../experiments/streets/python_v1/native_50/scene_checks.json),
[saved-buffer quality](../../experiments/streets/python_v1/gallery_quality_49/import_quality_checks.json)
and [fresh reopening](../../experiments/streets/python_v1/reopen_50/reopen_checks.json).
It has 336,122 convex quads, no overlaps, maximum ratio 7.97745 and all 296
previous parts exact. The fresh checks include both new fixtures' actual rise
Update and ownership reversal. All three real-network upgrades and saved-text
checks also pass on the unchanged `4cb19348…` source. The
[complete numerical rerun](../../experiments/streets/python_v1/suite_51/suite_checks.json)
passes 835 tests in 90 modules, zero failures/errors and identical source hashes.
The [core receipt audit](../../experiments/streets/python_v1/dev45_audit_checks.json)
verifies all eight terminal jobs and 1,179 Python files. The largest core worker
peaks at 606,846,976 bytes resident and 670,306,304 bytes sampled private memory.
The [independent supplement](../../experiments/streets/python_v1/dev45_supplement_checks.json)
also passes full bilinear-quad slope bounds, complete shared mouth edges,
disjoint world XY bounds for all 55 fixtures and the drawing from this combined
scene. The minimum fixture-box distance remains 10 m. Both new examples share
48 exact float32 mouth edges and have 48 remaining zero-height exterior edges.
The matching [dev45 ZIP](../../experiments/streets/python_v1/packages/dev45_v2/blosm_streets_python_dev45.zip)
passes [all four installed-copy jobs](../../experiments/streets/python_v1/package_dev45_checks.json).
The actual installation, complete addon registration, Python/GN choices,
external material assets, Update button, normal Canal import, saved-buffer
quality and fresh Connection-12 Update all pass. Installed runtime origins and
the final installed Python fingerprint remain exact. Registration peaks at
657,879,040 bytes resident; every worker remains below the 1 GiB limits.
The installed local edit takes 36.13 s, forced full 37.85 s and restoration
33.54 s, retaining 14/20 output objects and the unchanged other component.
These separate qualification timings are not a controlled comparison with
earlier versions. The [package record](../../experiments/streets/python_v1/packages/dev45_v2/package_checks.json)
verifies 17,117,902 bytes, SHA-256
`18c0889ea1486aeb3a988ade407ebdda836e55e906bbd10a165cd6842b47a184`,
1,222 payload files and all 1,179 Python files at the qualified source hash.
ZIP content is unchanged after testing; detached receipts give its completed
qualification. The first dev45 candidate retained the old quickstart scene
reference and was superseded before installation testing by this v2 package.

The initial [835-test run](../../experiments/streets/python_v1/suite_50/suite_checks.json)
has one failure: an old route test expects the newly supported rise change to
raise. The rejection test now checks a genuinely incompatible datum; the
[four focused route checks](../../experiments/streets/python_v1/dev45_routes_unit_01/unit_checks.json)
pass. The [source-preimage proof](../../experiments/streets/python_v1/dev45_route_test_source_proof.json)
reconstructs the preceding `c6f26095…` source exactly with that single test file
replaced. Runtime and all other test files are unchanged. A fresh qualification
queue ran the Blender jobs before the full numerical rerun (`suite_51`). That
complete rerun now passes; the failed predecessor remains as historical evidence.

The [independent saved-buffer audit](../../experiments/streets/python_v1/absorbed_profile_drawing_02/absorbed_profile_meshes.json)
passes both native examples: all 16 complete edges at each of the six mouths
match exactly, exterior curbs have zero physical error, and complete curb
chords are within 3.682 mm of the intended circular/tangent boundary under the
5 mm limit. The [actual-mesh drawing](STREETS_PYTHON_GEOMETRY.md#crown-rise-changes-owned-by-the-junction-dev45)
is generated from that export and has been visually inspected.

A separate [eight-case compound-ownership probe](../../experiments/streets/python_v1/compound_seam_candidate_01/compound_seam_candidate_checks.json)
remains unsuccessful. Both baseline and test-only candidate stop at partial
incoming profiles retained by compound clipping. Assigning an existing common
owner alone is insufficient; full physical entrance geometry must be reconciled
before filling that compound. The candidate is not integrated into the addon.

## Dev43 development increment: bounded sharp Exact seam fitting

The qualified dev43 Python-source fingerprint is
`709ff27409ac4e4be9bb3c0cf860f1d33b4b24de597f76a743f8fb2025c8f4dd`.
[`transition_fitting.py`](../way/street_network/python_backend/transition_fitting.py)
adds a bounded fallback after the established transition attempts fail. It
certifies curb chord error before trying fewer longitudinal rows, and may use
more of an explicitly free source end. It preserves mandatory bend/seam rows,
full mouth profiles, the exact centerline, crown/slope limits and direct quads.
Occupied far junctions, selection cuts and Bundle members do not receive this
free-end expansion. Previously accepted constructions retain their old branch.

All [20 signed-angle cases](../../experiments/streets/python_v1/exact_angles_02/angle_checks.json)
now pass, including the four previously rejected 150-degree cases at source
lengths 30/60 m. The final implementation also passes
[nine maintained focused tests](../../experiments/streets/python_v1/fitted_transitions_unit_03/unit_checks.json),
including redundant collinear source controls, allocation/error budgets,
resolutions 1/2/4, reversed authoring and Update/locality/restoration. That focused
test run preceded only the recorded sampling-piece metadata and native fixture
additions. The [complete current regression](../../experiments/streets/python_v1/suite_48/suite_checks.json)
passes **827 tests in 88 modules**, with zero failures/errors and identical
before/after source hashes. Test execution takes 841.31 s; the guarded worker
peaks at 137,564,160 bytes resident and 156,958,720 bytes sampled private memory.
These timings are not a controlled comparison with the previous run.

[Focused Blender checks on the current source](../../experiments/streets/python_v1/fitted_exact_native_01/fitted_native_checks.json)
pass for both new 150-degree fixtures. Each transition has 768 quads and all
16 complete profile edges at both mouths. Independently reconstructed analytic
curbs differ from published float32 boundaries by at most 3.760 mm, below the
unchanged 5 mm budget. Centerline error is below 0.94 micrometres. Physical
slopes, zero exterior edges, the highest-crown ceiling, real vertex/width edits,
furniture movement, local/full equality, exact restoration and cold-cache
identity all pass. The worker peaks at 238,751,744 bytes resident.

The 30 m example consumes 28.5 m of each free Section, leaving 1.5 m with its
unchanged full profile. The 60 m example consumes 32 m per side. These lengths
are explicit in the metadata and new gallery inputs; the source polylines and
their attributes are not rewritten. The fallback does not solve sharp Smooth
seams: the short Smooth probe still fails offset regularity and the longer one
still lacks the required transition reach. Both new fixtures deliberately use
Exact mode; existing fixtures and the global default retain Smooth.
See the [construction and error bound](STREETS_PYTHON_GEOMETRY.md#bounded-sharp-exact-seam-fitting).

The [combined 53-fixture scene](../../experiments/streets/python_v1/native_48/python_streets_checkpoint.blend)
now passes all seven gallery jobs:

- [generation, editing and restoration](../../experiments/streets/python_v1/native_48/scene_checks.json);
- [published mesh quality and exact old-buffer comparison](../../experiments/streets/python_v1/gallery_quality_47/import_quality_checks.json);
- [fresh reopening and editing of all 53 controllers](../../experiments/streets/python_v1/reopen_48/reopen_checks.json);
- [independent complete-quad physical slopes](../../experiments/streets/python_v1/slopes_gallery_dev43_01/slope_checks.json);
- [whole-edge joins and disjoint fixture bounds](../../experiments/streets/python_v1/gallery_interfaces_10/interface_checks.json);
- [authored circles, profiles and two-ended Bundle topology](../../experiments/streets/python_v1/extended_gallery_audit_03/extended_gallery_checks.json);
- [actual saved startup and demonstration scripts](../../experiments/streets/python_v1/scene_demo_dev43_01/demo_checks.json).

The gallery contains 326,458 convex direct quads, maximum edge ratio 7.97745,
zero GN modifiers and no published mesh overlaps. All 288 preceding mesh parts
retain exact float32 positions and topology. World XY fixture bounds are
disjoint, with a 10 m minimum gap. Generation/Update takes 353.62 s, peaking at
394,358,784 bytes resident; reopening/edit checks take 284.04 s, peaking at
404,643,840 bytes. These are qualification runs, not Update latency benchmarks.
All workers retain the 1 GiB resident/private process cap and 1 GiB system reserve.
The existing 51 fixtures retain Smooth; only the two named sharp regressions
are saved in Exact mode. The global authoring default stays Smooth.

The [completed eleven-job queue](../../experiments/streets/python_v1/dev43_checks.json)
retains the same source hash throughout. All three real-network upgrades keep
both prepared components, source data, mesh buffers and semantic keys exactly:

| Dev43 scene / upgrade receipt | Whole-scene quads | Worker time | Peak resident |
| --- | ---: | ---: | ---: |
| [Berlin](../../experiments/streets/python_v1/berlin_upgrade_dev43_01/python_street_import_upgraded.blend) / [check](../../experiments/streets/python_v1/berlin_upgrade_dev43_01/upgrade_checks.json) | 398,972 | 253.65 s | 599,498,752 bytes |
| [Canal Street](../../experiments/streets/python_v1/canal_upgrade_dev43_01/python_street_import_upgraded.blend) / [check](../../experiments/streets/python_v1/canal_upgrade_dev43_01/upgrade_checks.json) | 21,818 | 38.67 s | 268,746,752 bytes |
| [Nashville](../../experiments/streets/python_v1/nashville_upgrade_dev43_01/python_street_import_upgraded.blend) / [check](../../experiments/streets/python_v1/nashville_upgrade_dev43_01/upgrade_checks.json) | 263,213 | 71.35 s | 458,063,872 bytes |

These are cold version upgrades, not new local-edit latency measurements or
controlled speed comparisons. The final queue job exports the two new sharp
examples from the saved Blender buffers without changing the scene. Their
[topology and physical-height drawing](STREETS_PYTHON_GEOMETRY.md#bounded-sharp-exact-seam-fitting)
uses that export. The isolated ZIP checks below also pass; the goal remains active.

### Installable experimental checkpoint

The [local addon ZIP](../../experiments/streets/python_v1/packages/dev43_v1/blosm_streets_python_dev43.zip)
contains the complete addon source and required street libraries: 17,009,109
bytes, SHA-256 `662d0cc8434cbdd9d8dc60ca98e4eab17630fb28effd8bc21e9675717935485f`.
Its [build and qualification record](../../experiments/streets/python_v1/packages/dev43_v1/package_checks.json)
verifies every payload hash and reconstructs the same 1,176-file Python source
fingerprint from the ZIP. The original build receipt is retained separately;
qualification changes no ZIP bytes. This is a local experimental package, not
a public release or a completed goal.

All [four installed-package checks](../../experiments/streets/python_v1/package_dev43_checks.json)
pass in a separate Blender profile. The actual Install Add-on operator extracts
the ZIP; the complete addon registers and unregisters. Python and GN both
generate the same prepared T into 2,304 quads. The actual Python Update button
changes the dependent mesh, agrees with a full build and restores it exactly.
An external material-only library also works. All imported addon modules come
from the installed package, with no development-checkout imports, and its source
fingerprint remains unchanged.

The [installed normal Canal import](../../experiments/streets/python_v1/installed_import_dev43_02/installed_import_checks.json)
completes both components with 21,818 quads and no GN modifiers. Its
[independent saved-buffer check](../../experiments/streets/python_v1/installed_quality_dev43_02/import_quality_checks.json)
matches all 21 reference mesh parts exactly and passes convexity, edge ratio,
overlap and crown-envelope checks. The
[fresh installed-copy Update](../../experiments/streets/python_v1/installed_reopen_dev43_02/import_reopen_checks.json)
passes Connection 12's coordinated position/width edit, incremental/full equality,
exact restoration, 14/20 retained output objects and an unchanged other component.
Observed cold reconstruction is 34.67 s, local Update 37.51 s, forced full 38.92 s
and restoration 35.00 s. These measurements do not establish a speedup over
older runs: the local operation still costs almost a full rebuild. Profiling
that expense is the next performance task.

The first installation probe is retained as a
[harness failure](../../experiments/streets/python_v1/installed_addon_dev43_01/check_installed_python_addon_failed.py):
it treated Blender's virtual `[bl_ext]` module marker as a physical source path.
The corrected probe passes. A final
[installation-only repeat](../../experiments/streets/python_v1/installed_addon_dev43_03/installation_checks.json)
redirects all Blender user resources to the temporary profile and disables scene
thumbnail generation. It passes without the earlier out-of-profile thumbnail
write attempt. No addon code or ZIP content changes for these harness fixes.
This clean full-registration/backend/Update check peaks at 659,013,632 bytes
resident and 720,252,928 bytes sampled private memory; the installed real-import
and edit workers peak below 288 MB resident. All keep the 1 GiB cap/reserve.
The original Blender preferences and `buildings_gn` checkout remain unchanged.

## Qualified dev42 increment: Exact bent dimensional seams

Qualified dev42 Python-source SHA-256:
`5701c49b6aa6aba03a5fa70a25e97b7e7b9a2fa99ff9957d8c5cdddca979b6f0`.
[Seven focused tests](../../experiments/streets/python_v1/exact_transitions_unit_02/unit_checks.json)
pass on unchanged source. They cover authored controls inside the trim, mandatory
piece knots, complete profile rows, continuous sampling bounds, reverse source
directions, rotation/translation, offsets, unequal rises, datum 37, resolutions
1/2/4, zero display exaggeration, local/full equality and released-cache recovery.

[Focused native checks](../../experiments/streets/python_v1/exact_native_focus_01/exact_native_checks.json)
also pass for an ordinary dimensional seam and the originally failing bent-seam
two-ended Bundle. Actual published float32 coordinates stay within `4.88e-8` m
of the Exact centerline and share all 16 complete profile edges at each mouth.
Update/full equality, exact restoration and cold-cache identity pass. The
36.53 s worker peaks at 245,669,888 bytes resident under the unchanged 1 GiB
process cap and system reserve. These are qualification measurements, not
controlled performance comparisons.

The [complete numerical suite](../../experiments/streets/python_v1/suite_47/suite_checks.json)
passes **818 tests across 87 modules** with zero failures/errors and identical
`9038014b…` before/after fingerprints. Test execution takes 932.04 s; the
guarded worker peaks at 139,694,080 bytes resident and 156,229,632 bytes sampled
private memory. The subsequent
[source preimage proof](../../experiments/streets/python_v1/dev42_scene_instructions_proof.json)
reconstructs that exact source from the intermediate `29036626…` files by reverting only the
generated README string in `check_python_street_scene.py`. AST comparison proves
the Python logic unchanged; production files and numerical tests are identical.
This updates the inspection instructions without repeating unchanged numerical
tests. The nine-job `qualify_dev42.py` native queue ran on the intermediate
`29036626…` fingerprint. All nine jobs in its [completed queue](../../experiments/streets/python_v1/dev42_checks.json)
pass. The gallery checks are:

- [51-fixture scene and Update](../../experiments/streets/python_v1/native_46/scene_checks.json),
  with 398,585,856 bytes peak resident memory;
- [saved-buffer quality and comparison](../../experiments/streets/python_v1/gallery_quality_46/import_quality_checks.json):
  all 321,626 quads pass, with all 281 previous mesh parts exact;
- [fresh reopening and editing](../../experiments/streets/python_v1/reopen_47/reopen_checks.json)
  of all 51 controllers, including both new fixtures in Exact mode;
- [independent physical slopes](../../experiments/streets/python_v1/slopes_gallery_dev42_01/slope_checks.json);
- [whole-edge interfaces and disjoint fixture bounds](../../experiments/streets/python_v1/gallery_interfaces_09/interface_checks.json);
- [independent authored-circle, profile and annulus audit](../../experiments/streets/python_v1/extended_gallery_audit_02/extended_gallery_checks.json).

Open the [current combined scene](../../experiments/streets/python_v1/native_47/python_streets_checkpoint.blend).
All source polylines retain Smooth as the saved default. Its Exact demonstration
text switches either new fixture and presses Update; the quickstart documents
the steps. The previous gallery remains available in the history below.

An actual-text check then exposed a missing `bpy` import in the saved optional
bootstrap (`dev42_scene_demo_01` fails). The
[bootstrap source proof](../../experiments/streets/python_v1/dev42_bootstrap_source_proof.json)
reconstructs the nine-job source exactly by reverting that one generated string;
all production files, numerical tests and verifier logic are unchanged. The
[corrected checkpoint 47](../../experiments/streets/python_v1/native_47/scene_text_checks.json)
is opened again and matches the previous scene's source, mesh, manifest and
transform fingerprint exactly: 51 controllers / 288 meshes / 321,626 quads.
The [actual saved bootstrap/demo check](../../experiments/streets/python_v1/dev42_scene_demo_02/demo_checks.json)
now passes on the current fingerprint for both fixtures in Exact and Smooth,
with exact restoration. This closes the user-facing startup defect without
claiming the entire suite was rerun for a text-only correction.

All three existing real-network scenes upgrade from dev41 with exact authored
source, mesh buffers and semantic keys. Both prepared components survive in
each file. These timings are cold version upgrades under the measured load,
not local-edit latency or a controlled comparison with dev41:

| Dev42 scene / exact upgrade | Quads | Worker time | Peak resident |
| --- | ---: | ---: | ---: |
| [Berlin](../../experiments/streets/python_v1/berlin_upgrade_dev42_01/python_street_import_upgraded.blend) / [check](../../experiments/streets/python_v1/berlin_upgrade_dev42_01/upgrade_checks.json) | 398,972 | 259.69 s | 593,608,704 bytes |
| [Canal Street](../../experiments/streets/python_v1/canal_upgrade_dev42_01/python_street_import_upgraded.blend) / [check](../../experiments/streets/python_v1/canal_upgrade_dev42_01/upgrade_checks.json) | 21,818 | 40.52 s | 276,873,216 bytes |
| [Nashville](../../experiments/streets/python_v1/nashville_upgrade_dev42_01/python_street_import_upgraded.blend) / [check](../../experiments/streets/python_v1/nashville_upgrade_dev42_01/upgrade_checks.json) | 263,213 | 68.49 s | 455,241,728 bytes |

The gallery generation/Update worker takes 344.20 s and fresh reopening/editing
takes 286.50 s. Both peak below 399 MB resident. All nine workers retain the
same 1 GiB resident/private cap and 1 GiB system reserve. The dev41 local-edit
timings below remain versioned evidence, not claimed new dev42 measurements.

The [20-case angle sweep](../../experiments/streets/python_v1/exact_angles_01/angle_checks.json)
passes both signs of 5/35/65/110 degrees for 30 m and 60 m source spans. The
four 150-degree cases still fail orientation or quad-quality construction.
They are retained valid-input failures. Mixed Exact/Smooth bent seams also
remain unfinished. A [24-candidate direct-atlas trial](../../experiments/streets/python_v1/exact_atlas_probe_01/atlas_probe_checks.json)
also fails to resolve them: transverse 1:3 refinement alone meets neither all
quad-quality nor all physical-slope guards. Its code is not in production.
The following dev41 checkpoint is historical.

## Qualified 49-fixture checkpoint (dev41)

The qualified dev41 all-Python source fingerprint is
`689cd567c9ba6d7c8475d3425974fa1abf17863f33b243bd7abae7a32796ff90`.
The [source proof](../../experiments/streets/python_v1/extended_fixture_source_proof.json)
compares it with the verified recovery archive: every production file and
existing numerical test is byte-identical. The fixture module retains its
entire old byte prefix and appends only two new function definitions. The
[three new tests](../../experiments/streets/python_v1/extended_fixtures_unit_01/unit_checks.json)
pass in 14.63 s at 61,300,736 bytes peak resident; only the native verifier files
changed afterward. This is an 808-test baseline plus three focused tests, not
a claimed combined 811-test run.

The [six-job native queue](../../experiments/streets/python_v1/extended_gallery_checks.json)
passes on the unchanged current fingerprint:
[combined scene](../../experiments/streets/python_v1/native_45/python_streets_checkpoint.blend),
[generation/Update](../../experiments/streets/python_v1/native_45/scene_checks.json),
[saved buffers](../../experiments/streets/python_v1/gallery_quality_45/import_quality_checks.json),
[fresh reopening](../../experiments/streets/python_v1/reopen_46/reopen_checks.json),
[physical slopes](../../experiments/streets/python_v1/slopes_gallery_extended_01/slope_checks.json),
[interfaces](../../experiments/streets/python_v1/gallery_interfaces_08/interface_checks.json)
and [independent new-case audit](../../experiments/streets/python_v1/extended_gallery_audit_01/extended_gallery_checks.json).
All 263 prior mesh parts remain exact in position and topology. The new cases
pass position/width edits, local/full equality, exact restoration, furniture
updates and cold-cache identity. Native generation peaks at 393,011,200 bytes
resident under the same 1 GiB process cap/reserve policy.

The independent audit derives circles and full crown heights from authored
inputs, without rebuilding the mesh. It finds 96 / 112 whole mouth edges for
six/seven arms, with at most 4.920 mm curve approximation error. The two-ended
Bundle is one connected annulus with two boundary loops and 32 whole internal
body edges at each cap. Its four exterior terminal corners retain zero height.
See the [actual saved-mesh drawing](STREETS_PYTHON_GEOMETRY.md#higher-valence-and-two-ended-bundle-examples).

## Core qualification history

The preceding fully qualified dev40 source fingerprint is
`cf6ec7323707da1f019e1583062b9db8a39723660ddcfdb9126717a603675ade`.
The [805-test numerical run](../../experiments/streets/python_v1/suite_45/suite_checks.json)
passes on `cd96f6009a7efba9f0f02232eae180bd283f5b4b2e1e17ef0661f80f86c4dbd1`,
with identical before/after hashes. It takes **777.55 s**, peaking at
**136,380,416 bytes resident**. A subsequent
[source preimage proof](../../experiments/streets/python_v1/dev40_native_check_fix.json)
shows that only `tests/check_python_street_scene.py` changed; all production
code and numerical tests remain byte-identical. The correction recaptures the
remote object's identity before each incremental Update, after a separate
forced full rebuild intentionally replaces all outputs. The original
`native_42` run generated all 46 fixtures but failed this verification assertion;
it is not a passing combined-scene receipt. The corrected
[focused native check](../../experiments/streets/python_v1/terminal_update_native_01/terminal_update_checks.json)
passes local/full equality, aligned-branch switching, restoration, retained
remote identity across local Updates, and cold-cache output identity.
The [complete native queue](../../experiments/streets/python_v1/dev40_checks.json)
then passes on the current source: **46 fixtures**, **295,216 quads**, all **260**
old mesh parts exact, fresh reopening and independent slopes/interfaces.
The [independent cap audit](../../experiments/streets/python_v1/terminal_cap_native_01/terminal_native_checks.json)
confirms a manifold quad disk, 32 internal whole body edges, zero-height outer
curbs, **0.190 mm** sampled curb error and the retained 15 mm edge. See the
[current combined scene](../../experiments/streets/python_v1/native_43/python_streets_checkpoint.blend)
and [fresh reopening receipt](../../experiments/streets/python_v1/reopen_44/reopen_checks.json).
See [the cap construction and formulas](STREETS_PYTHON_GEOMETRY.md#nonparallel-merged-free-terminal-cuts).

### Dev41 qualified small-edit checkpoint

The completed dev41 core-qualification Python-source fingerprint is
`f101d652a40f5ac08bfc3646ef7a1ac8b0d75fe9273c8bc85e0ac4d0de64b14f`.
The [eight focused tests](../../experiments/streets/python_v1/terminal_unions_unit_03/unit_checks.json)
pass with identical before/after fingerprints, taking 85.73 s with
84,606,976 bytes peak resident memory. They independently check bounded boundary
error, preserved outer corners, metric grading and allocation limits, along
with the existing whole-profile, quad-quality, slope, restoration and cache
tests. Small and intermediate edits are exercised in both interpolation modes.

The [focused Blender check](../../experiments/streets/python_v1/terminal_update_native_02/terminal_update_checks.json)
also passes both modes with float32 source coordinates: local/full equality,
retained unrelated objects, aligned-branch switching, exact restoration and
cold-cache identity. It takes 105.93 s with 256,442,368 bytes peak resident
memory. These are qualification runs, not controlled performance comparisons.

The [complete numerical regression](../../experiments/streets/python_v1/suite_46/suite_checks.json)
now passes **808 tests in 85 modules**, with zero errors/failures, identical
before/after fingerprints and 140,607,488 bytes peak resident memory. The
source-pinned `qualify_dev41.py` queue is now running the combined scene,
saved-buffer comparison, fresh reopening, independent surface/interface audits
and exact upgrades of the three real-network scenes. The first seven jobs now
pass, including the
[46-fixture dev41 scene](../../experiments/streets/python_v1/native_44/python_streets_checkpoint.blend),
[saved-buffer audit](../../experiments/streets/python_v1/gallery_quality_44/import_quality_checks.json),
[fresh reopening](../../experiments/streets/python_v1/reopen_45/reopen_checks.json),
[physical slopes](../../experiments/streets/python_v1/slopes_gallery_dev41_01/slope_checks.json),
[interfaces](../../experiments/streets/python_v1/gallery_interfaces_07/interface_checks.json)
and [terminal cap](../../experiments/streets/python_v1/terminal_cap_native_02/terminal_native_checks.json).
All 263 dev40 mesh parts are exactly preserved in position and topology. Native
generation/Update takes 74.90 s with 377,692,160 bytes peak resident memory.
The three real-scene upgrades also pass, completing the
[ten-job queue](../../experiments/streets/python_v1/dev41_checks.json) on the
unchanged fingerprint. Each upgrade preserves all authored source data,
published positions/connectivity and semantic output keys:

| Current-version scene | Quads | Upgrade worker peak resident |
| --- | ---: | ---: |
| [Berlin](../../experiments/streets/python_v1/berlin_upgrade_dev41_01/python_street_import_upgraded.blend) | 398,972 | 603,336,704 bytes |
| [Canal Street](../../experiments/streets/python_v1/canal_upgrade_dev41_01/python_street_import_upgraded.blend) | 21,818 | 279,855,104 bytes |
| [Nashville](../../experiments/streets/python_v1/nashville_upgrade_dev41_01/python_street_import_upgraded.blend) | 263,213 | 455,913,472 bytes |

Fresh current-version seam edits, incremental/full equality and restoration
also pass in the [three-job edit queue](../../experiments/streets/python_v1/dev41_real_edits.json).
Each tested edit moves both copies of the seam vertex and changes both adjacent
Section widths. The other prepared component stays unchanged. These receipts
are distinct from the prior dev39 runs below.

| Dev41 fresh-edit check | Cold reconstruction | Local Update | Forced full | Restore | Retained outputs | Peak resident |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| [Berlin / Connection 50](../../experiments/streets/python_v1/berlin_dev41_reopen_01/import_reopen_checks.json) | 39.24 s | 1.50 s | 67.36 s | 1.73 s | 319/327 | 643,780,608 bytes |
| [Canal Street / Connection 12](../../experiments/streets/python_v1/canal_dev41_reopen_01/import_reopen_checks.json) | 11.16 s | 11.68 s | 12.20 s | 10.96 s | 14/20 | 277,164,032 bytes |
| [Nashville / Connection 3](../../experiments/streets/python_v1/nashville_dev41_reopen_01/import_reopen_checks.json) | 12.36 s | 0.89 s | 16.17 s | 2.57 s | 173/178 | 479,510,528 bytes |

All workers keep the same 1 GiB process cap and system reserve. These are
qualification measurements under current system load, not a controlled speedup
comparison with dev39. The quickstart points to the audited dev41 gallery.

### Additional construction coverage after the frozen regression

The [higher-valence check](../../experiments/streets/python_v1/higher_valence_04/valence_checks.json)
passes six and seven arms in both Exact and Smooth modes, with unequal widths
and crown rises. The six-arm junction uses 2,016 quads; the seven-arm junction
uses 1,646 / 1,620 respectively. Complete mouth edges, a manifold quad disk,
the highest-crown ceiling and the existing ratio-8 limit pass. Dense analytical
curve-to-published-edge checks measure at most **4.920 mm**, within the 5 mm
contract. Maximum curb-height roundoff is $3.85\times10^{-11}$ m; it is recorded
against the existing $10^{-8}$ m numerical coordinate tolerance, not hidden by
a visual tolerance. These initial numerical checks preceded the maintained
tests and combined native fixtures qualified in the current checkpoint above.

The [two-merged-end check](../../experiments/streets/python_v1/two_merged_ends_02/two_end_checks.json)
passes both modes with straight internal dimensional seams: 4,214 paving plus
188 grass quads, two nonparallel merged end caps, whole body joins and one
closed separator island. The paving correctly has annular Euler characteristic
zero. The [initial probe](../../experiments/streets/python_v1/two_merged_ends_01/two_end_checks.json)
also passes Smooth mode when the width-changing seams bend, but Exact mode
stops at `Exact bent Section seams need an explicit common miter`. That valid
construction remains unfinished. Keeping the dimensional seams straight in the
second probe isolates the two-end cap coverage; it does not resolve or discard
the original Exact-mode failure. Runtime source remains unchanged in both checks.

Start with the [inspection and editing guide](STREETS_PYTHON_QUICKSTART.md).
The [geometry document](STREETS_PYTHON_GEOMETRY.md) explains the construction,
formulas and drawings. The [development record](STREETS_PYTHON.md) retains
earlier experiments, failed checks and unresolved construction branches.

## What is implemented

Blosm's normal import can select **Python (experimental)** or **Geometry Nodes**.
The Python backend consumes the shared prepared Street, Section, intersection,
Connection and Bundle records. It generates ordinary Blender meshes, not a
Python wrapper around a per-fixture GN tree. Source polylines and physical
parameters remain editable. Pressing **Update** regenerates their dependent
pavement, intersections, separators and demonstrated furniture anchors.

This record separates a passing example from general support. Neither a
diagnostic nor a retained previous mesh proves that the rejected configuration
was generated. The implementation is experimental; the limitations below remain
visible when assessing the initial deliverable and the continuing project.

## Requirements and inspected evidence

| Requirement | Implementation and evidence | Verification scope |
| --- | --- | --- |
| Preserve prepared identities | [adapter](../way/street_network/python_backend/adapter.py), [source records](../way/street_network/python_backend/model.py); `test_python_street_adapter.py`, `test_python_street_core.py` | Directed endpoint incidence, Street membership, ordered multi-Section Bundle traversal, serialization and duplicate/missing identity diagnostics. Connectivity is not inferred from proximity. |
| Source grouping independent of semantics | [native grouping check](../../experiments/streets/python_v1/grouping_dev39_02/grouping_checks.json) | Grouped, per-Street and mixed input objects generate identical results. Joining source objects does not itself change semantic output. |
| Reusable mesh parts | [PatternLibrary](../way/street_network/python_backend/patterns.py); `test_python_street_core.py`, `test_python_street_bundles.py`, `test_python_street_work_reuse.py` | Immutable direct-quad connectivity and unit parameter coordinates, cache hits/eviction, bounded allocations, exact comparison with uncached construction. Dimensional/world-space meshes are not substituted for reusable patterns. |
| Six intersection families | [49-fixture scene](../../experiments/streets/python_v1/native_45/python_streets_checkpoint.blend), [generation and edit records](../../experiments/streets/python_v1/native_45/scene_checks.json), [fresh reopening](../../experiments/streets/python_v1/reopen_46/reopen_checks.json) | Named base fixtures `T`, `three_way`, `oblique`, `five_way`, `unequal_widths`, `divided_T`, plus `six_way`, `seven_way`, bent, seamed, overlapping and dimensional regressions. All 49 controllers reopen and update. |
| Rounded zero curbs, complete profiles and direct quads | [saved-buffer audit](../../experiments/streets/python_v1/gallery_quality_45/import_quality_checks.json), [whole-edge interfaces](../../experiments/streets/python_v1/gallery_interfaces_08/interface_checks.json); `test_python_street_junctions.py`, `test_python_street_conditions.py`, `test_python_street_sampling.py` | All 313,576 gallery quads pass convexity, overlap, crown-envelope and edge-ratio checks. Interface checks compare complete float32 mouth edges, not just nearby points. Requested analytic curves remain the sampling reference. |
| Controlled surface and entrance slopes | [independent gallery audit](../../experiments/streets/python_v1/slopes_gallery_extended_01/slope_checks.json) and real-network audits below; `test_python_street_entrance_slopes.py`, `test_python_street_height_projection.py`, `test_python_street_drainage.py` | Published bilinear-face gradients checked independently, including explicit float32 uncertainty. The 3% guard is an **entrance longitudinal** bound on the finite-strip branch; it is not a global crossfall limit. Incoming narrow, high-crown profiles can themselves have much larger transverse slopes. |
| Dynamic separators and merged Bundle paving | [Bundle construction](../way/street_network/python_backend/bundles.py), [union tests](../tests/test_python_street_union.py), [attached union tests](../tests/test_python_street_attached_unions.py), [dimensional Bundle tests](../tests/test_python_street_dimensional_bundles.py) | Curved varying gaps, appearing/disappearing islands, overlapping pavement, reversed members, compatible Section chains, finite skew mouths and the documented free-end branches. Remaining configurations are listed below. |
| Explicit, dependency-driven Update | [coordinator](../way/street_network/python_backend/updates.py), [scene publication](../way/street_network/python_backend/scene.py), [Update panel](../way/street_network/python_backend/ui.py); real edit receipts below | Authored fingerprints, dependency closure, task signatures, introduced spatial collisions, incremental/full equivalence, unaffected-object identity and exact restoration. No continuous geometry observer. |
| Failure safety, persistence and Undo | [publication rollback](../../experiments/streets/python_v1/publication_dev39_02/publication_checks.json), [Undo/Redo](../../experiments/streets/python_v1/undo_dev39_02/undo_checks.json), [fresh reopening](../../experiments/streets/python_v1/reopen_42/reopen_checks.json) | Includes a late metadata-publication failure after allocation: old output and manifest survive, temporary objects/meshes are removed, corrected retry succeeds. Actual Blender history and fresh-process reconstruction are exercised. |
| Source-relative furniture dependencies | [Placement records](../way/street_network/python_backend/model.py), [FurnitureProvider](../way/street_network/python_backend/providers.py), fresh-reopen receipt above | Arrow empties demonstrate station/side/offset anchoring, position and width changes, surface-height queries and persistence. They do not represent a completed lamp/tree asset-placement system. |
| Extensible classes and provider lifecycle | [backend registry](../way/street_network/python_backend/backend.py), [native extension check](../../experiments/streets/python_v1/extension_dev39_01/extension_checks.json) | A registered custom mesh provider receives dependencies and a custom material, updates, saves/reopens and matches a full build. Missing extension code produces a diagnostic and preserves previous output; scene data never supplies executable import paths. |
| Universal asset loading and backend choice | [asset service](ASSET_LIBRARIES.md), [31 native library tests](../../experiments/streets/python_v1/assets_dev39_01/asset_checks.json), [actual UI/renderer/backend check](../../experiments/streets/python_v1/backend_dev39_03/backend_checks.json) | Base/Pro availability, bundled/external identity, datablock types, storage/recovery and materials. Both backends build the same prepared T input into 2,304 quads. Python also imports from a material-only external library with no GN trees. This does not claim complete GN-network qualification. |
| Bounded cache lifetime and modest-memory execution | [native lifetime check](../../experiments/streets/python_v1/lifecycle_dev39_02/lifecycle_checks.json), worker receipts and measurements below | Pattern LRU defaults to 8 MiB per live engine; published buffers are released from the numerical cache. Deleted controllers are pruned on Update; opening a file and unregistering release engines. Live caches remain reusable. |
| Normal real OSM imports | Import, saved-buffer and edit receipts below | Three independent files, both components complete in each. Runtime code does not select an algorithm by a sample filename or fixture label. |

## Exact source versions

The fully qualified **dev39** repository Python-source SHA-256 is
`ff6bd7b2dd09ba6ed99b8e164f2cb36b2f6882353a2e6360c6e2d1102c1cb33a`.
`tests/check_python_street_suite.py` hashes every repository `*.py` filename and
its bytes before and after its run. This fingerprint includes tests and adapters,
not only numerical construction. The [final broad rerun](../../experiments/streets/python_v1/suite_44/suite_checks.json)
passes **800 tests across 84 modules**, with zero errors/failures and equal
before/after fingerprints. Its worker takes **255.35 s**, peaking at
**137,256,960 bytes resident** and **150,421,504 bytes sampled private memory**.
The [ten-job final qualification queue](../../experiments/streets/python_v1/final_dev39_checks.json)
also completes on that unchanged fingerprint. The earlier
[800-test run](../../experiments/streets/python_v1/suite_43/suite_checks.json)
passed on source
`0b27e00c8ad2b9fe1e1ea19b11a8a9fb6a36db1bcdbd50713d55548781a2d981`.

The difference is the runtime cache-lifetime hook in `scene.py`/`ui.py` and its
native test. It does not change geometry, source schemas or the implementation
version, which remains `0.1.0.dev39`. The current-source reopen, history,
publication, grouping, backend, lifetime and extension checks pass. The gallery
and three normal imports were constructed with the earlier fingerprint; this
distinction is retained rather than claiming they were regenerated after an
unchanged-geometry lifecycle fix.

## Real inputs, cost and memory

Each file completes both prepared components. Saved-buffer quality checks pass;
all 328 Berlin, 21 Canal and 185 Nashville mesh parts remain exact from their
respective prior qualified baselines. The edited Canal copy has 21,852 quads and
retains the requested 10 m corner.

| Input | Whole-scene quads | Import | Cold reconstruction | Local Update | Forced full | Restore | Retained output objects |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Berlin | 398,972 | 267.24 s | 149.65 s | 4.08 s | 247.84 s | 5.98 s | 319/327 |
| Canal Street | 21,818 | 40.27 s | 34.89 s | 39.32 s | 40.26 s | 36.42 s | 14/20 |
| Nashville | 263,213 | 73.22 s | 47.20 s | 3.25 s | 58.64 s | 9.03 s | 173/178 |

These are observed single-thread runs, not controlled performance comparisons
with older versions. The local edit changes both copies of an explicit Section
seam and both widths by 5 cm. Berlin uses Connection 50, Canal Connection 12,
Nashville Connection 3. The independent component is unchanged in every check.
The earlier dev37 Berlin receipt tests a broader **combined seam/Bundle edit**;
the latest seam-only timing must not be presented as that same workload.

- **Berlin:** [normal import](../../experiments/streets/python_v1/berlin_dev39_01/import_checks.json), [quality](../../experiments/streets/python_v1/berlin_dev39_quality_01/import_quality_checks.json), [edit/reopen](../../experiments/streets/python_v1/berlin_dev39_reopen_01/import_reopen_checks.json), [published slopes](../../experiments/streets/python_v1/slopes_berlin_dev39_01/slope_checks.json).
- **Canal:** [normal import](../../experiments/streets/python_v1/canal_dev39_01/import_checks.json), [quality](../../experiments/streets/python_v1/canal_dev39_quality_01/import_quality_checks.json), [edit/reopen](../../experiments/streets/python_v1/canal_dev39_reopen_02/import_reopen_checks.json), [published slopes](../../experiments/streets/python_v1/slopes_canal_dev39_01/slope_checks.json), [edited-copy slopes](../../experiments/streets/python_v1/slopes_canal_edited_dev39_01/slope_checks.json).
- **Nashville:** [normal import](../../experiments/streets/python_v1/nashville_dev39_01/import_checks.json), [quality](../../experiments/streets/python_v1/nashville_dev39_quality_01/import_quality_checks.json), [edit/reopen](../../experiments/streets/python_v1/nashville_dev39_reopen_01/import_reopen_checks.json), [published slopes](../../experiments/streets/python_v1/slopes_nashville_dev39_01/slope_checks.json).

Current qualification workers use a **1 GiB resident/private process cap** and
**1 GiB system reserve**, one worker at a time. The latest real-edit peaks are
643,743,744 / 276,520,960 / 480,481,280 bytes resident respectively. The gallery
reopen peaks at 392,318,976 bytes. Blender startup, unrelated scene contents,
external assets and the operating system require their own memory; this is
measured constrained execution, not a universal whole-computer memory guarantee.

Whole-scene counts are not per-Bundle counts. In Berlin, Bundle 11 uses 5,129
paving plus 2,712 grass quads, and Bundle 12 uses 2,231 plus 1,145: **7,841** and
**3,376** respectively. Resolution 1 is the coarsest supported setting; spacing
is `cell_size / resolution`, with local refinement where boundaries demand it.

## Packaging and reproducibility

The feature-independent [asset loader](../app/asset_library.py) remains shared
with GN and other Blender features. `app/asset_store.py` is the separate existing
building-part catalog and is unchanged. The Python renderer needs pavement and
grass materials; it does not need to append GN groups to create its meshes.

| Asset in `blosm/assets` | Current revision / role | SHA-256 |
| --- | --- | --- |
| `prochitecture_streets_v8.blend` | 8.14; default Base library, including the two Python material roles and the retained GN solver | `5b011f46b1ea6111342b2a05d03dbd47e6597815932de8d2750bc519fb7da7d9` |
| `prochitecture_streets_original_spans.blend` | 1.14; optional original-span GN preview | `a11946f361ae33584a8ea0f8dbdf3eccbf02ed451b16f9643e3c151a2ef35f91` |

The default blend must accompany an installable distribution, or **Street
library** must select a compatible external library. A source-only checkout
without either cannot use the default materials. Binaries and generated JSON
are intentionally unstaged under the user's repository policy. This record
does not stage, commit, merge or publish a distribution.

The [current dev43 recovery checkpoint](../../experiments/streets/python_v1/checkpoints/dev43_v1/streets_python_dev43_checkpoint.zip)
contains the complete current addon source, required libraries, documentation,
53-fixture scene, tested addon ZIP and selected qualification receipts. Its
[archive check](../../experiments/streets/python_v1/checkpoints/dev43_v1/checkpoint_checks.json)
verifies every payload hash and reconstructs the exact current Python-source
fingerprint. Extract into a new directory; keep `buildings_gn` separate. The
companion addon ZIP has passed the isolated installation checks above. The
recovery archive is a local handoff artifact, not a public release, and the full
historical evidence tree remains in the original workspace.

The historical [dev42 recovery checkpoint](../../experiments/streets/python_v1/checkpoints/dev42_v2/streets_python_dev42_checkpoint.zip)
contains that version's complete addon source, required libraries, documentation,
51-fixture scene and selected qualification receipts. Its
[archive check](../../experiments/streets/python_v1/checkpoints/dev42_v2/checkpoint_checks.json)
verifies every payload hash and reconstructs the exact current Python-source
fingerprint from the archived bytes. Extract into a new directory; keep
`buildings_gn` separate. This is a local recovery artifact, not an installed or
published release. It does not include the full historical evidence tree.

A historical [dev41 recovery checkpoint](../../experiments/streets/python_v1/checkpoints/dev41_v1/streets_python_dev41_checkpoint.zip)
preserves all 1,171 Python source files, tracked addon files, the required street
libraries, current documentation, the 46-fixture inspection scene and selected
qualification receipts. The [checkpoint check](../../experiments/streets/python_v1/checkpoints/dev41_v1/checkpoint_checks.json)
verifies every archived payload hash and reconstructs the exact `f101d652…`
all-Python fingerprint from the archive. It contains 1,296 archived source/asset/
evidence files, is 19,620,576 bytes, and has SHA-256
`d3dd5a4284966d6da7f7e167b6d534aee822b0a985f336a4e2307297f2d7c037`.
It is a recovery artifact outside Git, not an install-qualified or published
release. The full historical evidence and repository history remain in the
original workspace. Later documentation edits do not alter this frozen archive.

The [49-fixture recovery overlay](../../experiments/streets/python_v1/checkpoints/dev41_extended_v1/streets_python_dev41_fixture_increment.zip)
extends that baseline with the four changed/added test files, current
documentation/drawing, current inspection scene and selected new receipts.
Restore the baseline first, then overlay it. The
[overlay check](../../experiments/streets/python_v1/checkpoints/dev41_extended_v1/checkpoint_checks.json)
verifies all payload hashes and reconstructs the exact `689cd567…` Python
fingerprint from both archives together. The overlay is 3,763,620 bytes with
SHA-256 `85ee4d96d697ec260ec46bf80c25d8c9915d87013772f282276e8f9c217291f9`.
It is also a local recovery artifact, with no install or publication claim.

From the addon root, run the numerical suite with ordinary Python and NumPy:

```powershell
python -B tests/check_python_street_suite.py --output C:/dev/workspace/tmp/my_street_suite
```

For native checks use Blender 5.2 in background mode, one thread, and an output
directory outside the addon. `--python-exit-code 1` makes failed assertions fail
the process:

```powershell
& $blender --background --factory-startup --threads 1 --disable-autoexec --python-exit-code 1 --python tests/check_python_street_scene.py -- --output C:/dev/workspace/tmp/my_street_gallery
& $blender --background --factory-startup --threads 1 --disable-autoexec --python-exit-code 1 --python tests/check_python_street_reopen.py -- --scene C:/dev/workspace/tmp/my_street_gallery/python_streets_checkpoint.blend --output C:/dev/workspace/tmp/my_street_reopen
& $blender --background --factory-startup --threads 1 --disable-autoexec --python-exit-code 1 --python tests/check_python_street_lifecycle.py -- --output C:/dev/workspace/tmp/my_street_lifecycle
```

Set `$blender` to the installed executable first. Additional maintained native
scripts cover normal import, import quality, import reopening, source grouping,
Undo/Redo and publication. Their arguments are documented in their parsers and
the [development record](STREETS_PYTHON.md#inspect-and-reproduce). The workspace's
`prochitecture_gn/experiments/streets/blosm_integration/bounded_blender.py` wrapper
adds the measured cap/reserve checks. Set numerical thread variables to 1 and
run heavy checks serially. Qualification receipts and actual-mesh illustrations
live under `../experiments/streets/python_v1` relative to the addon folder;
they are development evidence, not addon runtime dependencies.

## Remaining construction and release work

The following are outstanding, not supported inputs disguised as errors:

- Nonmonotone overlapping Bundle charts and broader merged-terminal
  configurations. The qualified straight-seam two-end numerical case is not
  a proof for every bent dimensional-seam or terminal arrangement.
- Profile-changing or multiply claimed dimensional seams consumed by a junction,
  broader internal physical junction arrangements inside overlapping members,
  and general interior Street islands. Dev65 qualifies the paired curved-run
  crossing above; it does not establish every internal region graph or an
  overlapping incoming-mouth profile. The passing same-profile junction-owned
  width-seam example does not establish these broader branches either.
- Dev42 integrates Exact-mode bent width-changing seams into `SectionTransition`
  and Bundle curb handling. The captured failure and an ordinary seam now pass
  focused numerical/native checks, all 818 numerical regressions and the complete
  gallery/reopen audits. Dev43 resolves the four recorded 150-degree Exact
  failures and passes all 827 regressions plus 53-fixture gallery/reopen checks.
  Dev50/dev51 subsequently qualify the 60 m, 150-degree Smooth case, including
  native editing and fresh reopening. Other sharp Smooth configurations and
  route-spanning mixed Exact/Smooth seams still need construction work. Local
  mixed joins now have the dev69/dev71 construction and dev72 matching-profile
  short/sharp height fallback above; this does not establish every profile or
  arbitrary-angle configuration. The 30 m Smooth
  example fails its source offset-curvature proof; it must not be accepted by
  loosening that proof.
- General near-contact/tip grading and broader arbitrary-edit coverage. Input
  centerlines crossing or touching remain the explicitly agreed invalid case;
  ordinary pavement overlap is a construction problem to solve.
- Further reductions in complex Canal Update cost and wider independent real
  network coverage. The current local distribution and asset package are
  qualified above; measured Update timings do not establish universal
  responsiveness.
- Terrain/graded centerlines and isolated closed Section loops require separate
  providers. The initial provider explicitly reports nonzero-Z controls or a
  closed Section rather than silently generating an unrelated flat path.

The historical
[1 mm control-move failure](../../experiments/streets/python_v1/terminal_small_edit_01/small_edit_checks.json)
is resolved within the qualified single-terminal family by dev41's bounded
front sampling, local grading and checked interior initialization. The
[32-move development sweep](../../experiments/streets/python_v1/terminal_move_matrix_03/move_checks.json),
eight permanent tests, both-mode focused Blender checks and complete numerical,
combined-scene and reopening checks now pass. That closes the captured editing
gap; it does not establish all other Bundle configurations listed above.

### Historical construction gap and dev40 increment: merged free terminal directions

The [source-pinned probe](../../experiments/streets/python_v1/merged_terminal_edit_01/merged_terminal_edit.json)
starts with `one_sided_free_cap` and changes Section 2's penultimate control
from `(30, 4, 0)` to `(30, 4.05, 0)`. Widths, source IDs, endpoint positions and
profile rises do not change. Continuous centerline-separation certification
passes, and each member path passes its regular offset certification. The upper
terminal tangent becomes approximately `(0.9999875, -0.00499994)` while the lower
one remains `(1, 0)`.

Both Exact and Smooth rejected the resulting separator boundary under dev39.
The aligned merged-cap branch is no longer eligible; the separated fallback
cannot construct the pavement union and its disappearing separator. Previous
successful mesh buffers and manifest survive exactly, which is correct failure
handling but not a solution to this valid edit. The probe worker completes at
49,741,824 bytes peak resident, under the same 1 GiB process cap/reserve.

The required construction must retain the actual terminal cross-section segments
and lateral curbs, form their exterior union, and join that cap to the shared
Bundle body with direct quads. Rotating the source ends into a common plane,
silently trimming them to a common rectangle, lowering profile rises or simply
returning two overlapping road meshes would change the required geometry.
The footprint/height/quad checks and local/full equality must qualify this
branch before the captured edit can be called supported. Dev40 now has a
production fallback and passing focused numerical tests for that construction;
its native and full-suite qualification is recorded separately when complete.

The [bounded footprint prototype](../../experiments/streets/python_v1/merged_terminal_boundary_02/terminal_boundary.json)
now forms that exterior union for the captured Smooth input with the existing
`RegionUnion`, preserving each retained curb/terminal segment's source and
fractions. A computational cut at X = 35 bounds only the temporary end patch;
it does not trim the authored member. The union has one positive outer loop,
90.1160 square metres, and retains terminal segments from both Sections. It
contains a **15 mm exposed curb segment** where the nonparallel end cuts meet.
That small physical feature needs deliberate quad grading. The prototype peaks
at 33,574,912 bytes resident. It verifies footprint preparation only; profiled
surface construction, body seams, quads and Update integration remain undone.

The six-stage initial plan explicitly postpones full terrain alignment, later
transport features and finished furniture assets. They remain in the broader
street project. Both backend choices are preserved; this Python deliverable does
not turn the incomplete GN-network milestone into a finished one.

Source geometry editing currently preserves the number of controls and entity
identities. Keep duplicate incident endpoint controls coincident when moving a
junction or Section seam. Arbitrary topology changes need managed operations;
they are not safely inferred from vertex indices. The provider uses flat ground
and a common datum, retains whole entrance profiles and checks the highest
incoming crown. See the geometry document for the exact profile and slope rules.
