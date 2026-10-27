# Python streets: inspection and editing

The latest source is now on `buildings_gn` at `C:/dev/workspace/blosm`.
Runtime dev75 adds [shared rounded Section corners](STREETS_PYTHON_LOCAL_SMOOTHING.md#shared-corner-across-section-seams-dev75)
to the bounded local smoothing construction.
Open the [corrected Karl-Marx-Allee inspection scene](../../streets_dev/experiments/streets/python_shared_rounds/shared_rounds_inspection.blend).
It selects `Street 595` and opens directly at the repaired corner. The
[current 112-fixture gallery](../../streets_dev/experiments/streets/python_shared_rounds/gallery/python_streets_gallery.blend)
reproduces the previous fixture meshes and furniture exactly. Reload the current
addon before editing or pressing Update. See the local smoothing notes for the
focused regression and native verification scope.

## Previous qualified package and examples

The dev72 scenes and ZIP below remain the preceding qualified milestone;
reload the current addon and re-import or Update to use the revised curves.

Use **Blender 5.2** and the dev72 development checkout at
`C:/dev/workspace/streets_dev/blosm`, branch `streets_python`. The separate
checkout at `C:/dev/workspace/blosm` remains on `buildings_gn`. NumPy ships with
Blender; SciPy is used only by optional independent experiments.

Restart Blender or reload the addon after changing development code. The
[dev72 ZIP](../../streets_dev/experiments/streets/python_v1/packages/dev72_01/blosm_streets_python_dev72.zip)
has passed [isolated installation qualification](../../streets_dev/experiments/streets/python_v1/package_dev72_qualification_01.json).
It uses the addon ID `blosm`; use a separate Blender profile when comparing it
with a buildings-development installation. Enable the addon or run the scene's `OPTIONAL - Python street
editing.py` once to register the Streets panel.

<a id="current-development-gallery"></a>

## dev72 inspection scene

Open the [112-fixture dev72 gallery](../../streets_dev/experiments/streets/python_v1/gallery_generated_dev72_01/python_streets_gallery.blend).
It contains all earlier fixtures, arranged without overlaps. Resolution 1 and
exaggeration 10 are the saved inspection defaults. Search the Outliner for:

- `sharp_uniform_mixed_0_3` and `sharp_uniform_mixed_3_0`: short 110-degree
  Exact/Smooth joins with bounded interior height reconciliation. Each join has
  768 direct quads and preserves its entrance profiles, zero curbs and crest.
- `uniform_mixed_bundle` and `uniform_mixed_bundle_overlap`: separated
  carriageways with grass, and one pavement where the carriageways overlap.
- `attached_angular_multiple_3_internal`: curved Bundle members, four automatic
  grass islands and a transverse Street retaining its original Section IDs.
- `translated_exact_corner_seed`: the rounded-corner grading regression.

Suffix 0 denotes Exact polyline, 3 Smooth interpolation. The menu presents these
as two named options. Fixtures deliberately exercise both modes.

The [gallery audit](../../streets_dev/experiments/streets/python_v1/gallery_dev72_audit_01.json)
proves all 110 predecessor source/mesh/furniture/layout sets unchanged, actual
saved mesh quality, all 112 cold Updates and native edits of the two newest
examples. The [complete regression](../../streets_dev/experiments/streets/python_v1/suite_dev72_audit_01.json)
passes all 982 tests. These checks have distinct scopes; the full evidence and
remaining construction limits are in the [verification record](STREETS_PYTHON_ACCEPTANCE.md).

## Everyday workflow

For a streets-only first import:

1. In Blosm preferences, enable **experimental features** and set a writable
   **directory for data**.
2. Choose **OpenStreetMap**, import from **File**, and select **3D realistic**.
3. Enable roads; disable buildings, forests and the other feature layers for
   this demonstration. The bundled street library needs no external building
   or vegetation asset package.
4. Choose **Street backend → Python (experimental)** and press **Import**.

The [actual registered Import-operator check](../../streets_dev/experiments/streets/python_v1/installed_ui_import_dev72_01/operator_import_checks.json)
passes these settings from an isolated installation without replacing setup or
renderer functions. **Geometry Nodes** remains a
separate backend choice over the shared prepared Street/Section/junction/Bundle
data. No sample filename chooses a construction algorithm.

Open **Blosm Python street import report** in Blender's Text Editor. **COMPLETE**
means every prepared component generated. **PARTIAL** identifies unsupported
preparation or construction; inspect the listed reasons before treating the
network as complete.

Select a source centerline, move its vertices or use **Section properties** to
change width, offset, crown rise or interpolation. Then press **Streets →
Update**. Meshes remain unchanged until that action. Related carriageways,
intersections, separators and demonstration furniture arrows regenerate through
their recorded dependencies. Unrelated outputs retain their objects.

Keep incident endpoint copies coincident when moving a junction or Section
seam. Ordinary editing preserves control counts and persistent IDs. Split,
join, delete and reconnect operations need managed topology commands and are
outside this initial editing contract. Edit positions and parameters instead
of scaling source objects. A failed Update retains the last successful output;
correct or undo the source edit, then press Update again.

Use **Exaggeration = 1** for physical heights. **Resolution = 1** is the coarsest
setting; nominal spacing is **Cell size / Resolution**. Higher settings add
work and may reach the documented construction/allocation limits. Complex
affected regions can take seconds to regenerate. This is a manual Update
workflow; no continuously running geometry observer is installed.

## dev72 real-network examples

The [three-input audit](../../streets_dev/experiments/streets/python_v1/dev72_real_receipt_audit_01.json)
passes normal import, saved mesh quality, reopening, coordinated seam/width
editing, local/full equivalence and exact restoration. Open:

| Example | Editable Blender scene | Measured local Update |
| --- | --- | ---: |
| Karl-Marx-Allee | [Berlin](../../streets_dev/experiments/streets/python_v1/berlin_fresh_dev72_01/python_street_import.blend) | 1.66 s |
| Canal Street | [Canal](../../streets_dev/experiments/streets/python_v1/canal_fresh_dev72_01/python_street_import.blend) | 7.63 s |
| Nashville | [Nashville](../../streets_dev/experiments/streets/python_v1/nashville_fresh_dev72_01/python_street_import.blend) | 0.94 s |

These are individual qualification measurements with one numerical thread,
not a universal timing guarantee or controlled comparison. All nine workers
stay below 0.75 GiB resident/private usage with a 1 GiB system-memory reserve.
The files are examples, not a required input set or a restriction on imports.

## Materials and external assets

The default material library is `assets/prochitecture_streets_v8.blend`.
**Street library** can select a compatible external Blender library. Base/Pro
availability and material identity use the shared [asset service](ASSET_LIBRARIES.md).
An external Python material library can supply pavement and grass roles without
containing a GN solver. Demonstration furniture arrows verify dependency and
placement behavior; finished furniture assets are later providers.

## Reproduce the complete inspection gallery

The [shared catalogue](../tests/_python_street_catalogue.py) maintains authored
inputs. The [generator](../tests/generate_python_street_gallery.py) creates all
112 fixtures through the normal Update operator, without importing saved output
meshes. Use a fresh output directory. In PowerShell:

```powershell
$blender = 'C:/Program Files/Blender Foundation/Blender 5.2/blender.exe'
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
& $blender --background --factory-startup --threads 1 --disable-autoexec --python-exit-code 1 --python 'C:/dev/workspace/streets_dev/blosm/tests/generate_python_street_gallery.py' -- --output 'C:/dev/workspace/tmp/my_python_street_gallery'
```

Open `python_streets_gallery.blend` in the new output directory. For a smaller
inspection, repeat `--case`, for example `--case T --case sharp_uniform_mixed_0_3`.
The optional `--reference-catalogue` verifies development captures; ordinary
generation needs no earlier scene or captured record. Runtime caches are
released between fixtures, and layout uses actual fixture bounds.

See the [geometry and formulas](STREETS_PYTHON_GEOMETRY.md),
[development plan/history](STREETS_PYTHON.md) and
[qualified scope and limits](STREETS_PYTHON_ACCEPTANCE.md). Full terrain handling,
later transport features, production furniture and broader unsupported geometry
remain separate work. The historical records preserve earlier failed attempts;
the current scene and checks above are the inspection starting point.
