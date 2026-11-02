# Street Geometry Nodes backend: paused and removed

The experimental street GN backend was removed from `buildings_gn` on
2 November 2026 so current development can concentrate on Python streets.
The Legacy renderer was subsequently removed too. **Python is the only installed
street backend**; the [backend registry and lifecycle](STREET_BACKENDS.md) retain
a straightforward extension point for future GN and other implementations.
Building Geometry Nodes functionality remains available.

## What changed

- Removed the packaged street GN implementation under `way/gn`, its runtime
  registration/update handlers, native providers, source authoring and catalogs.
- Removed the street GN backend menu entry, construction/connected-junction
  preview controls and GN library selector. Python physical controls, its Update
  panel, material selector and the universal asset-storage controls remain.
- Removed the local `prochitecture_streets_v8.blend` and
  `prochitecture_streets_original_spans.blend` runtime libraries. They were
  untracked generated artifacts, so their removal is not a Git binary deletion.
- Moved the two preparation records/adapters still used by Python to
  `way/street_network/contracts/inputs.py` and `blosm_adapter.py`, preserving
  their numerical/topology behavior. Shared physical/source/ownership planners
  and their backend-independent tests remain; they do not construct GN trees.
- Reduced `way/street_network/assets.py` to material declarations. Retained
  `assets/prochitecture_street_materials.blend`, generic asset loading/storage,
  shared street preparation, Python mesh generation and all building GN code.
- Retained Python numeric selection 0. Old GN value 1 and Legacy value 2 are
  retired and migrate to Python. New backends must use new persistent numbers.
  Old `streetConstruction` / `streetAssetFile` data is no longer exposed by RNA.

Both retired street implementations remain in Git history. Shared building GN
functionality and generic geometry-node utilities were preserved.
Saved GN street objects are not converted into Python output by loading a file;
use the preserved addon revision for further editing with the old GN tooling.

## Recovery does not require a checkpoint folder or old Blender libraries

Start with the [verified GN restart and regeneration guide](../../prochitecture_gn/docs/streets/blosm_integration/GN_RESTART.md)
in the sibling source repository. It contains architecture/formula/diagram
links, exact library recipes, loading/qualification commands and known limits.

| Preserved source | Commit |
| --- | --- |
| `prochitecture_gn` generators, materials and recipes | `07720eac551ef55ac87e2b8b5eebe097e8e91f53` |
| Blosm integration before removal | `7c7e5d388fe27f0871d1d151e744d045627ae86b` |
| Git-only recovery tools | PGN `441e1c0793cdb7c108537ea571073365ee8b3cd1` |
| Complete regeneration/verification record | PGN `5881447` |

Git-only regeneration of v8 **8.19**, original spans **1.65** and shared
materials **1.0** passed exact full node/material/Text-catalog comparisons,
append/link contracts, Python material Update and script-disabled reopening
before deletion. The temporary checkpoint folder was then removed.
Keep these repository histories and the specified Blender 5.2 toolchain.
A normal commit removing the backend does not remove the earlier Git sources.

Recover and test the earlier integration in an isolated development checkout,
then reconcile it with the evolved Python contracts and GUI. Do not copy an old
addon tree over current building/Python development. The old GN seam limitation
and incomplete general-network scope remain documented in the restart guide.

## Removal verification

`tests/test_street_backend_registration.py` checks real restricted-context addon
registration, absent street GN files/modules/handlers and controls, preserved
building GN integration, stable Legacy IDs and retired-selection migration.
`tests/test_street_material_library.py` exercises the Python-only material
installation, Update and reopening without either street GN library.
The numerical suite and a fresh real-network import are checked independently;
the final results are also recorded in the Python acceptance notes.

Removal verification completed: 932 host tests, three native registration checks,
four material checks plus reopening, two ground-scope checks, and a fresh Berlin
import with exact prior geometry and successful reopen/local/full Update. See
[the detailed acceptance record](STREETS_PYTHON_ACCEPTANCE.md#python-only-backend-cleanup-validation--2-november-2026).
