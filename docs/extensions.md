# Blender extension development

The `buildings_gn` branch uses Blender's extension format and requires Blender
5.2 or newer, preserving this branch's existing minimum Blender version.
`blender_manifest.toml` is the source of version and package metadata.
Install the release ZIP through **Preferences > Get Extensions > Install from
Disk**, then open the Blosm tab in the 3D View sidebar.

For source development, link the checkout as `blosm` inside a configured local
extension repository, such as `extensions/user_default/blosm`. The manifest and
`__init__.py` must be at the root of that linked directory. After switching from
a legacy branch, restart Blender or refresh the local repository and add-on list,
then enable **Blosm**. A separate link under `scripts/addons` is not needed.

The branch includes PADL, Geometry Nodes building generation, Python streets,
and the bundled ANTLR 4.13.1 runtime. Its internal imports use the extension's
namespace; they do not rely on adding the checkout to `sys.path`.

`tests/check_development_extension.py` verifies repository discovery, isolated
imports, registration, PML compilation, custom setup scripts, data storage,
building/street controls, local 2D and both simple-3D backends, and unregistering.
Pass `--repository` pointing at a test repository containing a `blosm` link to
this checkout and `--artifacts` pointing outside the source tree. Use a separate
factory-startup Blender profile with all `BLENDER_USER_*` directories isolated.
The PADL adapter/parts checks also accept `BLOSM_TEST_PACKAGE` for an installed
extension; the parts check uses `BLOSM_TEST_ARTIFACTS` for its saved scene.

## Build and validate

Release tooling and Base/Pro inventories for the `release` branch live in the sibling
[`blosm_build` repository](../../blosm_build/README.md). This repository's
`blender_manifest.toml` contains runtime metadata only; do not add `[build]` here.
Direct builds with Blender's default exclusions are not release packages.
Those release inventories do not include all `buildings_gn` modules and assets;
do not use them to distribute this development branch without updating and
qualifying its inventories first.

Run these commands from the Blosm repository root with Python 3.11+.
Filesystem paths below are relative to that root; adjust the example Blender
executable location to your setup.

```powershell
$blender = '../tools/blender/blender.exe'

# Official build from committed source, with tests before publication:
python ../blosm_build/build_extension.py --source . --ref release --edition both --test --blender $blender

# Preview uncommitted changes:
python ../blosm_build/build_extension.py --source . --working-tree --edition both --test --blender $blender
```

Output defaults to `../artifacts/blosm-releases`. Each run has its
own logs, input hashes, provenance report and `packages/` directory. With
`--test`, standalone Blosm tests run before candidate packaging and the actual
Base and Pro ZIPs are installed and exercised in isolated Blender profiles.
No release ZIPs are published to `packages/` if any required check fails.
Without `--test`, manifest, ZIP and edition-inventory validation still run.

The install checker accepts `--edition base` or `--edition pro` (default: pro).
It checks preferences, namespaces, offline behavior, local OSM, GLB/B3DM, and
Pro PML when applicable. The broader release tests under `tests/` can additionally
exercise real datasets, UI and external Pro asset packs; see the
[General, Base and Pro test profiles](../tests/README.md).

## Runtime conventions

- Internal and vendored imports remain relative to the installed package.
  Blender's `sys.path` and unrelated global modules are not modified. The
  standalone download subprocess has its own isolated Python environment.
- Preferences use the root `__package__`, including the extension repository.
- Downloads check `bpy.app.online_access`. Revoking access closes an active
  download transport when it is next polled or submitted to. Local input and
  cached imagery remain usable offline.
- With no custom data directory, downloads go to
  `bpy.utils.extension_path_user`, outside the installed extension. User-selected
  asset directories remain external to the extension.
- 3D Tiles rebase private glTF node data before single-precision conversion;
  they do not patch Blender's shared glTF importer classes.
- Custom setup scripts load as `<extension>.setup._user_script`. Use relative
  imports such as `from ..app.blender import app` or
  `from ..manager import BaseManager`; old top-level imports such as
  `from app.blender import app` must be updated.
- Regenerated ANTLR parser files must retain relative imports of the bundled
  `pml.antlr4` runtime. No dependency installation is performed at runtime.

## Repository-wide import audit — 10 November 2026

The vendored graph package still contained absolute `lib.blosm_networkx` imports
after the outer addon modules were migrated. Removing the old sys.path setup
therefore exposed `ModuleNotFoundError: No module named 'lib.blosm_networkx'`.
Do not resolve this by installing NetworkX, restoring sys.path injection or
aliasing Blender's global `lib` module.

The audit covered all **819 runtime Python files**, with separate review of
standalone tests/tools, dynamic imports, resource lookup and generated PML.
Corrections include:

- 658 absolute imports in 222 NetworkX files and five CompGeom files now resolve
  relative to the importing package. Vendored implementation bodies were checked
  against the reference revision and preserved.
- Graph-atlas resource lookup uses `importlib.resources.files(__package__)`, so
  the data file is found inside any installed extension repository namespace.
- Seven street debug imports now reference the local `way.debug` module.
- The generated PML listener consistently imports its sibling parser relatively.
  Preserve this adjustment when regenerating the ANTLR files.
- `lib.pygeos.index_chain` imports Envelope/Quadrant from the existing `shared`
  module, replacing a reference to a nonexistent `constants` module.

The strengthened `tests/test_extension_contract.py` checks nested module names,
namespace packages without __init__.py and literal dynamic package references.
The previous root-package-only scan missed these cases. Absolute standard-library,
Blender and external-library imports remain appropriate. The download worker's
top-level `downloads` imports are intentional: it executes in a separate Python
process and does not add global modules or paths to Blender. Test and CLI entry
points have their own bootstrap code outside addon runtime. The bundled ANTLR
`_pygrun` CLI likewise executes its path setup only as a standalone tool.

One existing optional target remains absent: the legacy standalone
`app.command_line.initOverlay` references `overlay.command_line`. This is a
missing CLI implementation, not an extension namespace issue; this audit does
not add a command-line overlay renderer.

### Regression checks

`tests/test_vendor_imports.py` exercises the graph package under ordinary and
nested `bl_ext.<repository>.blosm` names with a conflicting global `lib` module.
It checks independent graph types across two extension repositories, graph
algorithms, serialization/pickling and the atlas resource without sys.path
changes. It can run with ordinary Python without Blender:

```powershell
python -B tests/test_vendor_imports.py
```

`tests/check_vendor_imports.py` runs those tests and the extension contract in
Blender, followed by the real extension-repository checks. Its native check now
imports WayMap and exercises its graph operations instead of only importing the
street-renderer facade. See the file's `--repository` and `--artifacts` options.
The configured repository must contain a `blosm` link to that checkout. Where
links are unavailable, qualify a copied candidate and run the test script from
inside that exact candidate; verify its source hashes first.

The isolated Blender 5.2 run passed **8 host contracts and 14 native checks**,
including discovery/registration, conflicting module names, PML compilation,
WayMap operations, atlas data, local 2D imports, GN/Python simple-3D imports in
both grouping modes, disable/re-enable and namespace cleanup. Source hashes,
logs and the audit report are under
`artifacts/streets/extension_relative_imports/`. No files in the reserved main
checkout or default Blender installation were changed during qualification.

This is technical extension compatibility. Submission to extensions.blender.org
is a separate publishing process with its own content and review requirements.

References: [extension manifest and packaging](https://docs.blender.org/manual/en/latest/advanced/extensions/getting_started.html),
[add-on migration](https://docs.blender.org/manual/en/latest/advanced/extensions/addons.html),
[Blender add-on guidelines](https://developer.blender.org/docs/handbook/addons/guidelines/).
