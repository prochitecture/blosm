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

This is technical extension compatibility. Submission to extensions.blender.org
is a separate publishing process with its own content and review requirements.

References: [extension manifest and packaging](https://docs.blender.org/manual/en/latest/advanced/extensions/getting_started.html),
[add-on migration](https://docs.blender.org/manual/en/latest/advanced/extensions/addons.html),
[Blender add-on guidelines](https://developer.blender.org/docs/handbook/addons/guidelines/).
