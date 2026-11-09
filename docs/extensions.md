# Blender extension release

The release branch uses Blender's extension format and requires Blender 4.5 or
newer. `blender_manifest.toml` is the source of version and package metadata.
Install the release ZIP through **Preferences > Get Extensions > Install from
Disk**, then open the Blosm tab in the 3D View sidebar.

## Build and validate

Build from Git-visible source to avoid including ignored local experiments.
The manifest excludes credentials, development scripts, tests, bytecode and
development assets; only the three shipped `.blend` files are included.
Keep output and validation profiles outside the repository, under
`C:/dev/workspace/artifacts/blosm-extension-compatibility` in this workspace.

```powershell
python tools/build_extension.py --blender 'C:/Program Files/Blender Foundation/Blender 4.5/blender.exe' --output-dir C:/dev/workspace/artifacts/blosm-extension-compatibility/package
blender --command extension validate
blender --command extension validate C:/dev/workspace/artifacts/blosm-extension-compatibility/package/blosm-2.7.28.zip
```

Use an isolated profile for the installation test; it installs and enables the
actual ZIP under `bl_ext.blosm_test.blosm`. This leaves normal Blender preferences
and installed add-ons untouched.

```powershell
$env:BLENDER_USER_RESOURCES = 'C:/dev/workspace/artifacts/blosm-extension-compatibility/test-profile'
$env:PYTHONDONTWRITEBYTECODE = '1'
blender --background --factory-startup --python-exit-code 1 --python tests/check_extension.py -- --package C:/dev/workspace/artifacts/blosm-extension-compatibility/package/blosm-2.7.28.zip --artifacts C:/dev/workspace/artifacts/blosm-extension-compatibility/install-test
python -B -m unittest discover -s tests -p 'test_*.py'
```

Download transport tests use localhost fixture servers and need loopback socket
access. The installation test covers registration, preferences, namespace
collisions, PML, local storage, offline guards, and GLB/B3DM geometry. The direct
glTF and progressive integration tests under `tests/` also verify placement,
attribution, cancellation and failure cleanup.

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
