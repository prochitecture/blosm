# Shared Blender asset loading

`app.asset_library` is independent of streets, buildings and nature. A feature
declares the libraries and datablocks it needs; the loader resolves availability,
loads them and returns the actual Blender IDs. The feature decides how to use or
instance those IDs.

```python
from .app.asset_library import AssetCatalog, AssetRef, LibrarySpec, BlenderAssetLoader

catalog = AssetCatalog([
    LibrarySpec('streets.native.v8', 'prochitecture_streets_v8.blend', revision='8.10'),
    LibrarySpec('trees.pro', 'trees_pro.blend', edition='PRO', revision='1'),
    LibrarySpec('paving.custom', '/my/library/paving.blend', source='external', revision='3'),
])
loader = BlenderAssetLoader.for_app(app, catalog)
material = loader.load(AssetRef('paving.custom', 'materials', 'Paving / Limestone'))
tree = loader.load(AssetRef('streets.native.v8', 'node_groups',
                            'Prochitecture Streets / Editable network v8'))
```

Bundled paths are relative to `blosm/assets`. They may not escape that directory.
`BASE` means available to both editions; `PRO` requires `app.isPro`. Product
packaging must include only the libraries distributed for that edition. A missing
library or unavailable Pro library produces a specific error, without silently
choosing an unrelated fallback. Asset availability does not change feature/UI
license rules.

External libraries use absolute paths or Blender `//` paths relative to the saved
scene. External source does not imply a different datablock API. A feature may
replace its catalog entry with a compatible external library under the same key.

Use Blender collection names such as `objects`, `collections`, `node_groups`,
`materials`, `meshes`, `images`, `worlds` or other kinds exposed by
`bpy.data.libraries.load`. A shader-node group is a `node_groups` asset; a material
with its shader is a `materials` asset. The loader does not link objects or
collections into the scene. Link/instance or copy them explicitly when appropriate.

`load_many()` resolves all libraries and checks all requested names before loading
missing IDs. It batches requests per library. Missing-name preflight is atomic;
an unexpected Blender/file I/O failure during the subsequent load is reported and
does not promise transactional rollback of Blender's internal dependency loading.

Direct construction defaults to appending. `for_app()` uses the import setting:
Automatic links certified, fingerprint-pinned revisions and embeds other
libraries; Embed in scene and Retain library revisions are explicit alternatives.
Per-call `link=True` or `link=False` overrides remain available.
Appended IDs carry a source/revision/kind/name tag, so a
later loader or reopened scene can reuse the right ID even if Blender renamed it
because of a collision. An unrelated local datablock with the same name is never
reused. Linked assets use `load(..., link=True)` and retain Blender library semantics;
their identity is the library path and original name. Publish linked revisions to
different versioned paths, or reload the library explicitly through Blender. The
service does not silently reload linked data that other scene objects are using.

## Retaining linked revisions within a modest memory budget

Linking a shared GN library can avoid the serialization peak of saving thousands
of appended nodes. The retained-library store is feature-independent. The
current public street catalogs remain uncertified, so Automatic still embeds
those libraries; the private candidate has not been promoted.

```python
catalog = AssetCatalog([
    LibrarySpec('street.programs', 'published_street_programs.blend',
                revision=published_revision, sha256=published_sha256,
                self_contained=True),
])
loader = BlenderAssetLoader.for_app(app, catalog, linked_root=retained_library_directory)
program = loader.load(AssetRef('street.programs', 'node_groups', program_name), link=True)
```

The publisher supplies the revision and SHA-256; `self_contained=True` must be
backed by an audit of that exact file. The development command
`audit_library_dependencies.py` opens the whole library in a fresh Blender worker
and rejects unpacked file paths or nested linked libraries. Copying a library
with unresolved relative textures or other file dependencies is not supported
by the single-file store. Such assets can continue using the ordinary loader at
their original location; a future package format must retain their dependencies.

The caller chooses a durable, writable store outside the replaceable addon
installation, or a project-local `libraries` directory. An empty app setting
selects `datafiles/blosm/libraries` under `bpy.utils.resource_path('USER')`.
A `//` directory requires a saved scene; other relative directories are rejected.
The layout is `<root>/<SHA-256>/library.blend`.
Files are copied and checked in 1 MiB chunks, then published atomically without
overwriting another revision. Requests in one batch share the resolution check.
Fingerprint-pinned linking requires this store so a changed source file cannot
silently reuse a stale linked datablock at the same mutable path.

- An addon update installs new bytes under a new fingerprint. Old scenes keep
  their existing library references; loading a new revision does not switch
  existing modifiers or reload their programs.
- An intact retained revision can be used after its original distribution file
  has moved or disappeared. Edition checks still apply when Blosm requests it.
- A missing revision is restored only from bytes matching its fingerprint.
  A damaged retained file produces an explicit error; the service does not
  replace it with the latest version. Restoring missing/corrupted files is an
  explicit maintenance action, not a render-time fallback.
- The store has no automatic eviction: older projects may still need any of
  its revisions. It is a persistent dependency directory, not a temporary cache.
- For a portable project, save relative library paths and copy the scene with
  its `libraries` directory. Moving only the scene does not carry its linked
  dependencies. The Asset libraries panel can restore the exact revision into
  an explicitly selected retained directory and relink the saved reference.

The byte store is implemented in `app.library_store`; the existing
`app.asset_store` building-part catalog is unchanged. Tests cover competing
installers, damaged/missing files, Base/Pro checks, source upgrades, different
revisions with the same datablock names, normal save/reopen, and relocation of a
project with its libraries. This policy addresses library lifetime, not a claim
that all network sizes fit the current process-memory budget.

## Exact-revision recovery

`app.asset_recovery` recognizes saved `<SHA-256>/library.blend` references. The
Blosm sidebar's **Asset libraries** panel offers **Restore exact revision** for a
missing file, and **Verify or restore revision** for a present file. Choose a
backup of the exact revision. By default it restores the saved location; **Use
current retained directory** chooses the configured directory and relinks the
existing Blender library. Save the project to persist a changed link.

Restoration streams through the same atomic verified store. A wrong fingerprint
is rejected before relinking. Existing damaged files are preserved; select a
different retained directory to recover there. Native reload failures are
reported, the old filepath is restored and the verified disk copy remains.
An unmanaged library with no retained digest layout is not guessed or upgraded.
No scene save/reopen or heavy background library check runs from panel drawing.

`tests/test_asset_library.py` includes 31 tests: normal loader/policy checks,
seven recovery boundary checks and an actual missing-library saved-scene
recovery through the registered operator. It also tests wrong-copy rejection,
relocation, a second save/reopen without the distribution copy, generated GN
geometry and preservation of a local editable object property. The native test
uses a disposable factory-startup worker, since it intentionally reopens scenes.

The current street library is generated by
`prochitecture_gn/experiments/streets/blosm_integration/build_assets.py`. Its runtime
code is packaged by `package_runtime.py` into `way/gn`. Installed Blender evaluation
does not read `prochitecture_gn` or any fixture scene. Rebuilding a library requires
a fresh background Blender 5.2 process. Changes to required interfaces must bump
the library revision and its feature-side validation.

Legacy utility loaders retain their previous behavior for existing features.
New/ported features should use this service; they should not introduce separate
street/building/nature copies of the library-loading mechanism.
