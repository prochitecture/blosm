# Street rendering backends

Updated 4 November 2026. **Python is the only installed street backend.**
Both the old Legacy renderer and the experimental GN renderer are removed from
`buildings_gn`. Building GN, shared street preparation, street materials and the
feature-independent asset service remain available.

The [terrain placement stage](STREETS_TERRAIN.md) now offers profile-preserving
Python placement or a small shared direct-projection GN modifier. This modifier
acts on Python-generated meshes; it is independent of street backend selection.

## Stable entry point and common preparation

`way/renderer_streets.py` retains the public `StreetRenderer(app, *, backend=None,
design=None, features=None, asset_loader=None)` entry point. It selects a factory
from `way/street_network/backends.py` and delegates the import lifecycle.
Import setup and network preparation have no Python-versus-GN branches.

All registered street backends receive the same prepared Streets, Sections,
Bundles, original centerline positions, provenance/connectivity audit and
boundary context. Built-in numerical street styles are available independently
of the selected renderer. Normal preparation disables pixel correction; an
explicit custom PML style remains a data-preparation choice rather than a
Legacy-renderer requirement. Street-only import no longer needs Legacy
`style/street/main.pml` or `asset_info/street.json` files.

| Lifecycle method | Responsibility |
| --- | --- |
| `prepare()` | Check dependencies, load required assets and initialize import resources |
| `render(manager, data)` | Read the shared prepared network, generate output and return the backend's context |
| `finalize()` | Finish publication of successful output |
| `cleanup()` | Release temporary backend resources |

Implementations inherit `StreetRenderBackend`. The factory receives `app` and
keyword arguments `design`, `features`, `asset_loader`. It owns backend-specific
assets, geometry, editing tools and caches; it must return a new renderer for
this import, not a cached mutable engine shared between imports. Existing
import transactions and per-backend rollback/publication policies remain the
implementation's responsibility. The facade forwards return values/errors and
the existing caller controls the lifecycle; it does not change failure handling.

## Adding a backend, including GN

1. Implement the four methods in a `StreetRenderBackend` subclass, for example
   under a future `way/street_network/gn_backend/` package.
2. Register a lazy factory at addon startup, before users select the backend.
   The common importer, preparation and GUI selector need no new dispatch branch.
3. Use the generic `app.asset_library.BlenderAssetLoader` and a backend-specific
   catalog for GN programs. Reuse the separate street material library. Keep
   edition/external-library/storage handling in the generic asset service.

Example registration, once the implementation exists:

```python
from .backends import StreetBackendDefinition, registry


def create_gn_backend(app, **options):
    from .gn_backend.renderer import GeometryNodesStreetRenderer
    return GeometryNodesStreetRenderer(app, **options)


registry.register(StreetBackendDefinition(
    key='geometry_nodes',
    label='Geometry Nodes',
    description='Generate streets using the shared GN programs',
    enum_value=3,
    factory=create_gn_backend,
))
```

`GeometryNodesStreetRenderer` in this example must inherit `StreetRenderBackend`
and implement its lifecycle. This is a future registration example; no GN
street-generation implementation is installed by this cleanup. The separate
terrain-projection modifier described above is available. The preserved GN
integration still needs reconciling with the current shared contracts and assets
when it is restored; registration alone does not establish geometry parity.

The registry imports no `bpy` or geometry implementation. Reading backend
metadata or drawing the selector never constructs an engine. Factories import
their implementation only when selected. Keys and numeric enum values must be
unique. Register/unregister providers on Blender's main thread. An optional
backend can call `registry.unregister(key)` when disabled; restore its active
UI selections first. The default Python provider cannot be unregistered.

## GUI and saved settings

The common street physical/material/storage controls stay visible. The backend
selector is hidden when only Python is registered and appears automatically
when a second backend registers. The RNA property `streetBackend` is retained
for saved files, scripts and future selection; its items come from the registry.

A definition may supply `draw_settings(layout, settings)` for additional controls
below the common options. Its property declarations and any specialized editing
panels belong to the backend's own registration; the hook should only draw UI.
The registry retains enum strings because Blender borrows callback strings.

Numeric identities are persistent file data:

- **0**: Python, the default.
- **1**: retired original GN selection; permanently reserved.
- **2**: retired Legacy selection; permanently reserved.
- **3+**: available for new backend definitions; assign stable unique values.

Startup/file loading migrates old GN and Legacy selections to Python. A new GN
provider can reuse the textual key `geometry_nodes` with a new value such as 3;
its active selection survives that migration. Existing generated meshes are
not converted to another backend by changing this import preference.

## Removed Legacy implementation and recovery

Removed `LegacyStreetRenderer`, its helper functions, `way/asset_store.py` and
`way/item_renderer/`. The shared building `app/asset_store.py` and the separate
root `item_renderer/` package were preserved. Historical Legacy-only style/
asset-info requirements and renderer setup were removed as well. Old external
asset packages and untracked prototype Blender files were not modified.

The preceding Blosm commit `5757f06` preserves the Legacy implementation.
For the experimental GN source and fully verified Git-only regeneration, see
[GN recovery](STREETS_GEOMETRY_NODES.md) and the sibling
[long-pause handoff](../../prochitecture_gn/docs/streets/blosm_integration/GN_RESTART.md).
No checkpoint directory is needed.

## Verification

`tests/test_street_backends.py` checks lazy construction, duplicate/reserved IDs,
provider removal, lifecycle/options/results, validation and settings hooks
without Blender. `tests/test_street_backend_registration.py` checks actual
restricted-context addon registration, old-selection migration, preserved
building GN/material controls, and registering/selecting/driving a temporary
GN-style backend through Blender RNA without adding importer branches.

Native preparation, material-only installation, normal real OSM import, exact
previous-mesh comparison and independent reopen/Update are checked separately.
See the current [acceptance record](STREETS_PYTHON_ACCEPTANCE.md) for completed
results and scope. The preceding 932-test GN-removal run remains its own dated
baseline; it is not relabeled as a new full run for this renderer cleanup.

Completed verification: seven host lifecycle/registry tests; four native
registration/extension checks; nineteen preparation checks; four material tests
plus reopening; exact Berlin import/quality and independent edit/Update/reopen.
All 328 imported mesh parts match the preceding Python baseline.
