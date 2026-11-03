# PADL simple 3D integration

In **3D simple**, choose the building backend **Geometry Nodes** (default) or
**Python**, and a building seed. PADL owns default roof/level/height policies.
**Single object** packs footprints into one input; disabling it creates dedicated
objects per complete building, keeping its outline, parts and holes together.
GN attaches the shared precompiled program. Python retains hidden
editable inputs and creates separate volume meshes. Realistic mode keeps its
existing implementation until its PADL operations/adapters are available.

## Assets

`assets/prochitecture_buildings.runtime.json` maps the role `simple_3d` to the
current `PADL simple_3d` tree in `prochitecture_buildings.blend`. Import code
does not contain the tree name. Rename/rebuild it and republish the catalog.
Loading checks library bytes, recursive graph fingerprints and program identity;
verified dependencies are reused across objects, service instances and reopen.
Edited/conflicting assets are refused.
One import batch verifies the shared program once and reuses a temporary status
reader across its building objects. The reader is removed when the batch ends;
later imports/independent loads perform their own verification.

`prochitecture_buildings.simple_3d.json` is extracted from the group's saved PADL
Text. Python reads it and verifies the same program hash without loading GN.
The compiler's roof-provider `.assets.json` catalog remains separate from the
runtime `.runtime.json` role catalog. Optional program `edition: "PRO"` restricts
that role in both backends. Base/Pro roof allocation remains a build-profile
choice. Publication derives `available_roofs` from the matching compiler
catalog and actual program dependencies; Python uses those same capabilities.
Ship the Blender file, catalogs and extracted document together.

Publish runtime sidecars after a build (Windows cmd, one line):

```bat
"C:\Program Files\Blender Foundation\Blender 5.2\blender.exe" --background --factory-startup "C:\dev\workspace\blosm\assets\prochitecture_buildings.blend" --python-exit-code 1 --python "C:\dev\workspace\blosm\tools\publish_building_assets.py" -- --tree "PADL simple_3d" --role simple_3d
```

Publication leaves the Blender library unchanged. Installed Blosm requires no
PADL, sample or GN source checkout paths.

## Mapping policy

`OsmFaceAttributeMapper` runs replaceable, ordered `TagRule` objects. Its sparse
logical map goes through shared `FootprintBindings`, writing FACE INT categories,
FLOAT heights and BOOLEAN presence flags. Add rules without changing backends.

| OSM tags | Footprint inputs / initial policy |
| --- | --- |
| `building:use` | Explicit present use takes precedence; unrecognized uses become `unknown` |
| `building=apartments` | `building_use=apartments` |
| `building=house/detached/semidetached_house/terrace/bungalow` | `building_use=house`; terraces are not subdivided |
| `building=office` | `building_use=office`, unless explicit use overrides it |
| `building=yes`, nonempty/non-`no`/non-`vacant` `office=*` | `office`, unless `amenity` supplies another function; an office POI does not reclassify an apartment block |
| Other `building`, `amenity`, `man_made` functions | `building_use=unknown`; extension rules, binding codes and PADL branches can add categories |
| `height` | `max_elevation` plus its presence flag; total height including roof, relative to T_min |
| `min_height` | `min_elevation`; absent becomes zero, without a flag |
| `roof:height` | `roof_height` plus its flag; already included in total height; ignored for explicitly flat roofs |
| `building:levels` | Positive integral `levels_wall_count` plus its flag; includes ground floor, excludes roof/basement floors |
| `roof:shape` | Flat/skillion/gabled/hipped; absent uses PADL default; unsupported supplied shapes are reported and omitted |

The current [OSM height definition](https://wiki.openstreetmap.org/wiki/Key:height)
measures from the lowest terrain contact point, matching PADL's T_min. Do not
add terrain relief or roof height. Lengths accept metres and feet/inches.
Invalid, negative, nonfinite, ranged or multi-valued heights and fractional counts
are reported and left absent. Explicit zero roof height stays supplied.

`ele`, `level`, generic `levels`, roof storeys and basement storeys are not
substitutes. Nonzero `building:min_level` needs a future part/raised-floor adapter.
Building parts generate independently with their own dimensions and holes.
Non-building amenity/office areas and `building=no/demolished/razed/construction/proposed`
outlines are excluded. Classification priorities and error handling are Blosm
policy rather than OSM requirements.

Initial import leaves `base_elevation`, `wall_min_height`, `top_height`,
`level_height`, `ground_floor_height` and `ground_floor_factor` absent. PADL
definitions/formulas supply them. Inferred, random and fitted outputs never
become supplied inputs automatically.

Official sources checked on 2026-10-31:
[building use](https://wiki.openstreetmap.org/wiki/Key:building:use),
[house types](https://wiki.openstreetmap.org/wiki/Tag:building%3Dhouse),
[office](https://wiki.openstreetmap.org/wiki/Key:office),
[height](https://wiki.openstreetmap.org/wiki/Key:height),
[min_height](https://wiki.openstreetmap.org/wiki/Key:min_height),
[roof height](https://wiki.openstreetmap.org/wiki/Key:roof:height),
[wall levels](https://wiki.openstreetmap.org/wiki/Key:building:levels),
[skipped levels](https://wiki.openstreetmap.org/wiki/Key:building:min_level).

## Architecture and edits

`PadlBuildingManager` connects the import lifecycle to `OsmFootprintSource` and
`BuildingGenerationService`. The source handles ways, multi-outer relations,
courtyard ownership and duplicate member suppression, preserving corner order
and winding. Full OSM way/relation keys determine samples independently of the
integer face association IDs.

Building outlines and parts use the [v2 footprint contract](../../padl/docs/FOOTPRINT_CONTRACT.md):
type 1 is a standalone building, type 3 a separate part, type 5 a building with
parts whose outline is context only, and type 6 an outline that also generates
its own part volume. Types 5 and 6 scan for type-3 faces with matching `parent_id`.
Type 6 generates once even when the scan is empty. No `has_parts`/`also_part`
attributes are stored. `has_holes` remains independent on outlines and parts.
Holes reference their immediate owning footprint, so they do not cut other parts.

Explicit `type=building` relation ownership takes precedence; otherwise a part
must be fully contained in exactly one outline (shared edges are allowed).
Ambiguous or orphan parts are reported and omitted. Distinct source elements
with coincident boundaries retain separate IDs; one element tagged both building
and part uses type 6. Overall building heights and level counts are not inherited
by parts. A generic part can inherit the building's use; an explicit part use
or `building:use` takes precedence. Unsupported child geometry never causes a
whole-building fallback volume. Volumes are joined without Boolean union or
automatic removal of internal/shared walls.

**Every part uses the whole building's T_min**, sampled from the building's outer
and hole boundaries. Its elevation inputs are offsets from that common datum,
including raised `min_height`. Its own lower terrain samples do not shift it.
Strict terrain coverage is checked for the reference rings and each part's own
rings. With no terrain, the common datum is world Z=0. GN compiler version 2 is
required for these roles and placement rules; rebuild older compiled programs.

`BuildingBackend` has GN and Python strategies consuming the same document,
inputs, identities and seed. `blosm/padl` remains a synchronized shared-core and
Python bundle. Blosm owns the precompiled-GN runtime adapter and verified asset
service; compilation, providers and source builders stay outside the addon.
Future realistic work can add program roles, mapping rules and operation/geometry
adapters while retaining the parser and import service.

Assigned terrain must cover every outer/hole boundary vertex. Both backends use
Blender's evaluated terrain triangulation and T_min; no terrain uses world Z=0.
Python gabled roofs support simple convex and concave footprints without holes.
They reuse the existing profile geometry algorithm, extracted into the synchronized
`padl/profile_geometry.py`; legacy `RoofProfile` uses that same kernel. PADL
provides dimensions, direction and placement explicitly, preserving source corners.
GN gabled roofs currently require convex, hole-free inputs. Geometry failures are reported
per footprint without silently choosing a different roof.

Flat roofs without holes are native single n-gons in both backends, built
directly from their footprint boundary, including concave/collinear corners.
Courtyard roofs keep their multi-face fill; no roof triangles are merged to
construct the hole-free cap. Source footprints remain editable.

Use the bundled API to set values and presence together:

```python
from blosm.padl.api import edit_parameter
edit_parameter(owner, footprint_id, 'levels_wall_count', 6)
```

GN reevaluates live, including branches. Supplied `max_elevation` remains fixed
and floor heights are fitted. Release H only by explicit caller/user action:
`edit_parameter(owner, footprint_id, 'max_elevation', clear=True)`.
For Python, replace the output mesh while retaining the editable owner:

```python
from blosm.building.padl_integration.backends import BuildingGenerationService
service = BuildingGenerationService(backend='python')
result = service.regenerate(owner, collection=output.users_collection[0], terrain=terrain, output_object=output)
```

Unhide Python input objects in the Outliner to inspect them. Each owner references
a Blender Text `PADL building import report` with mapping notices and backend
statuses. Failed footprints remain editable. Native GN warnings may additionally
describe rejected footprints; valid-only inputs must produce no warnings.

## Berlin acceptance

`tests/blender_padl_berlin.py` runs the actual addon import against
`C:\prochitecture\setup\bldg2_berlin_karl_marx_allee.blend`, its `Terrain` and
the local OSM extract. Copies go under `C:\dev\workspace\artifacts\padl_berlin`;
the original file hash is checked unchanged. `blender_padl_reopen.py` reopens
both outputs, checks reuse/parity and edits counts with H fixed and unconstrained.

Building-part acceptance on 2026-11-03 retains 188 footprints: 133 standalone
buildings, five context-only outlines and 50 associated parts. Both backends
generate 168 GN volumes (45 parts) and 183 Python volumes (50 parts).
Python now generates all fifteen concave gabled volumes rejected by GN, including
the main body and tower of St. Pius-Kirche. Mapping reports one unresolved part
association and twelve unsupported roof-shape notices. Compare vertical bounds
on their common generated footprint IDs. The Python run has zero generation
failures. GN has 1,965 faces; Python has 2,201 faces. These are
explicit current profile/solver limits, not complete source coverage.

Run either backend with `--dedicated` to test Single object disabled. Outputs
for the part acceptance runs are under `artifacts/padl_parts/packed` and
`artifacts/padl_parts/dedicated`. `tests/blender_padl_parts.py` additionally covers
type-6 outlines (with and without additional parts), separate courtyard ownership,
raised parts on sloping terrain, failed reference/part terrain coverage, input
edits, invalid ownership, packing equivalence and save/reopen.
