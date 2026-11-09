# Street terrain placement — dev76

**UV integration, dev77:** [UV maps](STREETS_UV.md) describe the prepared
surface and remain unchanged by vertical placement in either terrain mode.
All named face-corner maps survive terrain recovery and republishing. These
are plan-distance coordinates; they are not re-unwrapped over terrain relief.

Implemented on `buildings_gn`, 4 November 2026, for Blender 5.2. The Python
backend generates the street meshes. A separate placement stage offers the
two requested terrain modes; it does not restore the retired GN street solver.

## Using it

For a new OSM import, select the usual **Terrain** object in Blosm's import
panel. Under **Import roads and paths**, choose **Terrain placement**:

| Mode | Result | Update behavior |
| --- | --- | --- |
| **Preserve profiles (Python)** — default, proposal 3 | Fits a derived longitudinal elevation reference; translates complete cross-sections and coordinated junctions | Edit source controls, widths, terrain or placement settings, then press **Update** |
| **Direct projection (GN)** — proposal 1 | Projects each vertex vertically and restores its prepared relative height | Press **Update** to install/change the mode and validate coverage; the installed modifier subsequently evaluates terrain changes live |

For an existing generated network, select its controller, source Street, or
generated mesh. In the 3D View's **Streets → Python streets** panel, set
**Original terrain**, choose the mode, then press **Update**. This setting is
per connected network. Clearing Original terrain and updating restores the
prepared mesh elevation. Changing the import preference does not modify
previous imports. No terrain target means the existing construction behavior.

| Setting | Default | Meaning |
| --- | --- | --- |
| Clearance | 0.02 m | Vertical clearance above sampled original ground |
| Maximum reference grade | 0.10 = 10% | Maximum rise/run of the **added longitudinal reference**, not a bound on existing transverse crown or junction relief |
| Reference spacing | 4 m | Approximate distance between reference controls; seam contacts add compulsory controls |
| Terrain sample spacing | 1 m | Maximum requested spacing along a quad edge for terrain checks; every quad receives at least a 3 × 3 sampling |
| Maximum lift above ground | 25 m | Reject a fit whose vertex surface-to-ground gap exceeds this limit plus clearance |

References never introduce unconstrained controls finer than the published
road rows: a finer requested spacing uses the actual rows instead. This avoids
between-row dips and bounds control allocation by the existing mesh.

Source controls remain in their construction XY plane. Terrain placement
derives their vertical alignment; moving a source vertex's Z is not a manual
elevation override. Controller translation and Z rotation are supported;
controller scale/tilt is rejected. The terrain may be translated, rotated,
scaled and have its own independent modifiers. Scene unit scale is respected.

## What proposal 3 preserves

Ordinary carriageways retain every full cross-section, including its crown,
offset and prescribed profile rise. Their XY positions, rounded outlines,
quad connectivity, resolution and materials do not change. A junction's
existing curved surface receives one uniform vertical translation. All
connected mouths share its minimum, including the low rounded corners.
Each road's own crown rise remains intact. The additional reference has zero
derivative at a junction mouth; its analytical join is C1. The finite mesh
approximates that reference at its existing rows.

**Initial Bundle policy:** a complete Bundle region and any junctions sharing
its seams use one rigid elevation. This preserves divided profiles, separator
plateaus and overlap constructions. It is conservative on sloping or long
Bundles: it can leave substantial gaps above lower terrain, and the maximum
lift setting can reject it. Smoothly varying, coupled longitudinal references
inside Bundles are a later refinement, not implemented here. Ordinary roads
joining these protected regions use the fitted reference as usual.

The original terrain is never edited. No shoulders, cut/fill or earthworks are
generated. These belong to proposal 4. Sidewalk generation is also not included
yet; its footprint and height integration contract is implemented below.

## Complete corridor contract

The numerical contract is in
[terrain/alignment.py](../way/street_network/python_backend/terrain/alignment.py),
with the current mesh adapter in
[terrain/corridor.py](../way/street_network/python_backend/terrain/corridor.py).
It has no Blender or SciPy dependency; it uses NumPy.

- `ElevationReference(identity, stations)` names ordered longitudinal controls.
  A single station represents a rigid region.
- `CorridorComponent(identity, reference, xyz, quads, stations, datum, clearance, contacts)`
  supplies the **entire component footprint**, its
  already-designed heights, a reference station for every vertex, and any
  extra clearance. Contacts explicitly identify shared reference heights.
- `PlacementSettings` supplies fitting and sampling policy.
- `fit(references, components, terrain_sampler, settings)` returns references,
  placed coordinates and a report. `terrain_sampler(N×2 XY)` returns terrain Z,
  or raises on any missing sample.

All components using a reference constrain it together. A future curb or
raised sidewalk does not become the carriageway datum, and a material boundary
does not create a new height reference. Width, pavement materials and local
curb drops remain feature-specific geometry inputs.

A feature provider can publish a normal `Part` with this metadata:

```python
metadata['terrain_corridor'] = {
    'section_id': section.id,
    'stations': station_for_each_component_vertex,
    'datum': section.datum,
    'clearance': 0.0,  # additional to the network clearance
}
```

This binds a side component to the road's existing reference. Its prepared Z
already includes curb rise, crossfall or other prescribed relative heights.
Split such a component at junction-owned spans: one component must use one
surviving road reference. Junction components can instead share the rigid
junction reference through the numerical contract. A synthetic raised side
strip regression demonstrates that a component outside the carriageway can
raise the same fitted reference without changing the crown.

This keeps the future proposal-4 contact boundary at the outer edge of the
**complete corridor**, including sidewalks, rather than permanently assuming
that terrain must meet the carriageway edge.

## Construction and formulas

Let $z_0(v)$ be a prepared vertex height, $d$ its construction datum, $p(v)$
its relative height, and $T(x,y)$ the original terrain height:

$$
p(v) = z_0(v) - d.
$$

Direct projection evaluates

$$
z_1(v) = T(x_v, y_v) + p(v) + c.
$$

It preserves additive relief, **not** the complete cross-section on transverse
terrain slopes. Proposal 3 instead uses a common baseline $b(s)$:

$$
z_3(v) = b(s_v) + p(v).
$$

The full corridor provides sampled lower bounds. For each face $f$, use the
lower prepared interpolation from either possible quad diagonal, denoted
$p^-(x,y)$, and obtain

$$
F_f = \max_{(x,y)\in S_f} \left[T(x,y) + c + c_f - p^-(x,y)\right].
$$

Every reference control influencing a vertex of that face receives the bound
$F_f$. This deliberately conservative rule constrains the published faces as
well as the analytical reference. Contacts merge controls into shared graph
nodes. Let $L_i$ be their accumulated lower bounds. For adjacent controls a
distance $\Delta s$ apart, the graph majorant enforces

$$
\begin{aligned}
b_i &\geq L_i, \\
\lvert b_{i+1} - b_i \rvert &\leq \frac{1}{2}\,g_{\max}\,\Delta s.
\end{aligned}
$$

A maximum-priority queue propagates necessary raises. There is no iterative
optimization, and no smoothing followed by clipping. The factor one-half
reserves slope for the interpolating Hermite curve. Endpoint derivatives are
zero. Interior derivatives use a weighted harmonic mean of same-sign secants,
capped at twice the smaller neighboring secant magnitude; opposite signs give
zero derivative. This gives monotone segments without height overshoot.

For interval length $h$, normalized position $t$, heights $b_0, b_1$ and
derivatives $m_0, m_1$:

$$
\begin{aligned}
b(t) ={}& (2t^3 - 3t^2 + 1)b_0 + (t^3 - 2t^2 + t)h m_0 \\
        &+ (-2t^3 + 3t^2)b_1 + (t^3 - t^2)h m_1.
\end{aligned}
$$

The implementation checks the exact extrema of its quadratic derivative,
including interior stationary points, against the requested maximum grade.
It also checks lift above ground. A rigid junction uses constant $b = Z_j$.
An infeasible lift or missing terrain sample rejects the Update and leaves
previous generated objects and the committed manifest intact.

The guarantee is **sampled clearance**, not continuous clearance against every
terrain triangle. Narrow spikes/holes between samples can be missed. Decrease
sample spacing for rough terrain. Sampling is bounded to 64 subdivisions per
quad; exceeding this produces a diagnostic instead of an unbounded allocation.
This changes validation sampling, not road mesh resolution.

## Update, persistence and dependencies

[terrain/blender.py](../way/street_network/python_backend/terrain/blender.py)
samples an evaluated copy of the original terrain in the controller's physical
coordinate frame. Its geometry fingerprint, placement settings, scale and
placement revision invalidate the fit. A source edit uses the existing
numerical dependency graph, then recomputes the connected network's placement.
Only changed output coordinates, newly rebuilt parts or changed placement
settings cause publication. This initial version resamples the network rather
than maintaining a spatial terrain cache. Unrelated network controllers remain
independent.

`st_terrain_prepared_z` stores the prepared physical Z on each mesh point.
Furniture stores its prepared position. Switching modes, clearing terrain,
reopening and repeated Updates therefore start from the original relief,
never the previous projected height. Furniture receives the corresponding
street/junction reference displacement. No viewport polling or depsgraph edit
handler is installed. Generated mesh buffers are released from Python after
publication and recovered transiently when placement needs them.

Saved settings are in `Object.blosm_street_terrain`; the controller retains
`st_py_terrain_mode`, `st_py_terrain_signature`, and `st_py_terrain_report`.
The raw numerical update is committed separately from the placed publication,
so terrain does not contaminate geometry caches. Errors use the existing
outdated-output notice and retain the last valid objects.

## Direct GN asset and regeneration

Asset: [assets/prochitecture_street_terrain.blend](../assets/prochitecture_street_terrain.blend).
It contains one shared 21-node group,
`Prochitecture Streets / Direct terrain projection`, revision `1.0`, and
**no materials or objects**. Mesh materials continue to come from
`prochitecture_street_materials.blend` through the existing material service.
The generic `BlenderAssetLoader` loads the projection group from the base asset
catalog. Scene modifiers reference this shared group rather than per-fixture
copies. The library is appended, so a saved scene evaluates without the
generator checkout or a script execution.

The readable generator and exact command are in the sibling repository's
[terrain placement experiment](../../prochitecture_gn/experiments/streets/terrain_placement/README.md).
It does not use any checkpoint or legacy GN street library.

The graph captures `Position.Z − Prepared Datum` before projection. Terrain
Object Info uses relative coordinates. Rays start above the terrain bounding
box and reach below its minimum. Set Position restores relative height plus
clearance exactly once. Python preflights complete footprint coverage on
Update. Afterwards a new GN ray miss retains that vertex's prepared position
and records `st_terrain_hit = false`; press Update to get the explicit coverage
error. The simple modifier does not perform Python's fitting/grade checks.

## Inspection and validation

The reproducible gallery, real OSM scene links, commands and measured results
are recorded in the sibling [terrain placement README](../../prochitecture_gn/experiments/streets/terrain_placement/README.md).
The illustrated [concept guide](../../prochitecture_gn/docs/streets/terrain_placement/README.md)
remains useful; proposal 4 in that guide is future work.

Tests include complete profiles, shared junction heights, unchanged topology,
separated Bundle profiles, a future raised side strip, missing terrain, lift
rejection, mode cycling, original-terrain retention, source/width edits,
furniture references, no-op output identity, rollback, reopening, nonzero
construction datum, centimetre scenes, transformed terrain and rotated
controllers. Karl-Marx-Allee is tested as a real OSM network over
**synthetic rolling terrain**, not a surveyed Berlin DEM. The implementation does not
depend on that file.
