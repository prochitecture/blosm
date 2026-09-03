# Transforms

This document describes how transformations from the OGC 3D Tiles specifications 1.0 and 1.1 are implemented in the Blosm code.

According to section 6.7.5.2 (glTF transforms) of the OGC 3D Tiles specification 1.0, the order of transformations is:
1. glTF node hierarchy transformations
2. glTF y-up to 3D Tiles z-up transform
3. Any tile format specific transforms.
    * A Batched 3D Model Feature Table may define `RTC_CENTER`, which translates model vertices.
    * An Instanced 3D Model Feature Table defines per-instance position, normals, and scales. These are used to create per-instance 4x4 affine transform matrices that are applied to each instance.
4. Tile transform

The order of transformations is the same in the OGC 3D Tiles specification 1.1, except that point 3 is not applicable there.

Below, assume:

* `C`: `centerCoords` used in the Blosm code. This is a WGS84 ECEF coordinate in the global 3D Tiles coordinate system.
* `M`: accumulated 3D Tiles tile `transformMatrix` used in the Blosm code.
* `R`: `RTC_CENTER`. This is in the post-glTF, post-y-up-to-z-up tile/content coordinate system. It is applied before `M`, so it is not generally a global/ECEF coordinate if `M` is not the identity matrix.
* `G`: imported glTF object transform / node hierarchy after Blender's glTF axis conversion.
* `v`: glTF mesh vertex.
* `T(X)`: translation matrix by vector `X`.

With column-vector notation, the 3D Tiles placement before Blosm's scene offset is:

```
M @ G @ v
```

without `RTC_CENTER`, and:

```
M @ T(R) @ G @ v
```

with `RTC_CENTER`.

## RTC_CENTER is not available

Blosm keeps the glTF importer offset patch enabled. The offset `-C` is used to reduce numerical errors caused by large coordinate values. The imported glTF transform is effectively:

```
T(-C) @ G
```

According to the specification, without the Blosm offset, the 3D Tiles transform would be:

```
M @ G
```

Since the offset `-C` was already applied by the importer, a compensation wrapper must be applied to `M`:

```
T(-C) @ M @ T(C)
```

This temporarily moves coordinates back by `+C`, applies `M`, then shifts the coordinates back by `-C`. This is done in `apply_transform_matrix(..)` in `threed_tiles/blender.py`.

The resulting transformation matrix:

```
T(-C) @ M @ T(C) @ T(-C) @ G = T(-C) @ M @ G
```

## RTC_CENTER is available

Blosm temporarily disables the glTF importer offset patch when importing B3DM content with `RTC_CENTER`. The imported glTF transform is therefore:

```
G
```
The Blosm scene-center offset `-C` is still applied, but this time in the method `process_rtc(..)` in `threed_tiles/blender.py`. The RTC translation `R` is also applied there:

```
T(R) @ T(-C) @ G = T(R-C) @ G
```

In the current code this is written directly to `obj.matrix_world`:

```
obj.matrix_world = Matrix.Translation(rtc_center - self.centerCoords) @ obj.matrix_world
```

Using matrix operations is important here because changing `obj.location` does not immediately refresh `obj.matrix_world`; direct matrix assignment avoids requiring `bpy.context.view_layer.update()` before `apply_transform_matrix(..)` reads `obj.matrix_world`.

As in the no-RTC case, `apply_transform_matrix(..)` applies the compensation wrapper around `M`:

```
T(-C) @ M @ T(C)
```

The resulting transformation matrix:

```
T(-C) @ M @ T(C) @ T(R - C) @ G = T(-C) @ M @ T(R) @ G
```

This matches the 3D Tiles order: glTF transforms, then `RTC_CENTER`, then the tile transform, then Blosm's global-to-local offset.

## Final Blender orientation

After all selected 3D Tiles content has been imported, `finalize(..)` rotates the ECEF-centered result into Blosm's local Blender coordinate system. The rotation matrix is:

```
Q = Matrix.Rotation(lat - pi/2, 4, 'X') @ Matrix.Rotation(radians(-90 - centerLon), 4, 'Z')
```

where `lat` is computed from `centerCoords`.

For joined 3D Tiles objects, Blosm joins the imported objects, sets the object origin, assigns:

```
matrix_local = Q
```

and then applies rotation to the mesh data. For non-joined objects, each object location is first transformed by `Q`, and then the same rotation is assigned per object.

