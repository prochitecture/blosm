"""Extensible numerical UV strategies for Python-generated street surfaces.

``UVMap`` uses metres in the prepared construction plane. A unit in U or V
therefore represents one metre, independent of scene units, road width, mesh
resolution, display exaggeration or terrain projection. Materials choose the
repeat size. This is plan-distance mapping, not a terrain-surface unwrap.

Strategies return an array per face corner, allowing deliberate texture seams
without splitting geometric vertices. Future providers can implement the same
``coordinates(mesh)`` method for crossing, sidewalk or other local mappings.
"""
from dataclasses import dataclass
from typing import Protocol
import math
import numpy as np


DEFAULT_UV_MAP = 'UVMap'


def publication_uv_coordinates(values):
    """Return finite contiguous float32 UVs, the representation Blender stores.

    UVs use physical metres and receive no scene-unit conversion. Checking the
    cast matters: an extension can produce finite float64 coordinates outside
    Blender's float32 range. Reject those before a replacement mesh is published.
    The returned array can also be used for exact saved-buffer comparisons.
    """
    with np.errstate(over='ignore', invalid='ignore'):
        result = np.ascontiguousarray(values, dtype=np.float32)
    if result.ndim < 1 or result.shape[-1] != 2 or not np.isfinite(result).all():
        raise ValueError('UV coordinates must be finite pairs within Blender float32 range')
    return result


class UVMapping(Protocol):
    def coordinates(self, mesh) -> np.ndarray:
        """Return finite coordinates shaped (number_of_quads, 4, 2)."""
        ...


@dataclass(frozen=True)
class PlanarMapping:
    """Physical XY projection, with an optional origin and frame rotation.

    The default origin belongs to the complete prepared network coordinate
    frame, never to an individual Part's bounding box. Neighboring surfaces
    consequently retain the same texture phase when their Parts are rebuilt.
    Rotation is in radians; it defines the prepared-XY direction of the U axis.
    """
    origin: tuple = (0., 0.)
    rotation: float = 0.

    def coordinates(self, mesh):
        origin = np.asarray(self.origin, dtype=float)
        if origin.shape != (2,) or not np.isfinite(origin).all() or not math.isfinite(self.rotation):
            raise ValueError('Planar UV mapping needs a finite XY origin and rotation')
        c, s = math.cos(self.rotation), math.sin(self.rotation)
        axes = np.array(((c, -s), (s, c)))
        uv = (mesh.xyz[:, :2] - origin) @ axes
        return uv[mesh.quads]


@dataclass(frozen=True)
class RoadMapping:
    """U is longitudinal station; V is signed distance from the source path.

    Supply one station and transverse coordinate for every mesh vertex. The
    transverse values include the Section's lateral offset, so changing width
    does not recenter or rescale a texture. Stations remain physical metres and
    need not follow vertex index order. Optional origins support an explicitly
    coordinated phase between Sections; continuity across arbitrary independent
    Sections is a separate network policy, not inferred from proximity.
    """
    stations: object
    transverse: object
    station_origin: float = 0.
    transverse_origin: float = 0.

    def coordinates(self, mesh):
        stations = np.asarray(self.stations, dtype=float)
        transverse = np.asarray(self.transverse, dtype=float)
        if (stations.shape != (len(mesh.xyz),) or transverse.shape != stations.shape
                or not np.isfinite(stations).all() or not np.isfinite(transverse).all()
                or not math.isfinite(self.station_origin) or not math.isfinite(self.transverse_origin)):
            raise ValueError('Road UV mapping needs finite per-vertex stations and transverse distances')
        uv = np.column_stack((stations - self.station_origin, transverse - self.transverse_origin))
        return uv[mesh.quads]


def apply_mapping(mesh, mapping, name=DEFAULT_UV_MAP, *, overwrite=False):
    """Fill a named UV map while preserving provider maps unless asked to replace.

    Default/fallback mapping can safely run after feature-specific providers:
    it returns the existing payload unchanged if that map already exists.
    """
    if name in mesh.uv_layers and not overwrite:
        return mesh
    return mesh.with_uv_layer(name, mapping.coordinates(mesh))
