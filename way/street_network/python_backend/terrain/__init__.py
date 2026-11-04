"""Terrain placement contracts and numerical fitting, independent of Blender.

Features contribute their entire footprint, prescribed relative heights and
clearance. A reference is shared by every component of a cross-section. A
one-station reference denotes a rigid junction or protected Bundle region.
The original terrain is read only; terrain earthworks are a separate feature.
"""
from .alignment import CorridorComponent, ElevationReference, PlacementSettings, fit

__all__ = ('CorridorComponent', 'ElevationReference', 'PlacementSettings', 'fit')
