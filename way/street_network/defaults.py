"""Numeric preparation styles available without an external asset package.

The styles do not select meshes, shaders or renderers. Category-specific lane
counts/widths still resolve through StreetDimensionPolicy and the source tags.
"""
from ...grammar.street import Street, Section, SideLane
from ...setup.realistic_streets import _categoryToStreetStyle


def styles():
    return {name: [Street(), Section(), SideLane(length=20.)]
            for name in set(_categoryToStreetStyle.values()) | {'residential'}}
