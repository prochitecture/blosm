"""Typed geometry-only junction areas, separate from original incidences."""


class GeneratedJunctionArea:
    """A local surface construction must not create original End contacts."""
    __slots__=()


def is_generated_junction_area(area):
    return isinstance(getattr(area.layout,'generated_junction_area',None),GeneratedJunctionArea)
