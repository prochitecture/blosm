"""Exercise registered native street and simple-building producers together.

Loaded through the public setup-script hook by check_street_progressive.py.
The UI currently selects these renderers in different modes; this fixture uses
their real managers, parsers, assets and backends without changing UI policy.
"""
from . import Setup
from ..action.generate_streets import StreetGenerator
from ..way.renderer_streets import StreetRenderer
from ..building.padl_integration.manager import PadlBuildingManager


def setup(app, osm):
    source = Setup(app, osm)
    source.skipWays()
    source.roadsAndPaths()
    ways = source.getWayManager()
    ways.addAction(StreetGenerator(ways.styleStore, getStyle=ways.getStyle,
                                  doPixelCorrection=False))
    ways.addRenderer(StreetRenderer(app))
    buildings = PadlBuildingManager(osm, app)
    osm.addCondition(buildings.accepts, None, buildings)
