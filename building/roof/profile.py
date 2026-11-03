"""
This file is a part of Blosm addon for Blender.
Copyright (C) 2014-2018 Vladimir Elistratov
prokitektura+support@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import math
from . import Roof
from ...padl.profile_geometry import ProfiledRoofGeometry, ProfiledVert, Slot
from util.osm import parseNumber


# Use https://raw.githubusercontent.com/wiki/vvoovv/blosm/assets/roof_profiles.blend
# to generate values for a specific profile
gabledRoof = (
    (
        (0., 0.),
        (0.5, 1.),
        (1., 0.)
    ),
    {
        "numSamples": 10,
        "angleToHeight": 0.5
    }
)

roundRoof = (
    (
        (0., 0.),
        (0.01, 0.195),
        (0.038, 0.383),
        (0.084, 0.556),
        (0.146, 0.707),
        (0.222, 0.831),
        (0.309, 0.924),
        (0.402, 0.981),
        (0.5, 1.),
        (0.598, 0.981),
        (0.691, 0.924),
        (0.778, 0.831),
        (0.854, 0.707),
        (0.916, 0.556),
        (0.962, 0.383),
        (0.99, 0.195),
        (1., 0.)
    ),
    {
        "numSamples": 1000,
        "angleToHeight": None
    }
)

gambrelRoof = (
    (
        (0., 0.),
        (0.2, 0.6),
        (0.5, 1.),
        (0.8, 0.6),
        (1., 0.)
    ),
    {
        "numSamples": 10,
        "angleToHeight": None
    }
)

saltboxRoof = (
    (
        (0., 0.),
        (0.35, 1.),
        (0.65, 1.),
        (1., 0.)
    ),
    {
        "numSamples": 100,
        "angleToHeight": 0.35
    }
)


class RoofProfile(ProfiledRoofGeometry, Roof):
    """Legacy OSM/renderer adapter around the shared profile geometry kernel."""
    defaultHeight = 3.

    def __init__(self, data):
        Roof.__init__(self)
        ProfiledRoofGeometry.__init__(self, data)

    def init(self, element, data, osm, app):
        self.projections.clear()
        super().init(element, data, osm, app)
        self.initProfile()


    def getRoofHeight(self):
        tags = self.element.tags

        h = parseNumber(tags["roof:height"]) if "roof:height" in tags else None
        if h is None:
            if not self.angleToHeight is None and "roof:angle" in tags:
                angle = parseNumber(tags["roof:angle"])
                if not angle is None:
                    self.processDirection()
                    h = self.angleToHeight * self.polygonWidth * math.tan(math.radians(angle))
            if h is None:
                # get the number of levels
                if "roof:levels" in tags:
                    h = parseNumber(tags["roof:levels"])
                h = h * self.levelHeight if h else self.defaultHeight
        return h
