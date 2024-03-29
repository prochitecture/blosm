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

from . import RoofRealistic
from building.roof.hipped import RoofHipped
from .flat import RoofFlatRealistic, RoofFlat
from .profile import RoofProfileRealistic, RoofProfile


class RoofHippedRealistic(RoofRealistic, RoofHipped):
    
    def renderRoofTextured(self):
        if self.makeFlat:
            return RoofFlatRealistic.renderRoofTextured(self)
        else:
            super().renderRoofTextured()
    
    def renderWalls(self):
        if self.makeFlat:
            if self.mrw:
                RoofFlatRealistic.renderWalls(self)
            else:
                RoofFlat.renderWalls(self)
        else:
            if self.mrw:
                RoofProfileRealistic.renderWalls(self)
            else:
                RoofProfile.renderWalls(self)