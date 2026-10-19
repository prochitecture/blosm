"""Resolve numeric carriageway defaults independently of Blender assets."""
import math
from .way_properties import wayCategoryProps


class StreetDimensionPolicy:
    def resolve(self, street, category, oneway):
        defaults = wayCategoryProps.get(category, wayCategoryProps['other_roadway'])
        names = ('totalLanesOneway', 'totalNumLanesOneway', 'totalLanes', 'totalNumLanes') if oneway else ('totalLanes', 'totalNumLanes')
        count = self._style_value(street, names, defaults['nrLanes'])
        width = self._style_value(street, ('laneWidth',), defaults['lane'])
        if not math.isfinite(float(count)) or int(count) != count or count < 1:
            raise ValueError(f'{category}: lane-count default must be a positive integer')
        if not math.isfinite(float(width)) or width <= 0:
            raise ValueError(f'{category}: lane-width default must be positive')
        return {'nrLanes': int(count), 'laneWidth': float(width)}

    @staticmethod
    def _style_value(street, names, default):
        for name in names:
            value = street.getStyleBlockAttr(name)
            if value is not None:
                return value
        return default
