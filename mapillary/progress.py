"""Image-work estimates from successfully counted coverage tiles, without bpy."""
import math


def area(bounds):
    west, south, east, north = bounds
    # Spherical rectangle area, omitting the common radius-squared factor.
    return math.radians(east-west)*(math.sin(math.radians(north))-math.sin(math.radians(south)))


def overlap_area(bounds, tile):
    zoom, x, y = tile
    size = 1 << zoom
    west, east = x/size*360.-180., (x+1)/size*360.-180.
    north = math.degrees(math.atan(math.sinh(math.pi*(1-2*y/size))))
    south = math.degrees(math.atan(math.sinh(math.pi*(1-2*(y+1)/size))))
    clipped = (max(bounds[0],west),max(bounds[1],south),min(bounds[2],east),min(bounds[3],north))
    return area(clipped) if clipped[0]<clipped[2] and clipped[1]<clipped[3] else 0.


class ImageProgress:
    def __init__(self, bounds, tile_count):
        self.bounds, self.tile_count = tuple(bounds), tile_count
        self.total_area = area(bounds)
        self.counted_area = 0.
        self.counted_tiles = self.candidates = self.processed = 0

    def counted(self, tile, unique_total):
        # Failed/malformed tiles never reach this method. Empty valid tiles do.
        self.counted_area += overlap_area(self.bounds,tile)
        self.counted_tiles += 1
        self.candidates = unique_total

    @property
    def estimated_total(self):
        if self.counted_tiles == self.tile_count:
            return self.candidates
        if not self.candidates or not self.counted_area:
            return None
        return max(self.candidates, math.ceil(self.candidates*self.total_area/self.counted_area))

    def percent(self, complete=False):
        if complete:
            return 100
        total = self.estimated_total
        return min(99,int(100*self.processed/total)) if total else None
