"""Carriageway design policy shared by Blender integration and offline readers."""
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class StreetDesign:
    profile_rise: float = .08
    cell_size: float = .7
    interpolation: int = 3
    resolution: int = 1
    exaggeration: float = 10.
    categories: tuple = ('motorway', 'motorway_link', 'trunk', 'trunk_link',
        'primary', 'primary_link', 'secondary', 'secondary_link', 'tertiary',
        'tertiary_link', 'residential', 'living_street', 'service', 'unclassified')

    def __post_init__(self):
        if not math.isfinite(self.profile_rise) or self.profile_rise < 0:
            raise ValueError('Profile rise must be finite and nonnegative')
        if not math.isfinite(self.cell_size) or self.cell_size <= 0:
            raise ValueError('Cell size must be positive')
        if self.interpolation not in (0, 3) or self.resolution not in (1, 2, 3, 4):
            raise ValueError('Unsupported interpolation or mesh resolution')
        if not math.isfinite(self.exaggeration) or self.exaggeration < 0:
            raise ValueError('Exaggeration must be finite and nonnegative')

    @classmethod
    def from_app(cls, app):
        return cls(profile_rise=getattr(app, 'streetProfileRise', .08),
                   cell_size=getattr(app, 'streetCellSize', .7),
                   resolution=getattr(app, 'streetResolution', 1),
                   exaggeration=getattr(app, 'streetExaggeration', 10.))
