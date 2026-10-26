"""Original endpoint aliases on a derived physical incoming profile.

Parallel carriageways sharing a physical mouth are not a source Route. Each
alias retains its original station and maps only its own half of that mouth.
"""
from dataclasses import dataclass
import numpy as np
from .model import End


@dataclass(frozen=True)
class SourceMouth:
    endpoint: End
    station: float
    physical_outer_u: float = 0.
    physical_peak_u: float = 1.
    outward_outer_u: float = 0.
    outward_peak_u: float = 1.

    def outward_parameter(self, physical_u):
        v = ((np.asarray(physical_u)-self.physical_outer_u)
             /(self.physical_peak_u-self.physical_outer_u))
        return self.outward_outer_u+v*(self.outward_peak_u-self.outward_outer_u)

    def contains(self, physical_u, tolerance=1e-10):
        low, high = sorted((self.physical_outer_u, self.physical_peak_u))
        return (np.asarray(physical_u) >= low-tolerance) & (np.asarray(physical_u) <= high+tolerance)


def source_mouths(approach, station):
    derived = getattr(approach, 'source_mouths', None)
    if derived is not None:
        return derived(station)
    if hasattr(approach.path, 'trims'):
        (sid, end), at = list(approach.path.trims(station).items())[-1]
    else:
        sid, end = approach.source.id, approach.end
        at = station if end == 0 else approach.length-station
    return (SourceMouth(End(sid, end), float(at)),)
