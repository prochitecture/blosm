"""Bounded regular arc rows that retain mandatory physical curve corners."""
import math
import numpy as np
from .bundle_terminal_caps import sample_chain


def landmark_indices(points, landmarks):
    points, landmarks = np.asarray(points), np.asarray(landmarks).reshape(-1,2)
    result = set()
    for landmark in landmarks:
        distances = np.linalg.norm(points-landmark, axis=1)
        at = int(np.argmin(distances))
        if distances[at] <= 1e-6:
            result.add(at)
    return result


class ProtectedPolylineSampler:
    def __init__(self, *, spacing, tolerance, max_intervals=4096):
        if (not math.isfinite(spacing) or spacing<=0. or not math.isfinite(tolerance)
                or tolerance<=0. or type(max_intervals) is not int or max_intervals<1):
            raise ValueError('Protected curve sampling needs finite positive limits')
        self.spacing, self.tolerance, self.max_intervals = spacing, tolerance, max_intervals

    def sample(self, points, protected=()):
        points = np.asarray(points,dtype=float)
        if points.ndim!=2 or points.shape[1]!=2 or len(points)<2 or not np.all(np.isfinite(points)):
            raise ValueError('Protected curve sampling needs a finite physical polyline')
        cuts = sorted({0,len(points)-1,*protected})
        if cuts[0]!=0 or cuts[-1]!=len(points)-1:
            raise ValueError('Protected curve sample is outside its source')
        result, owners, maximum = [], [], 0.
        for first,last in zip(cuts,cuts[1:]):
            part = points[first:last+1]
            distances = np.linalg.norm(np.diff(part,axis=0),axis=1)
            if np.any(distances<=0.):
                raise ValueError('Protected curve has a collapsed physical interval')
            stations = np.r_[0.,np.cumsum(distances)]
            count = max(1,int(math.ceil(stations[-1]/self.spacing)))
            for _ in range(16):
                if len(owners)+count>self.max_intervals:
                    raise ValueError('Protected curve sampling exceeds its interval budget')
                try:
                    sampled = sample_chain(part,count,self.tolerance)
                    break
                except ValueError:
                    count = max(count+1,int(math.ceil(1.25*count)))
            else:
                raise ValueError('Protected curve has no bounded-error regular sampling')
            # Check the complete original/replacement station partitions, not
            # only the selected samples. Linear interpolation bounds each span.
            times = np.linspace(0.,stations[-1],count+1)
            query = np.unique(np.r_[times,stations])
            original = np.column_stack([np.interp(query,stations,part[:,i]) for i in range(2)])
            actual = np.column_stack([np.interp(query,times,sampled[:,i]) for i in range(2)])
            maximum = max(maximum,float(np.linalg.norm(original-actual,axis=1).max()))
            if maximum>self.tolerance:
                raise ValueError('Protected rows exceed the complete curve-error budget')
            result.extend(sampled[:-1])
            owners.extend([first]*count)
        result.append(points[-1])
        return np.asarray(result), owners, dict(intervals=len(owners), maximum_station_error=maximum,
            tolerance=self.tolerance, mandatory_corners=len(cuts), interval_budget=self.max_intervals)
