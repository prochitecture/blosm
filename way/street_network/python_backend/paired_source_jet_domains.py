"""Declared one-sided chord domains of an actual paired incoming height field.

These are prediction pieces, not new physical boundaries or source intervals.
The default height field, source curves and sampling remain unchanged.
"""
from dataclasses import dataclass
import numpy as np
from .model import Key
from .source_identity import matches_current_section

SOURCE_PAIRED_JET_POLICY = '1-actual-chord-prediction-domains'


@dataclass(frozen=True)
class PairedSourceJetDomain:
    approach: object

    @classmethod
    def create(cls, context, approach):
        from .bundle_merged_entrances import PairedAngularEntrance
        from .bundle_angular_curves import AngularCrown, RayCurb
        if type(approach) is not PairedAngularEntrance or type(approach.field) is not AngularCrown:
            return None
        if (len(approach.paths) != 2 or len(approach.visits) != 2
                or len(approach.original_ends) != 2 or len(approach.field.curves) != 4
                or len(approach.run.members) != 2 or approach.end_index not in (0, 1)
                or approach.field.frame is not approach.frame_chart
                or approach.sign != (1. if approach.end_index == 0 else -1.)
                or not np.isfinite((approach.length, approach.origin_station)).all()
                or approach.length <= 0.):
            raise ValueError('Paired source jet lost its original finite incoming field')
        for role, (path, visit, end) in enumerate(zip(approach.paths, approach.visits, approach.original_ends)):
            profile = context[Key('profile', visit.section_id)].oriented(int(not visit.forward))
            member = approach.run.members[role]
            if (path is not context[Key('path', visit.section_id)]
                    or not matches_current_section(context, path.section)
                    or not member.visits or visit != member.visits[0 if approach.end_index == 0 else -1]
                    or end != member.ends[approach.end_index]
                    or end.end != (int(not visit.forward) if approach.end_index == 0 else int(visit.forward))
                    or end.section_id != visit.section_id
                    or profile.half != (1 if role == 0 else -1)
                    or profile.rise != approach.field.rise or profile.datum != approach.field.datum):
                raise ValueError('Paired source jet changed an original Path or shared profile')
        for curve in approach.field.curves:
            if (type(curve) is not RayCurb or curve.frame is not approach.frame_chart
                    or curve.stations.ndim != 1 or curve.points.shape != (len(curve.stations), 2)
                    or len(curve.stations) < 2 or len(curve.stations) > approach.maximum_samples
                    or not np.isfinite(curve.points).all() or not np.isfinite(curve.stations).all()
                    or np.any(np.diff(curve.stations) <= 1e-9)):
                raise ValueError('Paired source jet needs its declared finite original chord lattice')
        return cls(approach)

    def outward_knots(self):
        a = self.approach
        chart = np.concatenate([curve.stations for curve in a.field.curves])
        stations = a.sign*(chart-a.origin_station)
        return np.unique(np.r_[0., stations[(stations > 0.) & (stations < a.length)], a.length])

    def gradient(self, station, q, interval):
        a = self.approach
        interval = np.asarray(interval, dtype=float)
        station = np.asarray(station, dtype=float)
        if (interval.shape != (2,) or not np.isfinite(interval).all()
                or min(interval) < 0. or max(interval) > a.length or interval[0] == interval[1]
                or not np.isfinite(station).all() or not np.isfinite(q)):
            raise ValueError('Paired source jet needs one finite original station piece')
        margin = 128.*np.finfo(float).eps*max(1., a.length)
        if np.any(station < min(interval)-margin) or np.any(station > max(interval)+margin):
            raise ValueError('Paired source jet query leaves its declared station piece')
        xy = a.point(station, .5*a.width*q)
        return a.field.gradient(xy, chord_interval=tuple(map(float, a.chart_station(interval))))
