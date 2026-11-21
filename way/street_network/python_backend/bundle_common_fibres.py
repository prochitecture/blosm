"""Authoritative profile interfaces for a finite, source-owned Bundle body."""
import numpy as np


class CommonFibreCuts:
    """Validate two prepared transverse cuts without fitting or extending them.

    Keys are local assembly identities. A shared crest receives the same key
    on both pavement halves; separated peaks remain distinct. Source Section
    IDs and original junction/derived-cut ownership stay with the caller.
    """

    def __init__(self, frame, field, cuts, namespace):
        if set(cuts) != {(role, end) for role in (0, 1) for end in (0, 1)}:
            raise ValueError('A finite Bundle body needs both complete cuts of both members')
        self.cuts, self.keys = {}, {}
        stations, merged = [], []
        for end in (0, 1):
            values = []
            for role in (0, 1):
                points = np.array(cuts[role, end], dtype=float, copy=True)
                if (points.ndim != 2 or points.shape[1] != 3 or len(points) < 2
                        or not np.all(np.isfinite(points))):
                    raise ValueError('A finite body cut needs a complete finite profile')
                local = frame.local(points[:, :2])
                if np.ptp(local[:, 0]) > 1e-8:
                    raise ValueError('A finite body cut must be an actual common fibre')
                if np.any(np.diff(local[:, 1]) * (1. if role == 0 else -1.) >= -1e-10):
                    raise ValueError('A finite body cut must run from its original outer curb to its crest')
                station = float(local[:, 0].mean())
                for ray in field.curves:
                    if station < ray.stations[0]-1e-8 or station > ray.stations[-1]+1e-8:
                        raise ValueError('A finite body cut would extrapolate an original curb')
                if np.max(abs(points[:, 2]-field.values(points[:, :2]))) > 1e-8:
                    raise ValueError('A finite body cut changed its prepared physical profile')
                outer = field.curves[0 if role == 0 else 3](station)
                peak = field.peaks(station)[role]
                expected = frame.world(np.array((station, peak)))
                if max(np.linalg.norm(points[0, :2]-outer),
                       np.linalg.norm(points[-1, :2]-expected)) > 1e-8:
                    raise ValueError('A finite body cut lacks its actual curb or crown endpoint')
                self.cuts[role, end] = points
                self.keys[role, end] = [(namespace, 'body cut', end, role, i) for i in range(len(points))]
                values.append(station)
            if abs(values[0]-values[1]) > 1e-8:
                raise ValueError('Paired body profiles must share the same actual fibre')
            stations.append(.5*(values[0]+values[1]))
            first, last = (self.cuts[role, end][-1] for role in (0, 1))
            if np.linalg.norm(first-last) <= 1e-8:
                # Coordinate agreement alone does not select ownership. These
                # are already explicit halves of one supplied physical cut.
                if not np.array_equal(first, last):
                    raise ValueError('A merged body cut needs one identical authoritative crest sample')
                token = namespace, 'body cut', end, 'shared crest'
                for role in (0, 1):self.keys[role, end][-1] = token
                merged.append(end)
            elif frame.local(first[:2])[1] <= frame.local(last[:2])[1]:
                raise ValueError('Separate body profiles overlap without their common crest')
        if stations[1] <= stations[0]+1e-8:
            raise ValueError('A finite Bundle body needs two distinct ordered fibres')
        self.interval = tuple(stations)
        self.merged_ends = tuple(merged)
