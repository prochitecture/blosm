"""Direct finite island tips and an independently sampled quad body."""
import math
import json
import numpy as np
from .mesh import Mesh, Assembler
from .bundle_axis import _regular_rows, _simplify
from .bundle_angular_curves import RayCurb
from .bundle_terminal_caps import sample_chain
from .four_side_charts import FourSidedQuadChart, quality
from .quad_charts import FixedRimQuality
from .validation import validate_overlaps

class AngularIslandBuilder:
    """Keep actual island boundaries; terminate zero-width tips in direct quads."""

    def __init__(self, frame, patterns, owner, *, spacing, curve_tolerance, height,
                 max_cells=50000, boundary_error=None):
        (self.frame, self.patterns, self.owner) = (frame, patterns, owner)
        (self.spacing, self.curve_tolerance) = (spacing, curve_tolerance)
        (self.height, self.max_cells) = (height, max_cells)
        self.boundary_error = .4*curve_tolerance if boundary_error is None else float(boundary_error)
        if not math.isfinite(self.boundary_error) or not 0.<=self.boundary_error<curve_tolerance:
            raise ValueError('Island source boundary exceeds its curve-accuracy budget')
        available = curve_tolerance-self.boundary_error-1e-8
        if available<=0.:
            raise ValueError('Island source boundary leaves no positive sampling-error budget')
        self.cap_error_budget = min(.6*curve_tolerance,available)
        self.body_tolerance = .8*curve_tolerance if boundary_error is None else min(.8*curve_tolerance,available)

    def _trimmed(self, ray, first, last, smooth=False, simplification=0.):
        middle = (ray.stations > first + 1e-08) & (ray.stations < last - 1e-08)
        points = np.vstack((ray(first), ray.points[middle], ray(last)))
        if not smooth:
            if simplification:
                points = _simplify(points,simplification)
            return _regular_rows(points, self.spacing, multiple=1)
        count = max(2, int(math.ceil(np.linalg.norm(np.diff(points, axis=0), axis=1).sum() / self.spacing)))
        for _ in range(8):
            if count > 4096:
                raise ValueError('Regular island sampling exceeds its allocation budget')
            try:
                return sample_chain(points, count, self.body_tolerance)
            except ValueError:
                count = int(math.ceil(1.25 * count))
        raise ValueError('The island body has no bounded-error regular sampling')

    def _cap_accuracy(self, ray, first, last):
        points = np.vstack((ray(first), ray.points[(ray.stations > first) & (ray.stations < last)], ray(last)))
        (a, b) = points[[0, -1]]
        delta = b - a
        t = np.clip((points - a) @ delta / float(delta @ delta), 0.0, 1.0)
        return float(np.linalg.norm(points - a - t[:, None] * delta, axis=1).max())

    def _candidate(self, upper, lower, distance, shear, smooth=False, *, cap_error_budget=None,
                   body_simplification=0.):
        rays = [RayCurb(self.frame, p) for p in (upper, lower)]
        (start, finish) = (max((r.stations[0] for r in rays)), min((r.stations[-1] for r in rays)))
        (first, last) = (start + distance, finish - distance)
        if last <= first + self.spacing:
            raise ValueError('Finite island caps consume the body')
        error = max((self._cap_accuracy(r, a, b) for r in rays for (a, b) in ((start, first), (last, finish))))
        cap_error_budget = self.cap_error_budget if cap_error_budget is None else cap_error_budget
        if error > cap_error_budget + 1e-08:
            raise ValueError(f'Tip chords exceed the curve-accuracy budget: {error}')
        (top, bottom) = (self._trimmed(r, first, last, smooth,body_simplification) for r in rays)
        if (len(top) - len(bottom)) % 2:
            at = int(np.argmax(np.linalg.norm(np.diff(top, axis=0), axis=1)))
            top = np.insert(top, at + 1, 0.5 * (top[at] + top[at + 1]), axis=0)
        (top_uv, bottom_uv) = (self.frame.local(p) for p in (top, bottom))
        middle_first = 0.5 * (top_uv[0] + bottom_uv[0]) + np.array((distance * shear, 0.0))
        middle_last = 0.5 * (top_uv[-1] + bottom_uv[-1]) - np.array((distance * shear, 0.0))
        mids = self.frame.world(np.asarray((middle_first, middle_last)))
        rim = np.vstack((bottom, mids[1:2], top[::-1], mids[0:1]))
        corners = (0, len(bottom) - 1, len(bottom) + 1, len(bottom) + len(top))
        uv = self.frame.local(rim)
        builder = FourSidedQuadChart(self.patterns, maximum_ratio=8.0, max_cells=self.max_cells)
        try:
            (mapped, counts) = builder.map(uv, corners)
            initial = quality(mapped)
            improvement = None
            if not builder.accepted(initial):
                (mapped, improvement) = FixedRimQuality(iterations=60, power=2).improve(mapped)
                if not builder.accepted(quality(mapped)):
                    raise ValueError(f'The finite island body has no accepted grading: {quality(mapped)}')
        except ValueError:
            # A short authored rim interval can defeat the first atlas while
            # another cached direct-quad topology fits the identical boundary.
            # Preserve every rim vertex and the existing shape/cell limits.
            alternative = builder.build(uv)
            mapped = alternative.mesh
            counts = alternative.metadata['rim_intervals']
            initial = alternative.metadata['initial']
            improvement = dict(fixed_rim_atlas=alternative.metadata)
        world = self.frame.world(mapped.xyz[:, :2])
        world[:len(rim)] = rim
        body = Mesh(np.c_[world, np.full(len(world), self.height)], mapped.quads, 'grass')
        top_keys = [('island', 'upper', i) for i in range(len(top))]
        bottom_keys = [('island', 'lower', i) for i in range(len(bottom))]
        midpoint_keys = [('island', 'cut_middle', i) for i in (0, 1)]
        body_keys = bottom_keys + [midpoint_keys[1]] + top_keys[::-1] + [midpoint_keys[0]]
        body_keys.extend((('island', 'body', i) for i in range(len(rim), len(world))))
        assembly = Assembler()
        assembly.add(body.xyz, body.quads, body_keys)
        root = np.asarray((upper[0], upper[-1]))
        for side in (0, 1):
            if side == 0:
                points = np.vstack((root[0], bottom[0], mids[0], top[0]))
                keys = [('island_tip', 0), bottom_keys[0], midpoint_keys[0], top_keys[0]]
            else:
                points = np.vstack((root[1], top[-1], mids[1], bottom[-1]))
                keys = [('island_tip', 1), top_keys[-1], midpoint_keys[1], bottom_keys[-1]]
            assembly.add(np.c_[points, np.full(4, self.height)], np.array(((0, 1, 2, 3),)), keys)
        mesh = assembly.mesh('grass')
        report = mesh.validate()
        if report['maximum_edge_ratio'] > 8.0:
            raise ValueError(f"Finite island world ratio exceeds 8: {report['maximum_edge_ratio']}")
        validate_overlaps({self.owner: mesh}, cell_size=2.0)
        rims = [(np.vstack((root[:1], top, root[1:])), [('island_tip', 0)] + top_keys + [('island_tip', 1)]), (np.vstack((root[:1], bottom, root[1:])), [('island_tip', 0)] + bottom_keys + [('island_tip', 1)])]
        return (mesh, rims, dict(construction='direct one-quad physical tips and unequal-rim body', quality=report, quads=len(mesh.quads), cap_distance=distance, shear=shear, boundary_chord_error=error, body_counts=counts, initial=initial, improvement=improvement))

    def build(self, upper, lower, *, smooth=False):
        failures = []
        for distance in (0.5 * self.spacing, 5.0 / 7.0 * self.spacing, self.spacing, 9.0 / 7.0 * self.spacing):
            for shear in (1.0 / 6.0, 1.0 / 5.0, 0.25):
                try:
                    (mesh, rims, report) = self._candidate(upper, lower, distance, shear, smooth)
                    report['failures'] = failures
                    return (mesh, rims, report)
                except ValueError as error:
                    failures.append(dict(distance=distance, shear=shear, error=str(error)))
        # Preserve the established selection first. A caller with a tighter
        # proved source-boundary error may spend its remaining budget on the
        # tip chord, without changing the overall physical-curve tolerance.
        # Exact source polygons often need none of the old 40% reserve.
        remaining = self.curve_tolerance-self.boundary_error-1e-8
        if remaining>.6*self.curve_tolerance:
            for distance in (.5*self.spacing,5./7.*self.spacing,self.spacing,9./7.*self.spacing):
                for shear in (1./6.,1./5.,.25):
                    try:
                        mesh,rims,report = self._candidate(upper,lower,distance,shear,smooth,
                                                          cap_error_budget=remaining)
                        report.update(failures=failures,source_boundary_error=self.boundary_error,
                                      cap_error_budget=remaining,
                                      complete_cap_error_bound=self.boundary_error+report['boundary_chord_error'])
                        return mesh,rims,report
                    except ValueError as error:
                        failures.append(dict(distance=distance,shear=shear,error=str(error),
                                             cap_error_budget=remaining))
        if not smooth:
            # Microscopic almost-collinear source fragments can force very
            # short boundary cells. Remove only a proved sub-tolerance amount;
            # substantial Exact curb bends remain mandatory. The replacement
            # pavement/grass interface is generated once and shared exactly.
            simplification = min(.1*self.curve_tolerance,remaining)
            for distance in (.5*self.spacing,5./7.*self.spacing,self.spacing,9./7.*self.spacing):
                for shear in (1./6.,1./5.,.25):
                    try:
                        mesh,rims,report = self._candidate(upper,lower,distance,shear,False,
                            cap_error_budget=remaining,body_simplification=simplification)
                        report.update(failures=failures,source_boundary_error=self.boundary_error,
                            cap_error_budget=remaining,body_simplification=simplification,
                            complete_curve_error_bound=self.boundary_error+max(
                                simplification,report['boundary_chord_error']))
                        return mesh,rims,report
                    except ValueError as error:
                        failures.append(dict(distance=distance,shear=shear,error=str(error),
                                             body_simplification=simplification))
        from .bundle_angular_island_tiles import AngularIslandTiles
        try:
            mesh,rims,report = AngularIslandTiles(self).build(upper,lower,smooth=smooth)
            report['established_cap_failures'] = failures
            return mesh,rims,report
        except ValueError as error:
            raise ValueError('No qualified finite island cap: '+json.dumps(failures)
                             +'; local tiles: '+str(error)) from error
