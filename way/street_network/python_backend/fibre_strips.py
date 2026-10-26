"""Direct, width-aware quad strips for capped common-fibre separator islands.

The source curb polylines are fixed. Added rows subdivide their existing
segments and are returned to the pavement builder as coordinated samples.
"""
from dataclasses import dataclass
import math
import numpy as np
from .mesh import Assembler, cross


@dataclass
class FibreStripResult:
    mesh: object
    extra_first: tuple
    extra_second: tuple
    metadata: dict


class CappedFibreFrame:
    def __init__(self, plan, island, number, spacing):
        self.plan, self.island, self.number, self.spacing = plan, island, number, spacing
        a, b = island.first.points - plan.origin, island.second.points - plan.origin
        self.original = a @ plan.axis
        self.upper, self.lower = a @ plan.left, b @ plan.left
        if (np.any(np.diff(self.original) <= 0.) or
                np.max(abs(b @ plan.axis - self.original)) > 1e-7):
            raise ValueError('Capped fibre strips need coordinated monotone source rows')
        self.across = max(4, 4 * int(math.ceil(np.linalg.norm(
            island.first.points - island.second.points, axis=1).max() / (4. * spacing))))
        n = len(self.original) - 1
        split = set()
        if island.collapsed[0]:
            split.update((1, 2))
        if island.collapsed[1]:
            split.update((n - 3, n - 2))
        rows = [self.original[0]]
        for i, (x, y) in enumerate(zip(self.original, self.original[1:])):
            if i in split:
                rows.extend((x + (y - x) / 3., x + 2. * (y - x) / 3.))
            rows.append(y)
        self.rows = np.asarray(rows)

    @staticmethod
    def smooth(t):
        t = np.clip(t, 0., 1.)
        return t ** 3 * (10. + t * (-15. + 6. * t))

    def refine(self):
        self.rows = np.r_[np.column_stack((self.rows[:-1],
            (2. * self.rows[:-1] + self.rows[1:]) / 3.,
            (self.rows[:-1] + 2. * self.rows[1:]) / 3.)).ravel(), self.rows[-1]]

    def point(self, x, v):
        x, v = np.broadcast_arrays(np.asarray(x), np.asarray(v))
        source_x = x.copy()
        width = np.interp(x, self.original, self.upper - self.lower)
        fade = np.maximum(0., 1. - (width / (.6 * self.spacing)) ** 2)
        # The cap's three profile vertices must not be collinear. A wider
        # physical opening may already exceed the ordinary shear fade width
        # at its first row, particularly at fine resolution. Keep a bounded
        # interior shear there; physical curb positions still have zero shear.
        cap_row = ((x == self.rows[1]) if self.island.collapsed[0] else np.zeros(x.shape, dtype=bool))
        if self.island.collapsed[1]:
            cap_row |= x == self.rows[-2]
        fade = np.where(cap_row, np.maximum(fade, .5), fade)
        delta = np.zeros_like(x)
        length = self.rows[-1] - self.rows[0]
        if self.island.collapsed[0]:
            weight = (1. - self.smooth((x - self.rows[0] - .25 * length) / (.25 * length))
                      if self.island.collapsed[1] else np.clip((self.rows[-1] - x) / self.spacing, 0., 1.))
            delta += np.minimum((x - self.rows[0]) / 7.8, .28 * self.spacing) * weight
        if self.island.collapsed[1]:
            weight = (self.smooth((x - self.rows[0] - .5 * length) / (.25 * length))
                      if self.island.collapsed[0] else np.clip((x - self.rows[0]) / self.spacing, 0., 1.))
            delta -= np.minimum((self.rows[-1] - x) / 7.8, .28 * self.spacing) * weight
        x = x + delta * fade * (1. - abs(2. * v - 1.))
        upper, lower = np.interp(x, self.original, self.upper), np.interp(x, self.original, self.lower)
        result = self.plan.world(x, (1. - v) * upper + v * lower)
        # Physical rims use the original world-space polylines verbatim.
        # This also avoids a round trip through the local chart frame.
        for side, points in ((0., self.island.first.points), (1., self.island.second.points)):
            selected = v == side
            if np.any(selected):
                sampled = np.column_stack([np.interp(source_x[selected], self.original, points[:, c]) for c in (0, 1)])
                result[selected] = sampled
        for row, collapsed in ((0, self.island.collapsed[0]), (-1, self.island.collapsed[1])):
            if not collapsed:
                selected = source_x == self.original[row]
                if np.any(selected):
                    result[selected] = ((1. - v[selected, None]) * self.island.first.points[row] +
                                        v[selected, None] * self.island.second.points[row])
        return result

    def address(self, row, v):
        return ('fibre island', self.number, 'row', row, round(float(v), 12))

    def columns(self, row):
        first = 1 if self.island.collapsed[0] else 0
        last = len(self.rows) - 2 if self.island.collapsed[1] else len(self.rows) - 1
        if (row == first and self.island.collapsed[0]) or (row == last and self.island.collapsed[1]):
            return 2
        if row in (0, len(self.rows) - 1):
            return self.across
        width = float(np.interp(self.rows[row], self.original, self.upper - self.lower))
        return min(self.across, max(2, 2 * int(math.ceil(width / (2. * self.spacing)))))

    def extras(self, role):
        boundary = self.island.first if role == 0 else self.island.second
        result = []
        for x in self.rows:
            i = int(np.clip(np.searchsorted(self.original, x, side='right') - 1, 0, len(self.original) - 2))
            t = float((x - self.original[i]) / (self.original[i + 1] - self.original[i]))
            if t < 1e-11 or t > 1. - 1e-11:
                continue
            point = self.point(np.array((x,)), np.array((float(role),)))[0]
            result.append((point, ('axial extra', self.number, i, round(t, 12), role)))
        return tuple(result)


class CappedFibreStrips:
    def __init__(self, *, max_cells=50000, max_refinements=2):
        if (type(max_cells) is not int or max_cells < 1 or type(max_refinements) is not int or max_refinements < 0):
            raise ValueError('Capped fibre strips need finite cell and refinement budgets')
        self.max_cells, self.max_refinements = max_cells, max_refinements

    def assemble(self, frames, context, height):
        assembly, tips = Assembler(), []
        for frame in frames:
            def add(xy, faces, addresses):
                if len(assembly.quads) + len(faces) > self.max_cells:
                    raise ValueError('Capped fibre strips exceed their cell budget')
                assembly.add(np.column_stack((xy, np.full(len(xy), height))), faces, addresses)
            for start, collapsed in ((True, frame.island.collapsed[0]), (False, frame.island.collapsed[1])):
                if not collapsed:
                    continue
                row = 1 if start else len(frame.rows) - 2
                root = 0 if start else len(frame.rows) - 1
                profile = frame.point(frame.rows[row], np.array((0., .5, 1.)))
                tip = frame.point(np.array((frame.rows[root],)), np.array((.5,)))[0]
                order = (2, 1, 0) if start else (0, 1, 2)
                addresses = [frame.address(root, .5), *(frame.address(row, (0., .5, 1.)[j]) for j in order)]
                add(np.vstack((tip, profile[list(order)])), np.array(((0, 1, 2, 3),)), addresses)
                tips.append(addresses[0])
            first = 1 if frame.island.collapsed[0] else 0
            last = len(frame.rows) - 2 if frame.island.collapsed[1] else len(frame.rows) - 1
            for i in range(first, last):
                pattern = context.patterns.strip(1, frame.columns(i), frame.columns(i + 1))
                u, v = pattern.uv.T
                x = (1. - u) * frame.rows[i] + u * frame.rows[i + 1]
                addresses = [frame.address(i if s == 0. else i + 1, t) if s in (0., 1.) else
                    ('fibre island', frame.number, 'cell', i, vi) for vi, (s, t) in enumerate(pattern.uv)]
                add(frame.point(x, v), pattern.quads[:, ::-1], addresses)
        return assembly.mesh('grass'), np.array([assembly.lookup[key] for key in tips], dtype=np.int32)

    def build(self, context, plan):
        frames = [CappedFibreFrame(plan, island, number, context.network.settings.spacing)
                  for number, island in enumerate(plan.islands)]
        failures = []
        for attempt in range(self.max_refinements + 1):
            try:
                mesh, tips = self.assemble(frames, context, plan.datum + plan.rise * plan.settings.exaggeration)
                quality = mesh.validate()
                e = np.roll(mesh.xyz[mesh.quads, :2], -1, axis=1) - mesh.xyz[mesh.quads, :2]
                sine = cross(e, np.roll(e, -1, axis=1)) / (np.linalg.norm(e, axis=-1) * np.linalg.norm(np.roll(e, -1, axis=1), axis=-1))
                tip_mask = np.isin(np.roll(mesh.quads, -1, axis=1), tips)
                intrinsic = float(sine[tip_mask].min()) if len(tips) else 1.
                ordinary = float(sine[~tip_mask].min())
                if quality['maximum_edge_ratio'] > 8.:
                    raise ValueError('Capped fibre strips need better edge grading')
                # A physical zero-gap cusp fixes its own angle. Other corners
                # must meet the ordinary limit or be no sharper than that cusp.
                limit = min(.03, intrinsic)
                if ordinary < limit * (1. - 1e-8):
                    raise ValueError('Capped fibre strips introduce an avoidable acute corner')
                return FibreStripResult(mesh,
                    tuple(item for frame in frames for item in frame.extras(0)),
                    tuple(item for frame in frames for item in frame.extras(1)),
                    dict(quality=quality, minimum_corner_sine=float(sine.min()),
                         minimum_ordinary_sine=ordinary, physical_cusp_minimum_sine=intrinsic,
                         ordinary_sine_limit=limit, refinements=attempt, previous_failures=failures,
                         strips=[dict(original_rows=len(f.original) - 1, rows=len(f.rows) - 1,
                                      maximum_columns=f.across) for f in frames],
                         construction='direct width-aware fibre strips with one-quad zero-gap caps'))
            except ValueError as exc:
                failures.append(str(exc))
                if 'budget' in str(exc) or attempt == self.max_refinements:
                    raise ValueError('Capped fibre separator remains unqualified: ' + str(exc)) from exc
                for frame in frames:
                    frame.refine()
