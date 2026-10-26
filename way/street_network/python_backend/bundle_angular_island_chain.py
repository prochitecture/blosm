"""Ordered separator islands joined by shared crown ridges, with direct quads."""
import copy
import numpy as np
from .bundle_axis import _simplify
from .bundle_angular_road_tiles import AngularRoadTiles
from .bundle_angular_surfaces import AngularSurfaceBuilder
from .mesh import Assembler


class AngularIslandChain:
    """Compose independently bounded islands without reusing their seam IDs."""
    def __init__(self, surface_builder, island_builder):
        self.surface, self.island = surface_builder, island_builder

    def build(self, loops, field, first, last, *, smooth=False):
        frame = self.island.frame
        pieces = []
        for loop in loops:
            upper, lower = self.surface.island_halves(loop, frame)
            if not smooth:
                upper, lower = (_simplify(points, 1e-9) for points in (upper, lower))
            if (np.linalg.norm(upper[0]-lower[0]) > 1e-8
                    or np.linalg.norm(upper[-1]-lower[-1]) > 1e-8):
                raise ValueError('Angular island halves disagree at their physical tips')
            stations = frame.local(upper[[0, -1]])[:, 0]
            if stations[1] <= stations[0]+1e-8:
                raise ValueError('An angular separator island has no positive station interval')
            pieces.append((stations, upper, lower))
        pieces.sort(key=lambda piece: piece[0][0])
        ends = frame.local(np.asarray((first, last)))[:, 0]
        previous = ends[0]
        for stations, _, _ in pieces:
            if stations[0] <= previous+1e-8 or stations[1] >= ends[1]-1e-8:
                raise ValueError('Angular islands need disjoint ordered intervals inside the free ends')
            previous = stations[1]
        if len(pieces) < 2:
            raise ValueError('An angular island chain needs at least two islands')

        assembly, rims, reports = Assembler(), [], []
        for index, (stations, upper, lower) in enumerate(pieces):
            builder = copy.copy(self.island)
            builder.max_cells -= len(assembly.quads)
            if builder.max_cells < 3:
                raise ValueError('Angular islands exceed their complete surface cell budget')
            mesh, boundary, report = builder.build(upper, lower, smooth=smooth)
            # Each builder deliberately uses local parametric keys. Prefix
            # the exported keys once; islands are distinct physical regions.
            boundary = [(xy, [(self.island.owner, 'island', index, key) for key in keys])
                        for xy, keys in boundary]
            assembly.add(mesh.xyz, mesh.quads)
            rims.append(boundary)
            reports.append(dict(index=index, station_interval=list(map(float, stations)), chart=report))

        high, low = [], []
        previous_point, previous_key = first, None
        for index, boundary in enumerate(rims):
            target, target_keys = boundary[0]
            ridge = self.surface.ridge(frame, field, previous_point, target[0],
                                       (self.island.owner, 'ridge', index))
            if previous_key is not None:
                ridge[1][0] = previous_key
            ridge[1][-1] = target_keys[0]
            high.extend((ridge, boundary[0]))
            low.extend((ridge, boundary[1]))
            previous_point, previous_key = target[-1], target_keys[-1]
        tail = self.surface.ridge(frame, field, previous_point, last,
                                  (self.island.owner, 'ridge', len(rims)))
        tail[1][0] = previous_key
        high.append(tail)
        low.append(tail)
        mesh = assembly.mesh('grass')
        return (mesh, self.surface.join_chains(high), self.surface.join_chains(low),
                dict(construction='ordered finite-tip islands and common crown bridges',
                     islands=len(rims), crown_bridges=len(rims)+1, ordered_disjoint_intervals=True,
                     whole_interfaces_shared=True, charts=reports, quality=mesh.validate()))


class MultiIslandSurfaceBuilder(AngularSurfaceBuilder):
    """Keep island interfaces and use bounded local charts for the long body."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extended_grading = True

    def _direct_single(self, frame, field, outer, seam, name, *, grass=False):
        if not grass and name.endswith('/body'):
            try:
                mesh, keys, report = AngularRoadTiles(self, frame, field).build(outer, seam, name)
                report['selection'] = 'local charts for multiple-island body'
                return mesh, keys, report
            except ValueError as error:
                local_failure = str(error)
            mesh, keys, report = self._mapped_single(frame, field, outer, seam, name)
            report['local_body_failure'] = local_failure
            return mesh, keys, report
        return self._mapped_single(frame, field, outer, seam, name, grass=grass)

    def _mapped_single(self, frame, field, outer, seam, name, *, grass=False):
        try:
            return super()._mapped_single(frame, field, outer, seam, name, grass=grass)
        except ValueError as original:
            local = '/body/local ' in name
            if (grass or not (name.endswith(('/head', '/tail')) or local)
                    or (local and len(outer[0])+len(seam[0]) > 512)):
                raise
            failure = str(original)
        from .bundle_angular_caps import LayeredAngularPatch
        mesh, keys, report = LayeredAngularPatch(self, frame, field).build(outer, seam, name)
        report['four_side_failure'] = failure
        return mesh, keys, report
