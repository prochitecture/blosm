"""Whole merged-profile mouth transported to its run body cut."""
import math
from types import SimpleNamespace
import numpy as np
from .bundle_angular_open_islands import gap_intervals
from .entrance_slopes import EntranceSlopeGuard
from .mesh import Mesh
from .slope_projection import corner_matrices


class MergedProfileBridge:
    max_cells = 2048

    def build(self, context, approach, station, mouth):
        settings = context.network.settings
        start = float(station)
        finish = start+2.*settings.spacing
        if finish >= approach.length-1e-6:
            raise ValueError('Merged connector consumes its original finite approach')
        field = approach.field
        domain = SimpleNamespace(frame=field.frame, rays=field.curves, field=field)
        first, last = sorted(map(float, approach.chart_station(np.array((start, finish)))))
        regions, roots = gap_intervals(domain, first, last)
        if roots or any(positive for a, b, positive in regions):
            raise ValueError('Merged connector needs its split-profile event policy')
        initial = approach.cross_section(start)
        delta = mouth[-1, :2]-mouth[0, :2]
        q = 2.*((mouth[:, :2]-mouth[0, :2])@delta/(delta@delta))-1.
        q[[0, -1]] = -1., 1.
        if not np.all(np.diff(q) > 0.):
            raise ValueError('Merged connector needs an ordered whole mouth')
        sides = q <= initial.peak
        fractions = (q-initial.peak)/np.where(sides, 1.+initial.peak, 1.-initial.peak)
        across = len(q)-1
        anchors = [start, finish]
        if any(path.exact for path in approach.paths):
            # Retain actual Exact curb corners, not all collinear sample knots.
            for ray in (field.curves[0], field.curves[3]):
                edges = np.diff(ray.points, axis=0)
                direction = edges/np.linalg.norm(edges, axis=1)[:, None]
                turns = np.linalg.norm(np.diff(direction, axis=0), axis=1) > 1e-9
                values = (ray.stations[1:-1][turns]-approach.origin_station)/approach.sign
                anchors.extend(values[(values > start+1e-8) & (values < finish-1e-8)])
        anchors = np.unique(anchors)
        rows = np.r_[np.concatenate([np.linspace(a, b, max(1, int(math.ceil((b-a)/settings.spacing)))+1)[:-1]
                                    for a, b in zip(anchors, anchors[1:])]), finish]
        if (len(rows)-1)*across > self.max_cells:
            raise ValueError('Merged connector exceeds its bounded quad budget')
        pattern = context.patterns.grid(len(rows)-1, across)
        s = np.repeat(rows, across+1)
        peak = np.array([approach.cross_section(at).peak for at in rows])
        varying_q = peak[:, None]+fractions[None, :]*np.where(sides[None, :], 1.+peak[:, None], 1.-peak[:, None])
        xy = approach.point(s, .5*approach.width*varying_q.ravel())
        xyz = np.c_[xy, field.values(xy)]
        front = np.arange(across+1)
        back = np.arange(len(xyz)-across-1, len(xyz))
        error = float(abs(xyz[front]-mouth).max())
        if error > 1e-8:
            raise ValueError('Merged connector disagrees with its actual junction mouth')
        xyz[front] = mouth
        faces = pattern.quads[:, ::-1]
        mesh = Mesh(xyz, faces)
        quality = mesh.validate()
        published = Mesh(xyz.astype(np.float32), faces)
        published_quality = published.validate()
        if published_quality['maximum_edge_ratio'] > 8.:
            raise ValueError('Merged connector exceeds the publication edge ratio')
        gradient = np.einsum('fcij,fj->fci', corner_matrices(published), published.xyz[published.quads, 2]-field.datum)
        slope = float(np.linalg.norm(gradient, axis=-1).max())
        longitudinal = field.frame.longitudinal_bound(published)
        if slope > .1+1e-6 or longitudinal > .03+1e-8:
            raise ValueError('Merged connector exceeds physical slope limits')
        joins = []
        for ids, at in ((front, start), (back, finish)):
            velocity = np.broadcast_to(approach.frame(at)[1], (len(xyz), 2)).copy()
            local_q = varying_q[0 if at == start else -1]
            velocity[ids] = approach.jet(at, .5*approach.width*local_q)[1]
            joins.append(EntranceSlopeGuard().validate(published, ids, velocity, published.xyz[:, 2]-field.datum))
        if abs(mesh.xyz[:, 2]-field.values(mesh.xyz[:, :2])).max() > 1e-8:
            raise ValueError('Merged connector changed the actual dynamic crown')
        self.front_ids, self.back_ids = front, back
        self.body_cut = mesh.xyz[back].copy()
        return mesh, dict(quads=len(faces), quality=quality, published_quality=published_quality,
            exact_mouth=np.array_equal(mesh.xyz[front], mouth), physical_slope=slope,
            angular_longitudinal=longitudinal, entrance_bounds=joins, start=start, finish=finish,
            original_sections=[path.section.id for path in approach.paths], all_direct_quads=True,
            addon_task_integration=True)
