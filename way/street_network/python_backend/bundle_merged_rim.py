"""Whole profile samples and certified slices of a merged physical curb."""
from dataclasses import replace
import math
import numpy as np
from .bundle_merged_entrances import PairedAngularEntrance
from .compounds import SampledFootprint, RimSegment
from .regions import Region
from .mesh import cross
from .paths import right, unit


def source_pieces(approach):
    if hasattr(approach.path, 'trims'):
        return tuple((piece.path, piece.forward, approach.path.stations[index])
                     for index, piece in enumerate(approach.path.pieces))
    return ((approach.path, approach.end == 0, 0.),)


def protected_parameters(boundary):
    mandatory = [0., 1., *list(boundary.lengths[1:-1]/boundary.length)]
    for index, (approach, parameters, cumulative, distance) in enumerate(getattr(boundary, 'leads', ())):
        candidates = [0.]
        if isinstance(approach, PairedAngularEntrance):
            high = (distance < 0.) == (approach.sign > 0.)
            rays = (0, 2) if high else (1, 3)
            if not all(path.exact for path in approach.paths):
                rays = ()
            for ray_index in rays:
                ray = approach.field.curves[ray_index]
                edges = np.diff(ray.points, axis=0)
                bend = abs(cross(edges[:-1], edges[1:])) > 1e-8*np.linalg.norm(edges[:-1], axis=1)*np.linalg.norm(edges[1:], axis=1)
                for coordinate in ray.stations[1:-1][bend]:
                    station = (coordinate-approach.origin_station)/approach.sign
                    if min(parameters) < station < max(parameters) and np.linalg.norm(
                            approach.point(station, distance)-ray(coordinate)) < 1e-8:
                        candidates.append(float(station))
        else:
            for path, forward, offset in source_pieces(approach):
                if path.exact:
                    candidates.extend(offset+(path.knots if forward else path.length-path.knots))
        for station in candidates:
            if not min(parameters)+1e-10 < station < max(parameters)-1e-10:
                continue
            length = (np.interp(station, parameters[::-1], cumulative[::-1]) if parameters[0] > parameters[-1]
                      else np.interp(station, parameters, cumulative))
            mandatory.append(float((length+(0. if index == 0 else boundary.lengths[2]))/boundary.length))
    knots = np.asarray(sorted(set(mandatory)))
    return knots[np.r_[True, np.diff(knots)>1e-10]]


class RetainedCurbPiece:
    def __init__(self, original, first, last):
        self.original, self.start, self.stop = original, float(first), float(last)
        self.first, self.last = original(np.asarray(first)), original(np.asarray(last))
        self.length = original.length*(last-first)
        self.corner, self.continuation = original.corner, self

    def __call__(self, t):
        return self.original(self.start+np.asarray(t)*(self.stop-self.start))

    def intervals(self, spacing, tolerance):
        original = self.original
        if original.corner.center is None:
            return max(1, int(math.ceil(self.length/spacing)))
        middle = .5*(self.start+self.stop)*original.length
        if original.lengths[1]-1e-10 <= middle <= original.lengths[2]+1e-10:
            angle = self.length/original.corner.radius
            maximum_angle = 2.*math.acos(1.-min(tolerance, original.corner.radius)/original.corner.radius)
            return max(1, int(math.ceil(self.length/spacing)), int(math.ceil(angle/maximum_angle)))
        index = 0 if middle < original.lengths[1] else 1
        approach, parameters, cumulative, distance = original.leads[index]
        delta = np.diff(cumulative)
        ratio = float(np.max(np.divide(abs(np.diff(parameters)), delta,
            out=np.zeros_like(delta), where=delta > 1e-12)))
        curvature, source_speed = 0., 1.
        if isinstance(approach, PairedAngularEntrance):
            for ray in approach.field.curves:
                edge = np.diff(ray.points, axis=0)
                normal = right(unit(edge))
                a, b = ray.points[:-1]-approach.frame_chart.center, ray.points[1:]-approach.frame_chart.center
                height = abs(np.sum(normal*a, axis=1))
                cosine = np.minimum(abs(np.sum(normal*unit(a), axis=1)), abs(np.sum(normal*unit(b), axis=1)))
                source_speed = max(source_speed, float(np.max(height/(cosine*cosine*approach.frame_chart.scale))))
            probes = self(np.array((.25, .5, .75)))
            chord = self.last-self.first
            fraction = ((probes-self.first)@chord)/(chord@chord)
            deviation = np.linalg.norm(probes-(self.first+fraction[:, None]*chord), axis=1)
            if deviation.max() > 1e-8:
                if any(path.exact for path in approach.paths):
                    raise ValueError('Merged Exact curb slice lost an actual physical breakpoint')
                curvature = max(path.curvature_bound/(1.-(abs(path.section.offset)+.5*path.section.width)*path.curvature_bound)
                                for path in approach.paths)
                tolerance *= .8
        else:
            for path, forward, offset in source_pieces(approach):
                if path.exact:
                    derivative = np.diff(path.normals, axis=0)/path.chords[:, None]
                    source_speed = max(source_speed, 1.+abs(distance)*float(np.linalg.norm(derivative, axis=1).max()))
                else:
                    k = path.curvature_bound
                    curvature = max(curvature, k/(1.-abs(distance)*k))
                    source_speed = max(source_speed, 1.+abs(distance)*k)
        speed = max(1., ratio*source_speed)
        step = spacing/speed
        if curvature > 1e-12:
            step = min(step, math.sqrt(8.*tolerance/curvature)/speed)
        return max(1, int(math.ceil(self.length/step)))


def merged_footprint(owner, layout, settings):
    rim, labels, curves = [], [], []
    for side, approach in enumerate(layout.approaches):
        left = approach.profile.coordinate(-1, np.linspace(0., 1., approach.profile.intervals(-1, layout.radial)+1))
        right = approach.profile.coordinate(1, np.linspace(0., 1., approach.profile.intervals(1, layout.radial)+1))
        u = .5*(np.r_[left, right[::-1][1:]]+1.)
        u = u[::-1]
        rim.extend(approach.point(layout.outward_stations[side], approach.offset+approach.width*(u-.5))[:-1])
        labels.extend(RimSegment('mouth', side, float(a), float(b)) for a, b in zip(u, u[1:]))
        boundary = layout.boundaries[side]
        knots = protected_parameters(boundary)
        for a, b in zip(knots, knots[1:]):
            piece = RetainedCurbPiece(boundary, a, b)
            count = piece.intervals(settings.spacing, settings.curve_tolerance)
            count += count%2
            if count > 20000 or len(rim)+count > 40000:
                raise ValueError('Merged physical rim exceeds its finite sample budget')
            t = np.linspace(0., 1., count+1)
            rim.extend(piece(t)[:-1])
            labels.extend(RimSegment('curb', len(curves), float(first), float(last)) for first, last in zip(t, t[1:]))
            curves.append(piece)
    return SampledFootprint(Region(owner, np.asarray(rim)), tuple(labels), replace(layout, boundaries=tuple(curves)))
