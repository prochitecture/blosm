"""Continuous complete-quad containment in explicitly supplied source strips."""
import math
import numpy as np
from .finite_containment import segments_on_polygon_side
from .model import Key
from .prepared_source_classifier import SourceProofClassifier
from .regions import Region, RegionUnion, _classify, signed_area
from .source_strip_containment import _distance


def _certify_side(mesh, boundary, inside, band, classifier):
    other = np.roll(boundary, -1, axis=0)
    edge_low, edge_high = np.minimum(boundary, other), np.maximum(boundary, other)
    error = ('A pavement quad leaves its finite original source strips' if inside else
             'A pavement quad fills a real original source gap')
    batch = max(1, min(128, 32768 // len(boundary)))
    for start in range(0, len(mesh.quads), batch):
        quads = mesh.xyz[mesh.quads[start:start+batch], :2]
        low = np.nextafter(quads.min(axis=1)-band, -np.inf)
        high = np.nextafter(quads.max(axis=1)+band, np.inf)
        nearby = ((low[:,None,0] <= edge_high[:,0]) & (high[:,None,0] >= edge_low[:,0]) &
                  (low[:,None,1] <= edge_high[:,1]) & (high[:,None,1] >= edge_low[:,1])).any(axis=1)
        # With no boundary segment in the expanded box, the entire box lies
        # on one side. A quad surrounding a whole hole always takes the near path.
        if np.any(~nearby):
            side, _ = classifier(quads[~nearby,0], boundary, band)
            if np.any(side != inside):
                raise ValueError(error)
        indices = np.flatnonzero(nearby)
        if not len(indices):
            continue
        selected = quads[indices]
        first = selected.reshape(-1, 2)
        last = np.roll(selected, -1, axis=1).reshape(-1, 2)
        if not np.all(segments_on_polygon_side(first, last, boundary, inside=inside,
                                               tolerance=band, classifier=classifier)):
            raise ValueError(error)
        for i in indices:
            points = boundary[np.all(boundary >= low[i], axis=1) & np.all(boundary <= high[i], axis=1)]
            if not len(points):
                continue
            polygon = quads[i]
            sides = np.stack((polygon, np.roll(polygon, -1, axis=0)), axis=1)
            hit, near = _classify(points, polygon, 1e-9)
            if np.any((hit | (near >= 0)) & (_distance(points, sides) > band)):
                raise ValueError('A pavement quad encloses an original source boundary')


def certify_mesh(mesh, curves, *, source_error, mesh_error, maximum_error,
                 max_cells=16000, max_segments=20000):
    """Check finite source polygons, including a hole fully inside a quad.

    Curves must be the actual two sides of each explicitly owned visible source
    span, in matching traversal order. The caller supplies its independently
    certified chord and mesh approximation budgets. No buffer, hull, inferred
    junction, filled hole or permissive geometry repair is constructed.
    """
    if (not all(math.isfinite(v) for v in (source_error, mesh_error, maximum_error))
            or min(source_error, mesh_error) < 0. or maximum_error <= 0.
            or source_error+mesh_error > maximum_error+1e-12
            or len(curves) % 2 or len(curves) == 0 or len(mesh.quads) > max_cells):
        raise ValueError('Invalid finite source-mesh containment budget')
    mesh.validate()
    operands = []
    for index in range(0, len(curves), 2):
        first, last = (np.asarray(curve.points if hasattr(curve, 'points') else curve)
                       for curve in curves[index:index+2])
        xy = np.vstack((first, last[::-1]))
        operands.append(Region(Key('body_source_strip', index//2), xy if signed_area(xy) > 0. else xy[::-1]))
    loops = RegionUnion(max_segments=max_segments).build(operands)
    outer = [loop.xy for loop in loops if not loop.hole]
    holes = [loop.xy for loop in loops if loop.hole]
    if len(outer) != 1:
        raise ValueError('A source-owned pavement piece needs one connected finite strip union')
    band = max(1e-7, mesh_error)
    if source_error+band > maximum_error+1e-12:
        raise ValueError('Finite source-mesh numerical band exceeds its complete error allocation')
    classifier = SourceProofClassifier()
    _certify_side(mesh, outer[0], True, band, classifier)
    for hole in holes:
        _certify_side(mesh, hole, False, band, classifier)
    return dict(proof='continuous_complete_quad_edges_inside_finite_original_strip_union',
        source_chord_error=source_error, mesh_boundary_band=band,
        maximum_combined_error=source_error+band, maximum_error=maximum_error,
        source_strips=len(operands), holes=len(holes), complete_quads=len(mesh.quads))
