"""Source-certified coordinate candidates from actual paired mouth normals.

The established chart remains first. Intersecting two declared mouth normals
chooses coordinates only; it creates no source incidence or physical geometry.
Every complete original source offset must certify positive angular progress.
"""
import numpy as np
from .bundle_axis import AxialChartUnavailable
from .bundle_angular_curves import AngularCoordinates
from .merged_run_inventory import inventory_frame
from .model import Key


def _station(layout, approach):
    indexes = [i for i, incoming in enumerate(layout.approaches) if incoming is approach]
    if len(indexes) != 1:
        raise ValueError('Mouth-frame coordinates need an actual original prepared approach')
    return float(layout.outward_stations[indexes[0]])


def candidate_frames(context, bundle, approach, specification):
    """Yield the unchanged established chart, then at most 32 certified poles."""
    failure = None
    try:
        original = inventory_frame(context, bundle, approach)
    except AxialChartUnavailable as error:
        failure = str(error)
    else:
        yield original, dict(family='established')
    if not hasattr(approach.frame_chart, 'axis'):
        return
    layout = specification.group.physical_layout
    at = _station(layout, approach)
    root = np.asarray(approach.point(at, 0.), float)
    root_axis = np.asarray(approach.frame_chart.axis, float)
    root_normal = np.array((-root_axis[1], root_axis[0]))
    visits = tuple(v for member in approach.run.members for v in member.visits)
    if len(visits) != len({v.section_id for v in visits}):
        raise ValueError('Mouth-frame coordinates repeat an original run Section')
    paths = tuple((context[Key('path', v.section_id)], v.forward) for v in visits)
    points = np.concatenate([path.points for path, _ in paths])
    origin = .5*(points.min(axis=0)+points.max(axis=0))
    reach = max(abs(path.section.offset)+.5*path.section.width for path, _ in paths)
    count = 0
    for side in specification.paired_sources:
        for entrance in side.entrances:
            if count >= 32: return
            if not hasattr(entrance.frame_chart, 'axis'): continue
            station = _station(side.group.physical_layout, entrance)
            point = np.asarray(entrance.point(station, 0.), float)
            axis = np.asarray(entrance.frame_chart.axis, float)
            normal = np.array((-axis[1], axis[0]))
            matrix = np.column_stack((normal, -root_normal))
            if abs(float(np.linalg.det(matrix))) <= 1e-8: continue
            parameters = np.linalg.solve(matrix, root-point)
            center = point+parameters[0]*normal
            radius = np.linalg.norm(points-center, axis=1)
            if not np.isfinite(center).all() or not np.isfinite(radius).all() or radius.min() <= 2.*reach:
                continue
            frame = object.__new__(AngularCoordinates)
            frame.center = center
            frame.scale = float(radius.mean())
            radial = origin-center
            frame.angle = float(np.arctan2(radial[1], radial[0]))
            radial_root = root-center
            frame.direction = float(np.sign(radial_root[0]*root_axis[1]-radial_root[1]*root_axis[0]))
            if not frame.direction: continue
            try:
                proof = min(frame.certify(path, forward=forward) for path, forward in paths)
            except AxialChartUnavailable:
                continue
            count += 1
            yield frame, dict(family='actual_paired_mouth_normals', root_owner=specification.group.members[0].text(),
                side_owner=side.owner.text(), side_end=int(entrance.end_index),
                center=center.tolist(), scale=frame.scale, whole_original_progress_lower=float(proof),
                original_sections=[v.section_id for v in visits], physical_geometry_changed=False,
                established_chart_failure=failure)
