"""Pure-Python camera pose interpretation. See README.md for source conventions."""
from dataclasses import dataclass
import math


def number(value):
    return (float(value) if isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(value) else None)


def matrix_product(a, b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)) for i in range(3))


def transpose(a):
    return tuple(zip(*a))


def rotation_vector(vector):
    if not isinstance(vector, (list, tuple)) or len(vector) != 3 or any(number(v) is None for v in vector):
        raise ValueError('Invalid angle-axis rotation vector')
    angle = math.hypot(*vector)
    if angle < 1e-12:
        return ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.))
    x, y, z = (v/angle for v in vector)
    c, s, d = math.cos(angle), math.sin(angle), 1-math.cos(angle)
    return ((c+x*x*d, x*y*d-z*s, x*z*d+y*s),
            (y*x*d+z*s, c+y*y*d, y*z*d-x*s),
            (z*x*d-y*s, z*y*d+x*s, c+z*z*d))


def coordinates(value):
    if isinstance(value, dict):
        value = value.get('coordinates')
    if not isinstance(value, (list, tuple)) or len(value) < 2:
        return None
    lon, lat = number(value[0]), number(value[1])
    return (lon, lat) if lon is not None and lat is not None and -180 <= lon <= 180 and -90 <= lat <= 90 else None


def projected_basis(projection, lon, lat):
    # Direction of local geodetic east in Blosm's map projection. Correct grid
    # convergence without scaling the camera axes with the projection's Jacobian.
    west = projection.fromGeographic(lat, lon-1e-5)
    east = projection.fromGeographic(lat, lon+1e-5)
    x, y = east[0]-west[0], east[1]-west[1]
    length = math.hypot(x, y)
    if not math.isfinite(length) or length < 1e-10:
        raise ValueError('Cannot determine camera orientation in this geographic projection')
    x, y = x/length, y/length
    return ((x, -y, 0.), (y, x, 0.), (0., 0., 1.))


@dataclass
class Pose:
    position: tuple
    rotation: tuple
    properties: dict
    lens: float
    sensor_height: float
    camera_type: str


def interpret(metadata, coverage, projection, height_offset=0.):
    props = {}
    point = None
    for key, source in (('computed_geometry', 'computed'), ('geometry', 'original')):
        point = coordinates(metadata.get(key))
        if point is not None:
            props['position_source'] = source
            break
    if point is None:
        point = coordinates(coverage.get('coverage_geometry'))
        props['position_source'] = 'coverage'
    if point is None:
        raise ValueError('Image has no geographic position')
    lon, lat = point
    x, y, _ = projection.fromGeographic(lat, lon)
    altitude = None
    for key in ('computed_altitude', 'altitude'):
        altitude = number(metadata.get(key))
        if altitude is not None:
            props['altitude_source'] = key
            break
    z = altitude-height_offset if altitude is not None else 0.
    if altitude is None:
        props['altitude_source'] = 'missing_scene_plane'
    else:
        props['source_altitude'] = altitude
    props.update(source_longitude=lon, source_latitude=lat, height_offset=height_offset,
                 vertical_datum='provider_unspecified; no geoid conversion')

    orientation = metadata.get('exif_orientation')
    if isinstance(orientation, bool) or orientation not in (1, 2, 3, 4, 5, 6, 7, 8):
        orientation = 1
        props['exif_assumed'] = True
    mirrored = orientation in (2, 4, 5, 7)
    # Mapillary Transform._createBasicWorldToCamera applies Rz(exif)*R.
    # Invert it once, then convert Blender (+X right,+Y up,-Z forward) axes
    # to SfM (+X right,+Y down,+Z forward). Mirroring is not a proper rotation.
    raw = metadata.get('computed_rotation')
    full_rotation = False
    try:
        r = rotation_vector(raw)
        angle = {1: 0., 3: math.pi, 6: math.pi/2, 8: 3*math.pi/2}.get(orientation, 0.)
        c, s = math.cos(angle), math.sin(angle)
        upright = ((c, -s, 0.), (s, c, 0.), (0., 0., 1.))
        r = matrix_product(transpose(matrix_product(upright, r)), ((1.,0.,0.), (0.,-1.,0.), (0.,0.,-1.)))
        props['orientation_source'] = 'computed_rotation'
        full_rotation = True
    except (ValueError, TypeError):
        heading = None
        for key, value in (('computed_compass_angle', metadata.get('computed_compass_angle')),
                           ('compass_angle', metadata.get('compass_angle')),
                           ('coverage_compass_angle', coverage.get('compass_angle'))):
            heading = number(value)
            if heading is not None:
                props['orientation_source'] = key
                break
        if heading is None:
            heading = 0.
            props['orientation_source'] = 'missing_assumed_north'
        theta = math.radians(heading)
        c, s = math.cos(theta), math.sin(theta)
        # Upright, horizontal assumption. Heading supplies no raw-image roll:
        # do not apply a second EXIF correction to this already-upright fallback.
        r = ((c, 0., -s), (-s, 0., -c), (0., 1., 0.))
        props['source_heading'] = heading
    r = matrix_product(projected_basis(projection, lon, lat), r)
    props.update(exif_orientation=orientation, exif_mirrored=mirrored,
                 pose_complete=bool(full_rotation and altitude is not None and not mirrored and not props.get('exif_assumed')),
                 pose_approximate=bool(not full_rotation or altitude is None or mirrored or props.get('exif_assumed')))

    kind = metadata.get('camera_type') or ('spherical' if coverage.get('is_pano') else 'unknown')
    params = metadata.get('camera_parameters')
    params = params if isinstance(params, (list, tuple)) and len(params) <= 32 and all(number(v) is not None for v in params) else []
    width, height = number(metadata.get('width')), number(metadata.get('height'))
    width = width if width is not None and width > 0 else 0.
    height = height if height is not None and height > 0 else 0.
    uw, uh = (height, width) if orientation >= 5 else (width, height)
    lens, sensor_height = 50., 24.
    projection_status = 'proxy_missing_intrinsics'
    if params and params[0] > 0 and uw and uh:
        lens = 36. * params[0]*max(width, height)/uw
        sensor_height = 36.*uh/uw
        if not 1 <= lens <= 5000 or not 1 <= sensor_height <= 100:
            lens, sensor_height = 50., 24.
        else:
            projection_status = ('pinhole_requires_image_aspect' if kind == 'perspective'
                                 else 'proxy_unsupported_camera_model')
            if kind == 'perspective' and any(params[1:]):
                projection_status = 'pinhole_without_radial_distortion'
    camera_type = 'PANO' if kind in ('spherical', 'equirectangular') else 'PERSP'
    if camera_type == 'PANO':
        projection_status = 'equirectangular_requires_supported_renderer'
    if mirrored:
        projection_status = 'proxy_mirrored_exif'
    props.update(camera_model=str(kind), camera_parameters=list(params), width=width, height=height,
                 upright_width=uw, upright_height=uh, projection_status=projection_status)
    return Pose((x, y, z), r, props, lens, sensor_height, camera_type)
