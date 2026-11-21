"""Explicit shared original-station preparation before physical mouths qualify.

These functions return new preparation geometry. They never edit an accepted
mouth or Product. Source IDs and station constraints come from explicit
original Path inventories; coordinate proximity is not a source identity.
"""
from dataclasses import dataclass, replace
import hashlib
import numpy as np


def source_station_knots(path, values, *, protected_stations=()):
    """Retain original control knots and one sample per numerical duplicate.

    This is a preparation partition. Every retained station is an existing
    input or authored control station; no accepted source point is averaged.
    """
    values = np.asarray(tuple(values), dtype=float)
    if (values.ndim != 1 or not np.isfinite(values).all()
            or np.any(values < 0.) or np.any(values > path.length)):
        raise ValueError('Compulsory samples must remain on their original finite Path')
    protected = {0., float(path.length), *map(float, protected_stations)}
    if any(not 0. <= value <= path.length for value in protected):
        raise ValueError('Protected source cuts must remain on their original finite Path')
    if path.exact:
        protected.update(map(float, path.stations))
    ordered = sorted(set(map(float, values)) | protected)
    roundoff = min(1e-9, max(1e-10, 256.*np.spacing(max(1., path.length))))
    groups, current = [], []
    for station in ordered:
        if current and station-current[0] > roundoff:
            groups.append(current)
            current = []
        current.append(station)
    if current:
        groups.append(current)
    result = []
    for group in groups:
        mandatory = [s for s in group if s in protected]
        if len(mandatory) > 1:
            raise ValueError('Distinct mandatory source cuts have unresolved numerical incidence')
        result.append(mandatory[0] if mandatory else group[0])
    return np.asarray(result)


def sample_source_parameters(path, first, last, spacing, tolerance, compulsory, *, maximum_samples=20000,
                             protected_stations=()):
    """Idempotent supplied-knot sampling under the original offset-curve bound."""
    import math
    if (not math.isfinite(spacing) or spacing <= 0. or not math.isfinite(tolerance) or tolerance <= 0.
            or not 0. <= first < last <= path.length or type(maximum_samples) is not int or maximum_samples < 2):
        raise ValueError('Shared source sampling needs a bounded finite original interval')
    if not path.exact:
        reach = abs(path.section.offset)+.5*path.section.width
        if reach*path.curvature_bound >= 1.:
            raise ValueError('Shared source sampling requires regular original offset curvature')
        curvature = path.curvature_bound/(1.-reach*path.curvature_bound)
        speed = 1.+reach*path.curvature_bound
        if curvature > 1e-12:
            spacing = min(spacing, math.sqrt(8.*tolerance/curvature)/speed,
                          math.pi/(2.*curvature)/speed)
    knots = source_station_knots(path, (first, last, *compulsory), protected_stations=protected_stations)
    first = float(knots[np.argmin(abs(knots-first))])
    last = float(knots[np.argmin(abs(knots-last))])
    knots = knots[(knots >= first) & (knots <= last)]
    if len(knots) > maximum_samples:
        raise ValueError('Shared original source sampling exceeds its finite sample budget')
    values = []
    margin = max(64.*np.finfo(float).eps*max(1., path.length), spacing*1e-12)
    target = spacing-margin
    if target <= 0.:
        raise ValueError('Shared source sampling cannot resolve its positive metric spacing')
    for a, b in zip(knots, knots[1:]):
        count = 1 if b-a <= spacing else max(2, int(math.ceil((b-a)/target)))
        if len(values)+count+1 > maximum_samples:
            raise ValueError('Shared original source sampling exceeds its finite sample budget')
        values.extend(np.linspace(a, b, count+1)[:-1])
    return np.r_[values, last]


def shared_source_parameters(entrances, *, maximum_per_source=20000):
    if type(maximum_per_source) is not int or not 2 <= maximum_per_source <= 20000:
        raise ValueError('Shared source partitions need their bounded original sample budget')
    originals, stations, protected = {}, {}, {}
    for entrance in entrances:
        if len(entrance.paths) != len(entrance.source_curve_stations):
            raise ValueError('A prepared entrance lost its original source station inventory')
        for path, values in zip(entrance.paths, entrance.source_curve_stations):
            sid = path.section.id
            if sid in originals and originals[sid] is not path:
                raise ValueError('Shared source partitions require the same original Path object')
            originals[sid] = path
            protected.setdefault(sid, set()).update(getattr(entrance, 'protected_source_parameters', {}).get(sid, ()))
            values = np.asarray(values, dtype=float)
            if (values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all()
                    or values.min() < 0. or values.max() > path.length
                    or not (np.all(np.diff(values) > 0.) or np.all(np.diff(values) < 0.))):
                raise ValueError('Shared source partition leaves its finite original Path')
            stations.setdefault(sid, set()).update(map(float, values))
            if len(stations[sid]) > maximum_per_source:
                raise ValueError('Shared original station partition exceeds its finite sample budget')
    result = {}
    for sid, values in sorted(stations.items()):
        knots = source_station_knots(originals[sid], values, protected_stations=protected.get(sid, ()))
        if len(knots) > maximum_per_source:
            raise ValueError('Shared original station partition exceeds its finite sample budget')
        result[sid] = tuple(map(float, knots))
    return result


def verify_shared_source_phase(entrances):
    records = {}
    for entrance in entrances:
        for path, values in zip(entrance.paths, entrance.source_curve_stations):
            records.setdefault(path.section.id, []).append((path, np.sort(np.asarray(values, dtype=float))))
    checks = []
    for sid, evaluations in sorted(records.items()):
        for first, (path, values) in enumerate(evaluations):
            for other_path, other in evaluations[first+1:]:
                if other_path is not path:
                    raise ValueError('Shared phase certificate lost original Path identity')
                low, high = max(values[0], other[0]), min(values[-1], other[-1])
                if high <= low:
                    continue
                a = values[(values >= low) & (values <= high)]
                b = other[(other >= low) & (other <= high)]
                if not np.array_equal(a, b) or not len(a) or a[0] != low or a[-1] != high:
                    raise ValueError(f'Original Section {sid} needs canonical shared-prefix re-preparation')
                checks.append(dict(section=sid, first=float(low), last=float(high), samples=len(a),
                    station_digest=hashlib.sha256(a.tobytes()).hexdigest()))
    return dict(original_path_identity=True, overlapping_station_arrays_exact=True, shared_intervals=checks)


@dataclass(frozen=True)
class OneSidedPreparationRequest:
    site: object
    local_layout: object
    near_visits: tuple
    opposite_visit: object
    near_role: int


def canonical_reprepare(context, requests, *, additional_entrances=()):
    """Build fresh source-matched proposals; caller must hydrate fresh Products.

    Additional terminal entrances contribute compulsory source knots, but their
    caller must also rebuild them from the returned partitions before claiming
    a joint phase certificate with the returned one-sided preparations.
    """
    from .one_sided_paired_junctions import OneSidedPairedSpecification
    requests = tuple(requests)
    if not requests or len(requests) > 32:
        raise ValueError('Canonical one-sided preparation needs one to thirty-two explicit original sites')
    if len({request.site.key for request in requests}) != len(requests):
        raise ValueError('Canonical preparation repeats an original site')
    def prepare(request, parameters=None, protected=None):
        return OneSidedPairedSpecification.build(context, request.site, request.local_layout,
            request.near_visits, request.opposite_visit, request.near_role, source_parameters=parameters,
            protected_source_parameters=protected)
    initial = tuple(prepare(request) for request in requests)
    parameters = shared_source_parameters((*additional_entrances,
        *(entrance for spec in initial for entrance in spec.entrances)))
    protected = {}
    for entrance in (*additional_entrances, *(a for spec in initial for a in spec.entrances)):
        for sid, values in getattr(entrance, 'protected_source_parameters', {}).items():
            protected.setdefault(sid, set()).update(values)
    prepared = tuple(prepare(request, parameters, protected) for request in requests)
    proof = verify_shared_source_phase(tuple(entrance for spec in prepared for entrance in spec.entrances))
    parameters = shared_source_parameters((*additional_entrances,
        *(entrance for spec in prepared for entrance in spec.entrances)))
    prepared = tuple(replace(spec, proof=dict(spec.proof, shared_source_phase=proof)) for spec in prepared)
    return prepared, parameters, proof
