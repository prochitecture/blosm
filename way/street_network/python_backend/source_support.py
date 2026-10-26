"""Bounded original-source support for junctions close to an import rectangle.

This import preparation owns no mesh. It preserves source identities, keeps old
cut controls, and records every extra original interval. Only an unambiguous
same-way continuation may supply more geometry; proximity never connects roads.
"""
from dataclasses import dataclass, replace
import math
import numpy as np
from .model import End, Key
from .patterns import PatternLibrary
from .providers import PathProvider
from .profiles import ProfileProvider
from .routes import Route, RouteTopology, route_key
from .updates import UpdateCoordinator, BuildContext


@dataclass
class SourcePreparation:
    network: object
    report: dict
    error: str = ''


class SourceSupportError(ValueError):
    pass


class ClippedApproachSupport:
    """Recover finite same-way support only after ordinary contact fitting fails."""
    policy = 'BOUNDED_ORIGINAL_APPROACH_SUPPORT_V2'

    def __init__(self, snapshot, context, *, maximum_extension=40., clearance=1.4):
        if not all(math.isfinite(v) and v > 0 for v in (maximum_extension, clearance)):
            raise ValueError('Source support needs finite positive extension and clearance bounds')
        self.context = context or {}
        self.maximum_extension, self.clearance = maximum_extension, clearance
        self.sources = {}
        for street in snapshot.streets.values():
            for item in street['items']:
                if item['kind'] == 'Section':
                    if item['id'] in self.sources:
                        raise SourceSupportError('Original support has multiple owners for one Section')
                    self.sources[item['id']] = item

    @staticmethod
    def _error(network, report, error):
        report['status'] = 'UNRESOLVED'
        report['error'] = str(error)
        return SourcePreparation(network, report, str(error))

    def _cuts(self, network):
        if self.context and (self.context.get('version') != 1 or self.context.get('issues')):
            raise SourceSupportError('Original outside-source context is unresolved')
        attached = {e for j in (*network.intersections, *network.connections) for e in j.ends}
        sections = network.section_index
        result = {}
        for node in self.context.get('nodes', ()):
            for opened in node.get('open_ends', ()):
                end = End(opened['section_id'], opened['end'])
                if end.section_id not in sections or end in attached:
                    continue  # A recovered original junction already owns it.
                if end in result:
                    raise SourceSupportError('A selection cut has multiple original outside owners')
                source = self.sources[end.section_id]
                section = sections[end.section_id]
                records = source.get('source_segments', ())
                if not records:
                    raise SourceSupportError('Selection cut lacks its original source interval')
                record = records[0 if end.end == 0 else -1]
                expected = (str(opened['source_way_id']), *map(str, opened['source_nodes']))
                if (tuple(map(str, record[:3])) != expected or
                        abs(record[3 if end.end == 0 else 4] - opened['parameter']) > 1e-10):
                    raise SourceSupportError('Selection cut differs from its original directed source interval')
                cuts = [r for r in node['incidences'] if
                        (r['source_way_id'], *r['source_nodes']) == expected and
                        abs(r['parameter'] - opened['parameter']) <= 1e-10]
                if (len(cuts) != 1 or math.dist(section.points[0 if end.end == 0 else -1][:2],
                                              cuts[0]['clipped_point']) > .001):
                    raise SourceSupportError('Selection cut no longer matches its original clipped position')
                result[end] = node, opened, cuts[0]
        return result

    def _trim_context_legs(self, network, report):
        """Honor the explicit straight support trims prepared by boundary recovery.

        These full adjacent controls are context, not permission to publish the
        entire outside segment. Keep them in the preparation report. Curved
        metric cuts require their own source interpolation transport.
        """
        sections, furniture, details = list(network.sections), list(network.furniture), []
        for i, section in enumerate(sections):
            source = self.sources[section.id]
            trims = tuple(float(source.get(name, 0.)) for name in ('clip_start', 'clip_end'))
            if not any(trims):
                continue
            if not all(math.isfinite(v) and v >= 0 for v in trims) or len(section.points) != 2:
                raise SourceSupportError('Explicit context trim needs a finite straight original support leg')
            a, b = np.asarray(section.points)
            length = float(np.linalg.norm(b[:2] - a[:2]))
            if sum(trims) >= length - 1e-6:
                raise SourceSupportError('Explicit context trims consume the whole original support leg')
            points = (tuple(a + trims[0] / length * (b - a)), tuple(b - trims[1] / length * (b - a)))
            sections[i] = replace(section, points=points,
                selection_ends=tuple(sorted(set(section.selection_ends) | {end for end, t in enumerate(trims) if t})))
            furniture = [replace(p, station=p.station - trims[0]) if p.section_id == section.id else p for p in furniture]
            details.append(dict(section=section.id, original_controls=section.points,
                                active_controls=points, trims_m=trims))
        result = replace(network, sections=tuple(sections), furniture=tuple(furniture))
        report['context_trims'].extend(details)
        return result

    @staticmethod
    def _approaches(context, junction):
        topology = RouteTopology(context.network)
        from .intersections import Approach
        result = []
        for end in junction.ends:
            route = Route.build(context, topology.trace(end))
            context.products[route_key(end)] = route
            result.append(Approach(route, end.end, context[Key('profile', end.section_id)]))
        return sorted(result, key=lambda a: math.atan2(a.direction[1], a.direction[0]) % math.tau)

    def _lengths(self, junction, approaches):
        from .corner_radii import CornerRadiusPolicy
        # A style radius is a preferred maximum, not a reason to publish a long
        # outside tail. Recover finite support for its permitted fitting range;
        # the ordinary contact solver still chooses the largest fitting radius.
        # Strict overrides keep their complete requested reach.
        radii = [tuple(CornerRadiusPolicy.candidates(*junction.corner_radius(
            End(a.source.id, a.end))))[-1] for a in approaches]
        result = [a.width + 2. * r + self.clearance for a, r in zip(approaches, radii)]
        for i, a in enumerate(approaches):
            j = (i + 1) % len(approaches)
            b = approaches[j]
            gap = math.atan2(float(a.direction[0] * b.direction[1] - a.direction[1] * b.direction[0]),
                             float(a.direction @ b.direction)) % math.tau
            if gap <= 1e-6 or gap >= math.pi - 1e-6:
                continue
            ha, hb = a.width * .5 + abs(a.offset), b.width * .5 + abs(b.offset)
            radius = max(radii[i], radii[j])
            result[i] = max(result[i], (hb + radius + math.cos(gap) * (ha + radius)) / math.sin(gap) + self.clearance)
            result[j] = max(result[j], (ha + radius + math.cos(gap) * (hb + radius)) / math.sin(gap) + self.clearance)
        return {End(a.source.id, a.end): value for a, value in zip(approaches, result)}

    def _extension(self, section, end, cut, needed):
        node, opened, incidence = cut
        source = self.sources[section.id]
        if node['topology_kind'] == 'JUNCTION':
            raise SourceSupportError('Outside junction support requires its complete original incidence plan')
        if node['topology_kind'] not in ('CONTINUATION', 'ENDPOINT'):
            raise SourceSupportError('Outside source topology is not classified')
        if opened['owner_scope'] != 'ACTIVE':
            raise SourceSupportError('Retained physical boundary support needs coordinated ownership')
        if type(incidence['source_end']) is not int or incidence['source_end'] not in (0, 1):
            raise SourceSupportError('Outside source end has no valid direction')
        outside = incidence['source_end']
        nodes = tuple(map(str, opened['source_nodes']))
        if str(node['osm_node_id']) != nodes[outside]:
            raise SourceSupportError('Outside node differs from its directed original segment')
        incoming = [leg for leg in node['legs'] if tuple(map(str, leg['source_nodes'])) == nodes
                    and str(leg['source_way_id']) == str(opened['source_way_id'])]
        if len(incoming) != 1:
            raise SourceSupportError('Outside cut has no unique original incident segment')
        incoming, = incoming
        if (incoming['source_end'] != outside or str(incoming['other_node_id']) != nodes[1 - outside]):
            raise SourceSupportError('Original outside leg has conflicting node directions')
        if incoming['category'] != section.category or incoming['tags'] != source.get('tags', {}):
            raise SourceSupportError('Outside source dimensions cannot cross a category or tag change')
        cut_xy = np.asarray(incidence['clipped_point'], dtype=float)
        node_xy = np.asarray(node['position'], dtype=float)
        endpoints = [None, None]
        endpoints[outside] = node_xy
        endpoints[1 - outside] = np.asarray(incoming['other_position'], dtype=float)
        parameter = opened['parameter']
        if (not 0. < parameter < 1. or
                np.linalg.norm((1. - parameter) * endpoints[0] + parameter * endpoints[1] - cut_xy) > .001):
            raise SourceSupportError('Outside coordinates do not reproduce the exact original source cut')
        controls, intervals, remaining = [], [], min(needed, self.maximum_extension)
        original_budget = remaining
        legs = [(node_xy, nodes, float(opened['parameter']), float(outside))]
        onward = [leg for leg in node['legs'] if tuple(map(str, leg['source_nodes'])) != nodes]
        if node['topology_kind'] == 'CONTINUATION' and len(onward) == 1:
            leg, = onward
            if (str(leg['source_way_id']) == str(opened['source_way_id']) and
                    leg['category'] == section.category and leg['tags'] == source.get('tags', {}) and
                    str(leg['source_nodes'][leg['source_end']]) == str(node['osm_node_id'])):
                legs.append((np.asarray(leg['other_position'], dtype=float), tuple(map(str, leg['source_nodes'])),
                             float(leg['source_end']), float(1 - leg['source_end'])))
        current = cut_xy
        for target, ids, start, stop in legs:
            length = float(np.linalg.norm(target - current))
            if not math.isfinite(length) or length <= 1e-6:
                raise SourceSupportError('Original outside support has a zero-length segment')
            distance = min(remaining, length)
            fraction = distance / length
            point = current + fraction * (target - current)
            controls.append(tuple(map(float, point)) + (0.,))
            finish = start + fraction * (stop - start)
            intervals.append((str(opened['source_way_id']), *ids,
                              finish if end == 0 else start, start if end == 0 else finish))
            remaining -= distance
            if remaining <= 1e-8:
                break
            current = target
        if original_budget - remaining <= 1e-6:
            raise SourceSupportError('No additional original outside support is available')
        points = tuple(reversed(controls)) + section.points if end == 0 else section.points + tuple(controls)
        true_endpoint = (node['topology_kind'] == 'ENDPOINT' and len(controls) == 1 and
                         abs(intervals[0][3 if end == 0 else 4] - outside) <= 1e-10)
        selection = set(section.selection_ends)
        selection.discard(end) if true_endpoint else selection.add(end)
        return replace(section, points=points, selection_ends=tuple(sorted(selection))), dict(
            section=section.id, end=end, source_way=str(opened['source_way_id']),
            original_cut=section.points[0 if end == 0 else -1],
            extra_controls=controls, source_intervals=list(reversed(intervals)) if end == 0 else intervals,
            added_polyline_length_m=original_budget - remaining, requested_extra_m=needed,
            maximum_extension_m=self.maximum_extension,
            output_end_kind='SOURCE_ENDPOINT' if true_endpoint else 'SELECTION_CUT')

    def prepare(self, network):
        from .intersections import IntersectionProvider, CornerConstructionError, ApproachExtentError
        report = dict(policy=self.policy, status='NO_EXTENSION_NEEDED', selection_cuts=[],
                      context_trims=[], extensions=[], provenance_phase='INITIAL_IMPORT',
                      style_support='FIT_WITHIN_ALLOWED_RADIUS_RANGE',
                      scope='Original source/contact preparation; meshes unqualified')
        try:
            cuts = self._cuts(network)
            sections = tuple(replace(s, selection_ends=tuple(sorted(set(s.selection_ends) |
                             {end.end for end in cuts if end.section_id == s.id}))) for s in network.sections)
            network = replace(network, sections=sections)
            network = self._trim_context_legs(network, report)
        except ValueError as exc:
            return self._error(network, report, exc)
        report['selection_cuts'] = [dict(section=e.section_id, end=e.end) for e in sorted(cuts)]
        if not cuts or not network.intersections:
            return SourcePreparation(network, report)
        patterns = PatternLibrary(network.settings.pattern_cache_bytes)
        coordinator = UpdateCoordinator((PathProvider(), ProfileProvider()), patterns)
        attempted = set()
        try:
            update = coordinator.prepare(network)
            coordinator.commit(update)
            while True:
                context = BuildContext(network, dict(update.products), patterns)
                failed = None
                for junction in network.intersections:
                    if len(junction.ends) < 3:
                        continue
                    approaches = self._approaches(context, junction)
                    try:
                        IntersectionProvider.footprint(context, junction)
                    except (CornerConstructionError, ApproachExtentError) as exc:
                        failed = (junction, approaches, exc)
                        break
                if failed is None:
                    report['status'] = 'CONTACTS_READY'
                    return SourcePreparation(network, report)
                junction, approaches, error = failed
                required = self._lengths(junction, approaches)
                replacements, additions, prepends = {}, [], {}
                bundled = {sid for bundle in network.bundles for sid in bundle.sections}
                for approach in error.approaches:
                    if not isinstance(approach.path, Route):
                        continue
                    last = approach.path.pieces[-1]
                    end = End(last.section.id, int(last.forward))
                    if end not in cuts or end in attempted:
                        continue
                    attempted.add(end)
                    if last.section.id in bundled:
                        raise SourceSupportError('Bundle boundary extension needs a coordinated member support plan')
                    root = End(approach.source.id, approach.end)
                    target = max(required[root], getattr(error, 'required_stations', {}).get(root, 0.) + self.clearance)
                    extra = target - approach.length
                    if extra <= 1e-6:
                        continue
                    replacement, provenance = self._extension(last.section, end.end, cuts[end], extra)
                    if last.section.id in replacements:
                        raise SourceSupportError('A shared terminal Section needs coordinated support ownership')
                    replacements[last.section.id] = replacement
                    prepends[last.section.id] = len(provenance['extra_controls']) if end.end == 0 else 0
                    provenance['junction'] = junction.key.text()
                    additions.append(provenance)
                if not replacements:
                    return self._error(network, report, error)
                candidate = replace(network, sections=tuple(replacements.get(s.id, s) for s in network.sections))
                prepared = coordinator.prepare(candidate)
                furniture = []
                for placement in network.furniture:
                    if placement.section_id in replacements:
                        old = update.products[Key('path', placement.section_id)]
                        new = prepared.products[Key('path', placement.section_id)]
                        shift = new.knots[prepends[placement.section_id]]
                        parameter = old.parameter(placement.station) + shift
                        placement = replace(placement, station=float(new.station(parameter)))
                    furniture.append(placement)
                network = replace(candidate, furniture=tuple(furniture))
                coordinator.commit(prepared)
                update = prepared
                report['extensions'].extend(additions)
        except ValueError as exc:
            return self._error(network, report, exc)
