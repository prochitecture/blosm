"""Visible Bundle separator runs following accepted junction curb chains."""
from dataclasses import dataclass
import numpy as np
from .mesh import Assembler
from .model import End, Key
from .ownership import PLAN, surface_key, CurbChainUnavailable
from .providers import Part


def compound_owners(context, bundle):
    plan = context.products.get(PLAN)
    if plan is None:
        return ()
    sections = set(bundle.sections)
    return tuple(group.members[0] for group in plan.groups if any(
        sid in sections for footprint in group.footprints for sid, end in footprint.layout.mouths))


def has_internal_junction(network, bundle):
    """A member's intermediate endpoint belongs to a multi-arm intersection."""
    internal = {End(t.section_id, 1 if t.forward else 0)
                for member in bundle.members for t in member.path[:-1]}
    return any(len(j.ends) > 2 and any(e in internal for e in j.ends) for j in network.intersections)


def has_consumed_sections(context, bundle):
    from .transitions import separator_range
    for member in bundle.members:
        for traversal in member.path:
            sid = traversal.section_id
            first, last = separator_range(context[Key('path', sid)], context[Key('road', sid)].metadata)
            if last <= first + 1e-7:
                return True
    return False


@dataclass(frozen=True)
class Terminal:
    endpoint: End
    owner: object
    attached: bool

    @property
    def token(self):
        return ('surface', self.owner.kind, self.owner.id) if self.owner is not None else ('attached' if self.attached else 'free',)


@dataclass
class Run:
    member: object
    terminals: tuple
    coordinates: tuple = ()


class CompoundBundleBuilder:
    @staticmethod
    def terminal(context, endpoint, station, metadata):
        covered = [span for span in (*metadata.get('compound_spans', ()), *metadata.get('junction_spans', ()))
                   if span['start'] - 1e-7 <= station <= span['end'] + 1e-7]
        owner = Key.parse(covered[0]['owner']) if covered else None
        attached = any(endpoint in j.ends for j in (*context.network.intersections, *context.network.connections))
        return Terminal(endpoint, owner, attached)

    @classmethod
    def runs(cls, context, member, role):
        from .bundles import Member, _join, _road_rim
        result, paths, rims, terminals, coordinates = [], [], [], [], []

        def finish():
            if paths:
                result.append(Run(Member(tuple(paths), _join(*rims),
                    (terminals[0].endpoint, terminals[-1].endpoint)), (terminals[0], terminals[-1]),
                    (coordinates[0], coordinates[-1])))
                paths.clear()
                rims.clear()
                terminals.clear()
                coordinates.clear()

        for traversal in member.path:
            sid = traversal.section_id
            path, road = context[Key('path', sid)], context[Key('road', sid)]
            from .transitions import separator_range
            a, b = separator_range(path, road.metadata)
            if b <= a + 1e-7:
                owners = {span['owner'] for span in (*road.metadata.get('compound_spans', ()),
                                                     *road.metadata.get('junction_spans', ()))}
                if not (terminals and terminals[-1].owner is not None and owners == {terminals[-1].owner.text()}):
                    finish()
                continue
            side = (1 if role == 0 else -1) * (1 if traversal.forward else -1)
            rim = _road_rim(context, sid, side)
            start, end = (a, b) if traversal.forward else (b, a)
            endpoints = (End(sid, 0 if traversal.forward else 1), End(sid, 1 if traversal.forward else 0))
            first = cls.terminal(context, endpoints[0], start, road.metadata)
            last = cls.terminal(context, endpoints[1], end, road.metadata)
            if not traversal.forward:
                rim = rim[::-1]
            if rims:
                previous = terminals[-1]
                connected = np.linalg.norm(rims[-1][-1] - rim[0]) <= 1e-7
                if previous.owner is not None or first.owner is not None:
                    connected = False
                    if previous.owner is not None and previous.owner == first.owner:
                        try:
                            bridge = context[surface_key(first.owner)].cap(
                                (previous.endpoint.section_id, previous.endpoint.end),
                                (first.endpoint.section_id, first.endpoint.end))
                        except CurbChainUnavailable:
                            pass  # Another mouth interrupts this median boundary.
                        else:
                            if max(np.linalg.norm(bridge[0] - rims[-1][-1]),
                                   np.linalg.norm(bridge[-1] - rim[0])) > 1e-7:
                                raise ValueError('Junction curb continuation does not meet its actual Bundle rims')
                            rims.append(bridge)
                            connected = True
                if not connected:
                    finish()
            from .bundles import CenterlineRange
            paths.append(CenterlineRange.build(path, a, b))
            rims.append(rim)
            if not terminals:
                terminals.append(first)
                coordinates.append(side)
            terminals.append(last)
            coordinates.append(side)
            # Defer an owned terminal until the next visible Section: the
            # accepted junction surface may supply an uninterrupted inner curb.
            if last.owner is None and ((b < path.length - 1e-7) if traversal.forward else (a > 1e-7)):
                finish()
        finish()
        return result

    @classmethod
    def regions(cls, context, bundle, runs):
        """Pair visible runs and retain an explicit crossing-to-free remainder."""
        if len(runs[0]) == len(runs[1]):
            return list(zip(*runs)), []
        from .bundle_regions import CrossingFreeContinuation, interior_crossings
        crossings = interior_crossings(context.network, bundle)
        longer = 0 if len(runs[0]) > len(runs[1]) else 1
        other = 1 - longer
        if len(runs[longer]) != len(runs[other]) + 1 or not runs[other]:
            raise ValueError('Unequal Bundle visibility needs an explicit region partition')
        plan = context.products.get(PLAN)
        def owner(site):
            return (plan.owner(site) if plan is not None else None) or site
        candidates = []
        for end in (0, 1):
            extra = runs[longer][0 if end == 0 else -1]
            terminal = runs[other][0 if end == 0 else -1].terminals[end]
            free = extra.terminals[end]
            if free.attached or free.owner is not None:
                continue
            for crossing in crossings:
                if (extra.terminals[1 - end].owner == owner(crossing.sites[longer]) and
                        terminal.owner == owner(crossing.sites[other])):
                    candidates.append((end, extra, CrossingFreeContinuation.build(
                        context.network, bundle, crossing, other, end)))
        if len(candidates) != 1:
            raise ValueError('Unequal Bundle events need one topology-proven crossing-to-free continuation')
        end, extra, continuation = candidates[0]
        boundary = cls.runs(context, continuation.member, other)
        if len(boundary) != 1 or boundary[0].terminals[end].attached:
            raise ValueError('Crossing continuation needs one complete visible run to its free boundary')
        paired = [list(side) for side in runs]
        paired[other].insert(0 if end == 0 else len(paired[other]), boundary[0])
        return list(zip(*paired)), [continuation.record()]

    @staticmethod
    def caps(context, first, second):
        from .bundles import BundleProvider, Cap
        caps = []
        for end, (a, b) in enumerate(zip(first.terminals, second.terminals)):
            q = first.coordinates[end], second.coordinates[end]
            if a.owner is None or a.owner != b.owner:
                caps.append(BundleProvider.cap(context, first.member, second.member, end, q))
                continue
            product = context[surface_key(a.owner)]
            ends = ((a.endpoint.section_id, a.endpoint.end), (b.endpoint.section_id, b.endpoint.end))
            outward = tuple(.5 * (v * (1 if e[1] == 0 else -1) + 1.) for v, e in zip(q, ends))
            curve = product.edge_cap(*ends, *outward)
            expected = (first.member.rim[0 if end == 0 else -1], second.member.rim[0 if end == 0 else -1])
            if max(np.linalg.norm(curve[i] - point) for i, point in ((0, expected[0]), (-1, expected[1]))) > 1e-7:
                raise ValueError('Compound Bundle cap does not meet its actual member rims')
            caps.append(Cap.rounded(curve))
        return caps

    @classmethod
    def build(cls, context, bundle):
        from .bundles import BundleProvider, CenterlineRange, certify_disjoint
        paths = []
        for member in bundle.members:
            visible = []
            for traversal in member.path:
                road = context[Key('road', traversal.section_id)]
                from .transitions import separator_range
                a, b = separator_range(context[Key('path', traversal.section_id)], road.metadata)
                if b > a + 1e-7:
                    visible.append(CenterlineRange.build(context[Key('path', traversal.section_id)], a, b))
            paths.append(tuple(visible))
        proof = certify_disjoint(paths[0], paths[1])
        profiles = [context[Key('profile', t.section_id)] for member in bundle.members for t in member.path]
        if any(abs(profile.datum - profiles[0].datum) > 1e-7 for profile in profiles):
            raise ValueError('Compound Bundle members require a common minimum height')
        runs = [cls.runs(context, member, role) for role, member in enumerate(bundle.members)]
        regions, continuations = cls.regions(context, bundle, runs)
        assembly, reports = Assembler(), []
        for first, second in regions:
            region_proof = certify_disjoint(first.member.path, second.member.path)
            datum = profiles[0].datum
            ceiling = datum + profiles[0].rise * context.network.settings.exaggeration
            if any(np.min(run.member.rim[:, 2]) < datum - 1e-7 or
                   np.max(run.member.rim[:, 2]) > ceiling + 1e-7 for run in (first, second)):
                raise ValueError('Separator continuation rim exceeds the original Bundle height envelope')
            caps = cls.caps(context, first, second)
            part = BundleProvider.region(context, bundle, (first.member, second.member), caps, proof + region_proof)
            assembly.add(part.mesh.xyz, part.mesh.quads)
            reports.append(part.metadata)
        if not reports:
            return None  # Every member span belongs to the compound pavement.
        mesh = assembly.mesh('grass')
        return Part(bundle.key, mesh, dict(quality=mesh.validate(), construction='junction-capped Bundle separator runs',
            runs=reports, derived_continuations=continuations,
            compound_owners=[owner.text() for owner in compound_owners(context, bundle)],
            physical_slopes=mesh.slope_statistics(context.network.settings.exaggeration)))
