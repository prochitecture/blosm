"""Curved Bundle unions partitioned at declared transverse Street crossings.

Derived runs keep the original Section identities and prepared physical mouths.
Existing separated construction remains authoritative when no run overlaps.
"""
from dataclasses import replace
import numpy as np
from . import model
from .bundle_axis import AxialChartUnavailable
from .bundle_angular_ports import AngularAttachedDomain, AngularPort
from .bundle_angular_curves import AngularCoordinates
from .bundle_angular_open_islands import gap_intervals
from .bundle_angular_attached import AttachedAngularPlan
from .bundle_cap_sampling import SourceRimSampling
from .bundle_regions import interior_crossings
from .bundle_visible_runs import VisibleBundleRuns
from .bundle_visible_ports import VisibleMouth, VisibleCaps, ProfileSamplingPatterns
from .bundles import BundleProvider, Member
from .mesh import Assembler
from .providers import Part
from .routes import Route, Piece, RouteTopology
from .separator_ports import PORTS


def shared_frame(context, bundle):
    """Choose coordinates from authored controls without refitting any path."""
    routes = []
    for member in bundle.members:
        pieces = [Piece(context[model.Key('path', visit.section_id)], visit.forward) for visit in member.path]
        route = Route(pieces)
        points = [piece.path.points if piece.forward else piece.path.points[::-1] for piece in pieces]
        route.points = np.concatenate([points[0], *(part[1:] for part in points[1:])])
        route.maximum_reach = max(abs(piece.section.offset) + .5*piece.section.width for piece in pieces)
        routes.append(route)
    return AngularCoordinates(routes)


class SampledAngularPort(AngularPort):
    def __init__(self, role, end, path, forward, profile, mouth):
        super().__init__(role, end, path, forward, mouth.station, len(mouth.parameters)-1, profile)
        self.mouth = mouth

    def coordinates(self, v):
        xy, tangent, heights = super().coordinates(v)
        xy, heights = np.array(xy, copy=True), np.array(heights, copy=True)
        values = np.asarray(v).reshape(-1)
        parameters = self.mouth.parameters
        indices = np.clip(np.searchsorted(parameters, values), 0, len(parameters)-1)
        exact = parameters[indices] == values
        prescribed = self.mouth.physical_xyz[indices[exact]]
        direct = np.c_[xy.reshape(-1, 2)[exact], heights.reshape(-1)[exact]]
        if len(prescribed) and np.max(abs(direct-prescribed)) > 1e-8:
            raise ValueError('Prepared compound profile disagrees with its source carriageway')
        xy.reshape(-1, 2)[exact] = prescribed[:, :2]
        heights.reshape(-1)[exact] = prescribed[:, 2]
        return xy, tangent, heights


class VisibleAngularDomain(AngularAttachedDomain):
    def __init__(self, context, bundle, frame):
        self._frame = frame
        super().__init__(context, bundle)

    def make_frame(self, paths):
        return self._frame

    def make_port(self, context, sites, role, end, visit, path, profile):
        endpoint = model.End(visit.section_id, end if visit.forward else 1-end)
        mouth = VisibleMouth.from_compound(context, endpoint, role, visit.forward)
        return (SampledAngularPort(role, end, path, visit.forward, profile, mouth) if mouth is not None
                else super().make_port(context, sites, role, end, visit, path, profile))

    def patterns_for(self, port):
        return (ProfileSamplingPatterns(self.context.patterns, port.mouth.parameters)
                if isinstance(port, SampledAngularPort) else self.context.patterns)


class VisibleAngularPlan(AttachedAngularPlan):
    @staticmethod
    def cap(context, first, second, end, *, rim_source=None):
        actual = VisibleCaps.common_compound(context, first, second, end)
        return (actual if actual is not None else
                BundleProvider.cap(context, first, second, end, rim_source=rim_source))

    @classmethod
    def build_domain(cls, context, bundle, domain):
        members = []
        for member in bundle.members:
            first, last = member.path[0], member.path[-1]
            ends = (model.End(first.section_id, 0 if first.forward else 1),
                    model.End(last.section_id, 1 if last.forward else 0))
            members.append(Member((), np.empty((0, 3)), ends))
        sampling = SourceRimSampling(RouteTopology(context.network).attachments(context.network))
        caps = [cls.cap(context, *members, end, rim_source=sampling) for end in (0, 1)]
        if any(cap is None for cap in caps):
            raise AxialChartUnavailable('An internal angular run lost its actual physical cap')
        domain.cap_parities = tuple((len(cap.rim)-1) % 2 for cap in caps)
        plan = cls.from_domain(context, bundle, domain)
        plan._pavement.metadata['actual_mouths'] = [dict(role=role, end=end,
            section=port.mouth.endpoint.section_id, source_end=port.mouth.endpoint.end,
            owner=port.mouth.owner.text(), station=port.mouth.station,
            parameters=port.mouth.parameters.tolist()) for (role, end), port in domain.ports.items()
            if isinstance(port, SampledAngularPort)]
        return plan


class InteriorAngularPlan:
    mode = 'union'
    transition_owners = ()
    requires_separator_query = True
    max_cells = AttachedAngularPlan.max_cells

    @classmethod
    def build(cls, context, bundle):
        crossings = interior_crossings(context.network, bundle)
        runs = VisibleBundleRuns.build(context.network, bundle, crossings)
        if not crossings or any(len(member.visits) != 1 for run in runs for member in run.members):
            raise AxialChartUnavailable('Internal angular runs need their multi-Section entrance adapter')
        if any(context[PORTS].resolve(endpoint) is None
               for run in runs for member in run.members for endpoint in member.ends):
            raise AxialChartUnavailable('Internal angular runs need actual attached mouths at both ends')
        frame = shared_frame(context, bundle)
        prepared, overlapping = [], False
        for run in runs:
            derived = replace(bundle, members=tuple(model.BundleMember(member.visits) for member in run.members),
                              interior_section_ids=())
            domain = VisibleAngularDomain(context, derived, frame)
            try:
                regions, roots = gap_intervals(domain, domain.head, domain.tail)
            except ValueError as error:
                raise AxialChartUnavailable(str(error)) from error
            overlapping |= bool(roots) or any(not positive for a, b, positive in regions)
            prepared.append((run, derived, domain))
        if not overlapping:
            return None
        children = tuple((run, VisibleAngularPlan.build_domain(context, derived, domain))
                         for run, derived, domain in prepared)
        return cls(context, bundle, children)

    def __init__(self, context, bundle, children):
        self.bundle, self.children = bundle, children
        self.source_ids = tuple(visit.section_id for member in bundle.members for visit in member.path)
        self._separator = None
        assembly = Assembler()
        for run, plan in children:
            mesh = plan.pavement(context).mesh
            assembly.add(mesh.xyz, mesh.quads)
        mesh = assembly.mesh()
        if len(mesh.quads) > self.max_cells:
            raise ValueError('Internal angular Bundle exceeds its complete cell budget')
        self._pavement = Part(bundle.key, mesh, dict(
            construction='visible angular pavement unions at original interior crossings',
            source_sections=list(self.source_ids), quality=mesh.validate(),
            runs=[dict(source_sections=run.source_sections, boundaries=run.boundaries,
                       pavement=plan.pavement(context).metadata) for run, plan in children],
            combined_cell_budget=self.max_cells, all_direct_quads=True))

    def pavement(self, context):
        return self._pavement

    def separator(self, context):
        if self._separator is None:
            assembly, reports = Assembler(), []
            for run, plan in self.children:
                part = plan.separator(context)
                assembly.add(part.mesh.xyz, part.mesh.quads)
                reports.append(dict(source_sections=run.source_sections, boundaries=run.boundaries,
                                    separator=part.metadata))
            mesh = assembly.mesh('grass')
            if len(mesh.quads) + sum(len(plan.query.quads) for run, plan in self.children) > self.max_cells:
                raise ValueError('Internal angular Bundle exceeds its complete cell budget')
            self._separator = Part(self.bundle.key, mesh, dict(
                construction='original-junction-capped angular separator runs',
                runs=reports, quality=mesh.validate(), all_direct_quads=True))
        return self._separator

    def height(self, xy):
        return np.maximum.reduce([plan.height(xy) for run, plan in self.children])

    def release_geometry(self):
        for run, plan in self.children:
            plan.release_geometry()
