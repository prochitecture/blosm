"""Attached curved Bundle unions retaining actual profiled junction mouths.

The established axial/finite providers remain first. This branch restricts no
source controls, changes no corner radius and publishes ownership only after
the replacement pavement has passed its geometric guards.
"""
from dataclasses import replace
import numpy as np
from .bundle_axis import AxialChartUnavailable
from .bundle_angular_ports import AngularAttachedDomain
from .bundle_angular_attached_body import AngularAttachedBody
from .bundle_ends import AttachedAxialPlan
from .mesh import Mesh, Assembler
from .model import Key
from .providers import Part
from .surface_query import QuadHeightField
from .updates import BuildContext
from .validation import validate_overlaps


class AttachedAngularPlan:
    mode = 'union'
    transition_owners = ()
    max_cells = 100000
    requires_separator_query = True
    margin_factors = (2., 2.5, 3., 1.5, 1., .75, .5, 4.)

    @classmethod
    def build(cls, context, bundle):
        from .bundle_compounds import compound_owners
        if compound_owners(context, bundle):
            raise AxialChartUnavailable('Compound-owned mouths require their visible-region Bundle construction')
        domain = AngularAttachedDomain(context, bundle)
        from .bundle_cap_sampling import cap_parities
        domain.cap_parities = cap_parities(context, bundle)
        return cls.from_domain(context, bundle, domain)

    @staticmethod
    def cap(context, first, second, end, *, rim_source=None):
        from .bundles import BundleProvider
        return BundleProvider.cap(context, first, second, end, rim_source=rim_source)

    @classmethod
    def from_domain(cls, context, bundle, domain):
        failures = []
        for factor in cls.margin_factors:
            try:
                body = AngularAttachedBody(domain)
                pavement, grass, report = body.build(factor*context.network.settings.spacing)
                report['selection'] = dict(margin_factor=factor, prior_failures=failures)
                return cls(context, bundle, domain, body, pavement, grass, report)
            except ValueError as error:
                failures.append(dict(margin_factor=factor, reason=str(error)))
        raise ValueError('No qualified attached angular pavement: '+str(failures))

    def __init__(self, context, bundle, domain, body, pavement, grass, report):
        from .drainage import BoundaryDrainage
        self.bundle, self.settings = bundle, context.network.settings
        self.frame, self.field = domain.frame, domain.field
        self.datum, self.rise = self.field.datum, self.field.rise
        self.source_ids = tuple(visit.section_id for member in bundle.members for visit in member.path)
        self.proof_checks = domain.proof_checks
        self.ranges = {sid: tuple(sorted(port.station for port in domain.ports.values()
                                        if port.path.section.id == sid)) for sid in self.source_ids}
        self._bridge_rims = {key: interface['inner_rim'].copy()
                             for key, interface in body.bridge_interfaces.items()}
        self._cut_counts = {end: len(cut[0])-1 for end, cut in body.cut_interfaces.items()}
        edges = np.sort(np.stack((pavement.quads, np.roll(pavement.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        edges, counts = np.unique(edges, axis=0, return_counts=True)
        if counts.max() > 2:
            raise ValueError('Attached angular pavement contains a nonmanifold edge')
        boundary = edges[counts == 1]
        degrees = np.bincount(boundary.ravel(), minlength=len(pavement.xyz))
        if np.any(degrees[np.unique(boundary)] != 2):
            raise ValueError('Attached angular pavement has an unresolved boundary branch')
        fixed = np.zeros(len(pavement.xyz), bool)
        fixed[np.unique(boundary)] = True
        original = pavement.xyz[fixed].copy()
        pavement.xyz[:, 2], drainage = BoundaryDrainage(max_cells=self.max_cells).fill(pavement, fixed)
        if not np.array_equal(original, pavement.xyz[fixed]):
            raise ValueError('Attached angular drainage changed a prescribed boundary')
        self.maximum_slope = .1
        slopes = pavement.slope_statistics()
        longitudinal = self.frame.longitudinal_bound(pavement)
        if slopes['maximum'] > self.maximum_slope+1e-8 or longitudinal > .03+1e-8:
            raise ValueError('Attached angular pavement exceeded its physical slope limits')
        if len(pavement.quads)+len(grass.quads) > self.max_cells:
            raise ValueError('Attached angular Bundle exceeded its complete cell budget')
        self.query = QuadHeightField(pavement, datum=self.datum, ceiling=self.field.ceiling, max_cells=self.max_cells)
        self.grass_query = QuadHeightField(grass, datum=self.datum, ceiling=self.field.ceiling, max_cells=self.max_cells)
        self._grass_body = grass
        self._separator = None
        display = self._display(pavement)
        self._pavement = Part(bundle.key, display, dict(report,
            construction='attached angular common-fibre pavement union',
            source_sections=list(self.source_ids), attached_mouth_ranges=self.ranges,
            centerline_proof_checks=self.proof_checks, angular_certificate_lower=domain.curve_proofs,
            crown_ceiling=self.field.ceiling, physical_slopes=slopes,
            angular_longitudinal_bound=longitudinal, measured_slope_limit=self.maximum_slope,
            drainage=drainage, quality=display.validate(), resolution=self.settings.resolution,
            all_direct_quads=True))

    def _display(self, mesh):
        xyz = mesh.xyz.copy()
        xyz[:, 2] = self.datum+self.settings.exaggeration*(xyz[:, 2]-self.datum)
        return Mesh(xyz, mesh.quads.copy(), mesh.material)

    def pavement(self, context):
        return self._pavement

    def separator(self, context):
        if self._separator is not None:
            return self._separator
        from .bundles import BundleProvider, Member, Cap
        from .ownership import CompoundProduct, surface_key
        # Re-evaluate the actual rounded cap fields at physical exaggeration 1.
        # Road metadata supplies the same XY samples. This also preserves the
        # physical height query when the display exaggeration is zero.
        physical_network = replace(context.network, settings=replace(self.settings, exaggeration=1.))
        products = dict(context.products)
        for junction in context.network.intersections:
            if not any(end.section_id in self.source_ids for end in junction.ends):
                continue
            key = surface_key(junction.key)
            product = products.get(key)
            if isinstance(product, CompoundProduct):
                rims = {}
                for index, rim in product.rims.items():
                    heights = (product.field.datum+(rim[:, 2]-product.field.datum)/product.exaggeration
                               if product.exaggeration else product.field(rim[:, :2]))
                    rims[index] = np.column_stack((rim[:, :2], heights))
                products[key] = replace(product, rims=rims, exaggeration=1.)
        physical = BuildContext(physical_network, products, context.patterns)
        members = tuple(Member.build(physical, member, role) for role, member in enumerate(self.bundle.members))
        joins, caps = [], []
        for end in (0, 1):
            rims = [self._bridge_rims[role, end] if end == 0 else self._bridge_rims[role, end][::-1]
                    for role in (0, 1)]
            restricted = tuple(replace(member, rim=rim) for member, rim in zip(members, rims))
            flat = Cap.free(rims[0][-1 if end == 0 else 0], rims[1][-1 if end == 0 else 0],
                            self._cut_counts[end])
            actual = self.cap(physical, *members, end)
            if actual is None:
                raise ValueError('Attached angular separator lost its actual rounded junction cap')
            cap = BundleProvider.region(physical, self.bundle, restricted,
                [actual, flat] if end == 0 else [flat, actual], self.proof_checks)
            caps.append(cap)
            joins.append((flat.rim[0], flat.rim[-1], len(flat.rim)-1))
        assembly = Assembler()
        assembly.add(self._grass_body.xyz, self._grass_body.quads,
                     AttachedAxialPlan.seam_keys(self._grass_body, joins))
        for cap in caps:
            assembly.add(cap.mesh.xyz, cap.mesh.quads, AttachedAxialPlan.seam_keys(cap.mesh, joins))
        mesh = assembly.mesh('grass')
        edges = np.sort(np.stack((mesh.quads, np.roll(mesh.quads, -1, axis=1)), axis=-1).reshape(-1, 2), axis=1)
        unique, counts = np.unique(edges, axis=0, return_counts=True)
        islands = len(self._pavement.metadata['islands'])
        if counts.max() > 2 or len(mesh.xyz)-len(unique)+len(mesh.quads) != islands:
            raise ValueError('Attached angular separator caps do not form the source-derived islands')
        if len(mesh.quads)+len(self.query.quads) > self.max_cells:
            raise ValueError('Attached angular Bundle exceeded its complete cell budget')
        slopes = mesh.slope_statistics()
        if slopes['maximum'] > self.maximum_slope+1e-8:
            raise ValueError('Attached angular separator exceeded its physical slope limit')
        validate_overlaps({Key('separator', self.bundle.id): mesh}, cell_size=max(1., 2.*self.settings.spacing))
        self.grass_query = QuadHeightField(mesh, datum=self.datum, ceiling=self.field.ceiling, max_cells=self.max_cells)
        display = self._display(mesh)
        self._separator = Part(self.bundle.key, display, dict(construction='attached angular separator with original rounded caps',
            islands=islands, attached_caps=[True, True], caps=[cap.metadata for cap in caps],
            centerline_proof_checks=self.proof_checks, quality=display.validate(), physical_slopes=slopes,
            crown_ceiling=self.field.ceiling, all_direct_quads=True))
        self._grass_body = None
        return self._separator

    def height(self, xy):
        xy = np.asarray(xy, dtype=float)
        shape = xy.shape[:-1]
        points = xy.reshape(-1, 2)
        heights = np.maximum(self.query.evaluate(points, default=self.datum),
                             self.grass_query.evaluate(points, default=self.datum))
        return (self.datum+self.settings.exaggeration*(heights-self.datum)).reshape(shape)

    def release_geometry(self):
        self._grass_body = None
