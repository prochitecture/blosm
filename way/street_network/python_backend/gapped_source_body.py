"""Explicit body tiles around real, separated one-sided junctions.

Source inventory owns the topology. A near-side junction retains its complete
profiles and published curb; only its opposite original source band crosses
that station interval. Private common fibres create no source Ends or aliases.
There is deliberately no provider/discovery hook in this module.
"""
from types import SimpleNamespace
from dataclasses import replace
import numpy as np

from .bundle_angular_island_chain import MultiIslandSurfaceBuilder
from .bundle_angular_open_islands import clip_ray, gap_intervals
from .bundle_angular_source_strips import SourcePreservingAngularBridge
from .bundle_cap_sampling import SourceRimSampling
from .bundle_merged_body import MergedRunBody
from .bundle_merged_bridge import MergedProfileBridge
from .bundle_visible_ports import ProfileSamplingPatterns, VisibleCaps
from .bundles import BundleProvider, Member, Cap
from .drainage import BoundaryDrainage
from .fixed_rim_bands import BandedSeparatorChart, GradedRimSeparatorChart
from .four_side_charts import quality
from .gapped_source_inventory import GappedSourceInventory, GappedSourceUnavailable
from .mesh import Assembler, Mesh
from .model import Key
from .ownership import surface_key
from .profiles import plateau_extension
from .regions import signed_area
from .routes import RouteTopology
from .separator_surfaces import FixedRimSeparator
from .validation import validate_overlaps
from .bundle_merged_runs import assign_declared_samples
from .paired_source_body import PairedSourceBodyPiece
from .separator_ports import PORTS
from .updates import BuildContext
from .gapped_source_validation import complete_width
from .providers import CarriagewayProvider


class _Domain:
    def __init__(self, context, inventory, ports):
        self.context, self.bundle = context, inventory.bundle
        self.frame, self.field = inventory.frame, inventory.field
        self.rays, self.ports = self.field.curves, ports
        values = {end: [self.frame.local(port.coordinates(np.array((0., 1.)))[0])[:, 0]
                       for (role, at), port in ports.items() if at == end] for end in (0, 1)}
        self.head = max((float(row.max()) for row in values[0]), default=-np.inf)
        self.tail = min((float(row.min()) for row in values[1]), default=np.inf)

    def sides(self, role):
        return (self.rays[0], self.rays[1]) if role == 0 else (self.rays[3], self.rays[2])

    def patterns_for(self, port):
        return ProfileSamplingPatterns(self.context.patterns, port.mouth.parameters)


class _CapSampling(SourceRimSampling):
    """Keep actual unconsumed crossroad curb provenance for whole-edge audit."""
    def __init__(self, context):
        super().__init__(RouteTopology(context.network).attachments(context.network))
        self.context, self.sections = context, set()

    def __call__(self, context, section_id, side):
        points = super().__call__(self.context, section_id, side)
        self.sections.add(section_id)
        return points


def _ordered(points, keys, frame):
    if frame.local(points[[0, -1], :2])[0, 0] > frame.local(points[[0, -1], :2])[1, 0]:
        return points[::-1].copy(), list(keys)[::-1]
    return points.copy(), list(keys)


def _join(chains, *, closed=False):
    points, keys = [], []
    for xy, names in chains:
        xy = np.asarray(xy)
        if len(xy) != len(names) or not len(xy):
            raise ValueError('A gapped body chain lost its explicit shared identities')
        if points:
            if not np.array_equal(points[-1], xy[0]):
                raise ValueError('Gapped body chains changed an authoritative complete endpoint')
            if keys[-1] != names[0]:
                raise ValueError('Gapped body chains have different declared endpoint identities')
            points.extend(xy[1:]); keys.extend(names[1:])
        else:
            points.extend(xy); keys.extend(names)
    if closed:
        if not np.array_equal(points[-1], points[0]) or keys[-1] != keys[0]:
            raise ValueError('A gapped body disk has no exact closed physical boundary')
        points.pop(); keys.pop()
    return np.asarray(points), keys


class GappedSourceBodyPiece:
    """A bounded explicit consumer; ordinary successful providers stay first."""
    max_cells = 50000

    @staticmethod
    def _balance_bank_ends(bank, tolerance):
        """Redistribute a clipped sampling interval on a private smooth bank."""
        if len(bank) < 3:
            return None
        balanced = bank.copy()
        changed = False
        for row in (balanced, balanced[::-1]):
            a, b, c = row[:3, :2]
            short = np.linalg.norm(b-a)
            long = np.linalg.norm(c-b)
            if short <= 0. or long <= 8.*short:
                continue
            # The replacement point is on the longer original segment.
            # A monotone projection of the skipped corner onto the new chord
            # bounds the deviation in both directions over the full interval.
            middle = b + (.5*(long-short)/long)*(c-b)
            chord = middle-a
            length2 = chord @ chord
            if length2 <= 0.:
                continue
            fraction = (b-a) @ chord / length2
            if 0. < fraction < 1. and np.linalg.norm(b-a-fraction*chord) <= tolerance:
                row[1, :2] = middle
                changed = True
                if len(bank) < 5:
                    break
        return balanced if changed else None

    @staticmethod
    def make_body():
        from .split_merged_body import BoundedGappedRunBody
        return BoundedGappedRunBody()

    @staticmethod
    def _validate(context, inventory):
        """Recheck original identity and actual published interfaces at use."""
        for role, spans in enumerate(inventory.spans_by_role):
            original = tuple(inventory.bundle.members[role].path)
            for span in spans:
                visits = tuple(span.visits)
                starts = [i for i, visit in enumerate(original) if visit == visits[0]]
                if (span.role != role or len(starts) != 1
                        or original[starts[0]:starts[0]+len(visits)] != visits
                        or len(span.source_intervals) != len(visits)):
                    raise ValueError('Gapped source intervals changed their original ordered Bundle visits')
                for visit, piece, row in zip(visits, span.route.pieces, span.source_intervals):
                    path = context[Key('path', visit.section_id)]
                    if (piece.path is not path or piece.forward != visit.forward
                            or row['section'] != visit.section_id or row['forward'] != visit.forward
                            or not 0. <= row['first'] < row['last'] <= path.length):
                        raise ValueError('A gapped source interval changed its original finite Path')
                for end, node in enumerate((span.head, span.tail)):
                    if node.prepared_mouth is not None:
                        PairedSourceBodyPiece._mouth(context, node.prepared_mouth, end)
                        continue
                    if node.kind != 'PHYSICAL' or node.endpoint is None or node.registry_key is not None:
                        raise ValueError('A separated source endpoint needs its actual original physical incidence')
                    visit = visits[0 if end == 0 else -1]
                    port = node.port
                    if (port is None or port.path is not context[Key('path', visit.section_id)]
                            or (port.role, port.end, port.forward) != (role, end, visit.forward)
                            or node.source_station != port.station
                            or port.mouth.endpoint != node.endpoint or port.mouth.owner != node.owner
                            or context[PORTS].resolve(node.endpoint) is None
                            or not np.array_equal(node.physical_xyz, port.mouth.physical_xyz)):
                        raise ValueError('A gapped source endpoint changed its complete published original profile')
                    row = span.source_intervals[0 if end == 0 else -1]
                    station = row['first'] if (end == 0) == row['forward'] else row['last']
                    if station != port.station:
                        raise ValueError('A gapped source interval does not reach its actual complete profile')
        for gap in inventory.gaps:
            if (not any(gap.before is span for span in inventory.spans_by_role[gap.role])
                    or not any(gap.after is span for span in inventory.spans_by_role[gap.role])
                    or not any(gap.opposite is span for span in inventory.spans_by_role[1-gap.role])
                    or gap.bridge.owner != gap.before.tail.owner
                    or gap.bridge.owner != gap.after.head.owner
                    or not np.array_equal(gap.bridge.physical_xyz[[0, -1]],
                        np.array((gap.before.tail.physical_xyz[-1], gap.after.head.physical_xyz[-1])))):
                raise ValueError('A separated source gap changed its real near-side junction curb')

    def _bridge(self, context, inventory, endpoints, end, margin, label, cap=None):
        domain = _Domain(context, inventory, {(role, end): node.port for role, node in endpoints.items()})
        if any(node.kind != 'PHYSICAL' or node.port is None or node.prepared_mouth is not None
               for node in endpoints.values()):
            raise GappedSourceUnavailable('Separated source bridges require actual original physical profiles')
        result = {}
        for role, node in endpoints.items():
            bridge = SourcePreservingAngularBridge(domain)
            mesh, report = bridge.build(role, end, margin, physical_output=True)
            if cap is not None and role == 1:
                parity = sum(len(item['interface']['inner_ids'])-1 for item in result.values())
                if (parity+len(bridge.interfaces['inner_ids'])-1+len(cap.rim)-1) % 2:
                    mesh, report = bridge.build(role, end, margin, physical_output=True, extra_rows=1)
            interface = bridge.interfaces
            keys = [(label, role, 'bridge', i) for i in range(len(mesh.xyz))]
            ids = [interface['keys'].index(name) for name in interface['body_cut'][1]]
            inner = np.asarray(interface['inner_ids'], np.int32)
            result[role] = dict(mesh=mesh, keys=keys, interface=interface, back=np.asarray(ids, np.int32),
                cut=mesh.xyz[ids].copy(), inner=inner, node=node, report=report)
        return result

    @staticmethod
    def _cut(inventory, station, role, count):
        frame, field = inventory.frame, inventory.field
        edges, peaks = field.edges(station), field.peaks(station)
        y = np.linspace(edges[0 if role == 0 else 3], peaks[role], count+1)
        xy = frame.world(np.c_[np.full(len(y), station), y])
        # Exact Cartesian curb evaluation is the authoritative endpoint.
        xy[0] = field.curves[0 if role == 0 else 3](station)
        return np.c_[xy, field.values(xy)]

    @staticmethod
    def _attach(info, cut, names):
        if not np.array_equal(info['mesh'].xyz[info['back']], cut):
            raise ValueError('A source bridge changed its complete private common fibre')
        for node, name in zip(info['back'], names): info['keys'][int(node)] = name

    def _disk(self, context, inventory, chains, label):
        boundary, names = _join(chains, closed=True)
        if len(boundary) % 2:
            raise ValueError('A fixed source grass rim needs coordinated body-side parity')
        if signed_area(boundary[:, :2]) < 0.:
            boundary, names = boundary[::-1].copy(), names[::-1]
        try:
            surface = FixedRimSeparator(context.patterns, spacing=context.network.settings.spacing,
                max_cells=min(self.max_cells, 12000)).build(boundary,
                    datum=inventory.field.datum, height=inventory.field.rise, exaggeration=1.)
            mesh, metadata = surface.mesh, surface.metadata
        except ValueError as previous:
            failures = [str(previous)]
            candidates = [(BandedSeparatorChart, count) for count in (2, 4, 8, 16)]
            candidates.extend((GradedRimSeparatorChart, count) for count in (2, 4))
            for builder, bands in candidates:
                try:
                    chart = builder(context.patterns, spacing=context.network.settings.spacing,
                        max_cells=min(self.max_cells, 12000), bands=bands).build(boundary[:, :2])
                    mesh = chart.mesh
                    mesh.xyz[:, 2], extension = plateau_extension(mesh.xyz[:, :2], boundary,
                        inventory.field.datum, inventory.field.rise)
                    mesh.xyz[:len(boundary)] = boundary
                    mesh.material = 'grass'
                    stored = Mesh(mesh.xyz.astype(np.float32).astype(float), mesh.quads, 'grass')
                    if max(mesh.validate()['maximum_edge_ratio'], stored.validate()['maximum_edge_ratio']) > 8.:
                        raise ValueError('A private source grass disk needs better stored grading')
                    metadata = dict(chart=chart.metadata, height=extension, prior_failures=failures)
                    break
                except ValueError as error:
                    failures.append(str(error))
            else:
                raise ValueError('No accepted bounded source grass disk: '+str(failures)) from previous
        keys = [(label, 'grass', i) for i in range(len(mesh.xyz))]
        keys[:len(boundary)] = names
        return mesh, keys, dict(surface=metadata, physical_rim=len(boundary))

    def _opposite(self, context, inventory, role, head, tail, hk, tk, bank, label):
        settings = context.network.settings
        builder = MultiIslandSurfaceBuilder(context.patterns, spacing=settings.spacing,
            curve_tolerance=settings.curve_tolerance, max_cells=8000*settings.resolution)
        first, last = (float(inventory.frame.local(cut[:, :2])[:, 0].mean()) for cut in (head, tail))
        curve = clip_ray(inventory.field.curves[0 if role == 0 else 3], first, last, {})
        curve, kinds, names = builder.simplify_sample(curve, [('outer', role)]*(len(curve)-1),
            1e-9 if any(context[Key('path', v.section_id)].exact
                        for v in inventory.gaps[0].opposite.visits) else .6*settings.curve_tolerance,
            (label, 'opposite outer'))
        if max(np.linalg.norm(curve[0]-head[0, :2]), np.linalg.norm(curve[-1]-tail[0, :2])) > 1e-8:
            raise ValueError('An opposite source band changed its private outer endpoints')
        curve[[0, -1]] = head[0, :2], tail[0, :2]
        names[0], names[-1] = hk[0], tk[0]
        boundary = (np.vstack((head[::-1, :2], curve[1:], tail[1:, :2])),
            [('cut', 'head')]*(len(head)-1)+kinds+[('cut', 'tail')]*(len(tail)-1),
            hk[::-1]+names[1:]+tk[1:])
        seam = (bank[0][:, :2], bank[1])
        boundary = builder.balanced_physical(boundary, seam, str(label))
        mesh, keys, report = builder._direct_single(inventory.frame, inventory.field,
            boundary, seam, str(label)+'/body')
        prescribed = dict(zip(hk+tk, np.vstack((head, tail))))
        for i, key in enumerate(keys):
            if key in prescribed:
                if np.linalg.norm(mesh.xyz[i]-prescribed[key]) > 1e-8:
                    raise ValueError('An opposite source tile changed its complete common profile')
                mesh.xyz[i] = prescribed[key]
        return mesh, keys, report

    def build(self, context, inventory, *, namespace):
        if not isinstance(inventory, GappedSourceInventory):
            raise ValueError('A gapped body needs a certified original source inventory')
        hash(namespace)
        self._validate(context, inventory)
        source_context = context
        if context.network.settings.exaggeration != 1.:
            raise ValueError('Gapped source pieces require actual physical owner buffers at exaggeration one')
        authored_tolerance = context.network.settings.curve_tolerance
        # Source phase/Paths stay authoritative. Existing chart helpers get a
        # smaller construction allowance so their simplified chords leave the
        # independent complete-width proof its unchanged two-millimetre band.
        construction_tolerance = min(authored_tolerance, .0015)
        network = replace(context.network, settings=replace(context.network.settings,
            curve_tolerance=construction_tolerance))
        context = BuildContext(network, context.products, context.patterns)
        self.inventory = inventory
        frame, field, bundle = inventory.frame, inventory.field, inventory.bundle
        settings = context.network.settings
        paths = tuple(context[Key('path', sid)] for sid in inventory.source_sections)
        view = SimpleNamespace(frame_chart=frame, field=field, paths=paths, sign=1.)
        margin = 2.*settings.spacing
        endpoints = [{role: own[0].head for role, own in enumerate(inventory.spans_by_role)},
                     {role: own[-1].tail for role, own in enumerate(inventory.spans_by_role)}]
        outer, caps, merged = {}, {}, {}
        sampling = _CapSampling(source_context)
        members = tuple(Member((), np.empty((0, 3)), (endpoints[0][role].endpoint,
            endpoints[1][role].endpoint)) for role in (0, 1))
        for end in (0, 1):
            own = endpoints[end]
            if all(node.kind == 'PHYSICAL' for node in own.values()):
                cap = VisibleCaps.common_compound(context, *members, end)
                if cap is None: cap = BundleProvider.cap(context, *members, end, rim_source=sampling)
                if cap is None: raise GappedSourceUnavailable('Separated outer profiles need their actual original cap')
                caps[end] = cap
                outer[end] = self._bridge(context, inventory, own, end, margin,
                    (namespace, 'outer', end), cap)
            elif all(node.prepared_mouth is not None for node in own.values()):
                mouth = own[0].prepared_mouth
                if any(node.prepared_mouth.approach is not mouth.approach
                       or node.prepared_mouth.owner != mouth.owner for node in own.values()):
                    raise GappedSourceUnavailable('Outer merged profiles need one actual paired preparation')
                bridge = MergedProfileBridge()
                target = SimpleNamespace(frame=frame, field=field, source_sections=inventory.source_sections,
                    source_visits=inventory.source_visits, entrance=mouth.approach)
                mesh, report = bridge.build(context, mouth.approach, mouth.station, mouth.xyz,
                    target_inventory=target)
                halves = mouth.split_body_cut(bridge.body_cut,
                    crest_column=report.get('shared_crest_column'))
                merged[end] = dict(bridge=bridge, mesh=mesh, mouth=mouth, report=report,
                    halves=halves, keys=[(namespace, 'merged bridge', end, i) for i in range(len(mesh.xyz))])
            else:
                raise GappedSourceUnavailable('A mixed paired/separated outer interface needs its explicit policy')
        gap_bridges = []
        for index, gap in enumerate(inventory.gaps):
            left = self._bridge(context, inventory, {gap.role: gap.before.tail}, 1, margin,
                (namespace, 'gap', index, 'before'))[gap.role]
            right = self._bridge(context, inventory, {gap.role: gap.after.head}, 0, margin,
                (namespace, 'gap', index, 'after'))[gap.role]
            a, b = (float(frame.local(item['cut'][:, :2])[:, 0].mean()) for item in (left, right))
            regions, roots = gap_intervals(SimpleNamespace(frame=frame, field=field, rays=field.curves), a, b)
            if roots or any(not positive for first, last, positive in regions):
                raise GappedSourceUnavailable('Separated profile connectors cross a genuine overlap event')
            gap_bridges.append((left, right))
        boundaries = []
        for end in (0, 1):
            boundaries.append({role: (outer[end][role]['cut'] if end in outer else merged[end]['halves'][role])
                               for role in (0, 1)})
        gap_cuts = []
        for gap, (left, right) in zip(inventory.gaps, gap_bridges):
            rows = []
            for info in (left, right):
                station = float(frame.local(info['cut'][:, :2])[:, 0].mean())
                count = max(len(info['cut'])-1, 16*settings.resolution)
                rows.append({gap.role: info['cut'], 1-gap.role: self._cut(inventory, station, 1-gap.role, count)})
            gap_cuts.append(rows)
        starts = [boundaries[0]]+[rows[1] for rows in gap_cuts]
        finishes = [rows[0] for rows in gap_cuts]+[boundaries[1]]
        cores, report = [], dict(cores=[], gaps=[], outer=[], inventory=inventory.proofs)
        pavement, grass = Assembler(), Assembler()
        for index, (head, tail) in enumerate(zip(starts, finishes)):
            cuts = {(role, end): row[role] for end, row in enumerate((head, tail)) for role in (0, 1)}
            body = self.make_body()
            mesh, grassy, details = body.build_between(context, bundle, view, cuts,
                namespace=(namespace, 'core', index))
            pavement.add(mesh.xyz, mesh.quads, body.pavement_keys)
            if grassy is not None: grass.add(grassy.xyz, grassy.quads, body.grass_keys)
            cores.append(body); report['cores'].append(details)
        for end, infos in outer.items():
            core = cores[0 if end == 0 else -1]
            for role, info in infos.items():
                self._attach(info, core.cuts[role, end], core.cut_keys[role, end])
                pavement.add(info['mesh'].xyz, info['mesh'].quads, info['keys'])
            chains = []
            rim0 = (infos[0]['mesh'].xyz[infos[0]['inner']],
                    [infos[0]['keys'][i] for i in infos[0]['inner']])
            rim1 = (infos[1]['mesh'].xyz[infos[1]['inner']],
                    [infos[1]['keys'][i] for i in infos[1]['inner']])
            # Source bridges are ordered from their real mouth to body fibre.
            cut, keys = core.open_grass_cuts[end]
            flat = np.c_[cut, np.full(len(cut), field.ceiling)]
            if not np.array_equal(flat[0], rim0[0][-1]): flat, keys = flat[::-1], keys[::-1]
            cap = caps[end].rim.copy()
            if not np.array_equal(cap[0], rim1[0][0]): cap = cap[::-1]
            if not np.array_equal(cap[[0, -1]], np.array((rim1[0][0], rim0[0][0]))):
                raise ValueError('An actual outer cap changed its real complete source peaks')
            cap_keys = [(namespace, 'actual cap', end, i) for i in range(len(cap))]
            cap_keys[0], cap_keys[-1] = rim1[1][0], rim0[1][0]
            chains = [rim0, (flat, keys), (rim1[0][::-1], rim1[1][::-1]), (cap, cap_keys)]
            boundary, boundary_keys = _join(chains, closed=True)
            longitudinal = tuple(_ordered(info['mesh'].xyz[info['inner']],
                [info['keys'][i] for i in info['inner']], frame)[0] for info in infos.values())
            cap_members = tuple(Member((), rim, member.ends) for rim, member in zip(longitudinal, members))
            full_flat = flat if np.array_equal(flat[0], infos[0]['cut'][-1]) else flat[::-1]
            private = Cap(full_flat[:1], full_flat, full_flat[-1:], False)
            pair = [caps[end], private] if end == 0 else [private, caps[end]]
            try:
                surface = BundleProvider.region(context, bundle, cap_members, pair, 0)
                mesh = surface.mesh
                names = assign_declared_samples(mesh, boundary, boundary_keys, (namespace, 'cap', end))
                details = surface.metadata
            except ValueError:
                mesh, names, details = self._disk(context, inventory, chains, (namespace, 'cap', end))
            grass.add(mesh.xyz, mesh.quads, names)
            report['outer'].append(dict(end=end, connectors=[item['report'] for item in infos.values()], cap=details))
        for end, item in merged.items():
            bridge, mouth, mesh = item['bridge'], item['mouth'], item['mesh']
            core = cores[0 if end == 0 else -1]
            crest = item['report'].get('shared_crest_column', mouth.crest_column)
            columns = np.arange(len(bridge.back_ids))
            halves = (columns[:crest+1], columns[crest:][::-1])
            if mouth.approach.end_index == 1: halves = halves[::-1]
            for role, ids in enumerate(halves):
                if not np.array_equal(mesh.xyz[bridge.back_ids[ids]], core.cuts[role, end]):
                    raise ValueError('An actual merged source connector changed its complete body fibre')
                for node, key in zip(bridge.back_ids[ids], core.cut_keys[role, end]): item['keys'][int(node)] = key
            pavement.add(mesh.xyz, mesh.quads, item['keys'])
        for index, (gap, infos) in enumerate(zip(inventory.gaps, gap_bridges)):
            left, right = infos; before, after = cores[index:index+2]
            role, opposite = gap.role, 1-gap.role
            for info, core, end in ((left, before, 1), (right, after, 0)):
                self._attach(info, core.cuts[role, end], core.cut_keys[role, end])
                pavement.add(info['mesh'].xyz, info['mesh'].quads, info['keys'])
            head, tail = before.cuts[opposite, 1], after.cuts[opposite, 0]
            hk, tk = before.cut_keys[opposite, 1], after.cut_keys[opposite, 0]
            first, last = (float(frame.local(cut[:, :2])[:, 0].mean()) for cut in (head, tail))
            ray = field.curves[1 if opposite == 0 else 2]
            xy = clip_ray(ray, first, last, {})
            bank = np.c_[xy, np.full(len(xy), field.ceiling)]
            bank[[0, -1]] = head[-1], tail[-1]
            bk = [(namespace, 'opposite inner', index, i) for i in range(len(bank))]
            bk[0], bk[-1] = hk[-1], tk[-1]
            arim = _ordered(left['mesh'].xyz[left['inner']], [left['keys'][i] for i in left['inner']], frame)
            brim = _ordered(right['mesh'].xyz[right['inner']], [right['keys'][i] for i in right['inner']], frame)
            actual = gap.bridge.physical_xyz
            ak = [(namespace, 'actual junction curb', index, i) for i in range(len(actual))]
            ak[0], ak[-1] = arim[1][-1], brim[1][0]
            near = _join((arim, (actual, ak), brim))
            lc, lk = before.open_grass_cuts[1]; rc, rk = after.open_grass_cuts[0]
            lc, rc = (np.c_[cut, np.full(len(cut), field.ceiling)] for cut in (lc, rc))
            if not np.array_equal(rc[0], near[0][-1]): rc, rk = rc[::-1], rk[::-1]
            if not np.array_equal(lc[0], bank[0]): lc, lk = lc[::-1], lk[::-1]
            total = (len(near[0])-1)+(len(bank)-1)+(len(lc)-1)+(len(rc)-1)
            if total % 2:
                at = int(np.argmax(np.linalg.norm(np.diff(bank[:, :2], axis=0), axis=1)))
                bank = np.insert(bank, at+1, .5*(bank[at]+bank[at+1]), axis=0)
                bk.insert(at+1, (namespace, 'opposite inner parity', index))
            label = (namespace, 'opposite', index)
            try:
                mesh, names, details = self._opposite(context, inventory, opposite, head, tail,
                    hk, tk, (bank, bk), label)
            except ValueError:
                if any(context[Key('path', visit.section_id)].exact for visit in gap.opposite.visits):
                    raise
                tolerance = min(1e-6, .001*settings.curve_tolerance)
                balanced = self._balance_bank_ends(bank, tolerance)
                if balanced is None:
                    raise
                bank = balanced
                mesh, names, details = self._opposite(context, inventory, opposite, head, tail,
                    hk, tk, (bank, bk), label)
                details['shared_smooth_bank_sampling_tolerance'] = tolerance
            pavement.add(mesh.xyz, mesh.quads, names)
            disk, names, disk_report = self._disk(context, inventory,
                (near, (rc, rk), (bank[::-1], bk[::-1]), (lc, lk)), (namespace, 'junction grass', index))
            grass.add(disk.xyz, disk.quads, names)
            report['gaps'].append(dict(owner=gap.bridge.owner.text(), role=role,
                connectors=[info['report'] for info in infos], opposite=details, grass=disk_report))
        pave, grassy = pavement.mesh(), grass.mesh('grass')
        limit = self.max_cells*settings.resolution
        if len(pave.quads)+len(grassy.quads) > limit:
            raise ValueError('Gapped source bodies exceed their unchanged finite combined budget')
        self.pavement_keys = [None]*len(pave.xyz)
        self.grass_keys = [None]*len(grassy.xyz)
        for assembler, keys in ((pavement, self.pavement_keys), (grass, self.grass_keys)):
            for key, at in assembler.lookup.items(): keys[at] = key
        from .gapped_grass_recovery import recover
        grassy, self.grass_keys, recovery = recover(grassy, self.grass_keys, field,
            namespace=(namespace, 'assembled grass'), max_cells=limit-len(pave.quads))
        if recovery is not None:
            report['assembled_grass_recovery'] = recovery
        checks = []
        for mesh in (pave, grassy):
            mesh.validate(); double = quality(mesh)
            stored = Mesh(mesh.xyz.astype(np.float32).astype(float), mesh.quads, mesh.material)
            stored.validate(); published = quality(stored)
            if (max(double['maximum_edge_ratio'], published['maximum_edge_ratio']) > 8.
                    or (mesh.material != 'grass' and min(double['minimum_corner_sine'], published['minimum_corner_sine']) < .03)):
                raise ValueError('Gapped source tiles need better complete stored quad quality')
            slope, along = stored.slope_statistics()['maximum'], frame.longitudinal_bound(stored)
            if slope > .1+1e-6 or (mesh.material != 'grass' and along > .03+1e-8):
                raise ValueError('Gapped source tiles exceed the unchanged physical slope bounds')
            fixed = np.zeros(len(mesh.xyz), bool)
            fixed[np.unique(MergedRunBody._boundary_vertices(stored))] = True
            drainage = BoundaryDrainage(max_cells=limit).fill(stored, fixed)[1]
            if drainage['maximum_fill'] > 1e-6:
                raise ValueError('Gapped source tiles create an interior depression')
            checks.append(dict(material=mesh.material, physical_quality=double, stored_quality=published,
                maximum_slope=slope, longitudinal_bound=along, drainage=drainage))
        payloads = {Key('gapped_pave', 0): pave, Key('gapped_grass', 0): grassy}
        owners = {node.owner for own in inventory.spans_by_role for span in own for node in (span.head, span.tail)}
        owners.update(span.owner for span in inventory.owned_spans)
        # The overlap grid batches faces in insertion order. Use original
        # junction IDs to make its persisted report reproducible on reimport.
        for owner in sorted(owners):
            product = context[surface_key(owner)]
            payloads[Key('gapped_query_pave', owner.id)] = product.pavement.mesh
            if product.grass is not None: payloads[Key('gapped_query_grass', owner.id)] = product.grass.mesh
        for sid in sorted(sampling.sections):
            road = CarriagewayProvider.build(source_context, source_context.network.section_index[sid],
                sampling.attachments[sid])
            if road.mesh is None:
                raise ValueError('An original outer cap crossroad lost its actual published pavement')
            payloads[Key('gapped_query_crossroad', sid)] = road.mesh
        edge = lambda a, b: tuple(sorted((tuple(a), tuple(b))))
        pavement_edges = {edge(pave.xyz[a], pave.xyz[b]) for a, b in MergedRunBody._boundary_vertices(pave)}
        owner_edges = {edge(mesh.xyz[a], mesh.xyz[b])
                       for key, mesh in payloads.items() if key.kind in ('gapped_query_pave', 'gapped_query_crossroad')
                       for a, b in MergedRunBody._boundary_vertices(mesh)}
        shared, actual_curbs = 0, 0
        for a, b in MergedRunBody._boundary_vertices(grassy):
            token = edge(grassy.xyz[a], grassy.xyz[b])
            if token in pavement_edges: shared += 1
            elif token in owner_edges: actual_curbs += 1
            else: raise ValueError('A gapped source grass edge has no identical whole owned pavement edge')
        for mesh in (pave, grassy):
            rim = MergedRunBody._boundary_vertices(mesh)
            degree = np.bincount(rim.ravel(), minlength=len(mesh.xyz))
            if np.any(degree[np.unique(rim)] != 2):
                raise ValueError('Gapped source tiles contain an unresolved boundary branch')
        self.profile_ids = {}
        for role, spans in enumerate(inventory.spans_by_role):
            for i, span in enumerate(spans):
                for end, node in enumerate((span.head, span.tail)):
                    points = node.physical_xyz if node.prepared_mouth is None else node.prepared_mouth.xyz
                    ids = []
                    for point in points:
                        hits = np.flatnonzero(np.all(pave.xyz == point, axis=1))
                        if len(hits) != 1:
                            raise ValueError('Gapped source tiles lost an exact complete original physical profile')
                        ids.append(int(hits[0]))
                    self.profile_ids[role, i, end] = np.asarray(ids, np.int32)
        source_width = complete_width(source_context, inventory, pave,
            tuple(context[surface_key(owner)].pavement.mesh for owner in sorted(owners)))
        overlap = validate_overlaps(payloads, cell_size=2.)
        self.frame, self.field = frame, field
        report.update(construction='actual separated side-source tiles', stored_checks=checks,
            approximation_budget=dict(authored_curve_tolerance=authored_tolerance,
                existing_source_phase_unchanged=True, construction_curve_tolerance=construction_tolerance,
                complete_source_width_limit=.002),
            complete_source_width=source_width,
            total_quads=len(pave.quads)+len(grassy.quads), overlap=overlap,
            exact_source_profiles=len(self.profile_ids), shared_whole_grass_edges=shared,
            actual_owner_curb_edges=actual_curbs, boundary_manifold=True,
            actual_outer_crossroad_curbs=sorted(sampling.sections),
            new_source_aliases=0, junction_area_remeshed=False, runtime_discovery_enabled=False)
        return pave, grassy, report
