"""Source coupling admits real overlapping bodies without completed meshes."""
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import (Section, Network, Settings, Key, End,
    Traversal, Bundle, BundleMember, Intersection)
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.compounds import CompoundGroup, SampledFootprint, RimSegment, SharedClaim
from tested_python_streets.regions import Region, RegionUnion
from tested_python_streets.coupled_source_owners import (candidates, visible_intervals,
    analytic_body_overlap, side_footprint_overlap, prepare, _owners, _close, _close_steps, SourceSideTerminal)


class Line:
    def __init__(self, first, last):
        self.first, self.last = np.asarray(first), np.asarray(last)

    def __call__(self, t):
        return self.first+np.asarray(t)[..., None]*(self.last-self.first)

    def certified_parameters(self, sampler, first, last):
        return sampler.regular(first, last, sampler.spacing/np.linalg.norm(self.last-self.first))


def fixture(*, separated=False, reverse=False, angle=0., shift=(0., 0.)):
    rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
    def points(radius, first, last):
        theta = np.linspace(first, last, 9)
        xy = np.c_[radius*np.cos(theta), radius*np.sin(theta)] @ rotation.T+shift
        return tuple((x, y, 0.) for x, y in xy)
    inner = 90. if separated else 94.
    sections = [Section(1, 1, points(inner, 0., .8), width=8.),
                Section(2, 2, points(100., 0., .6), width=8.),
                Section(3, 2, points(100., .6, .8), width=8.)]
    # Original terminal owners each contain both actual terminal junctions.
    for sid, a, b in ((10, sections[0].points[0], sections[1].points[0]),
                      (20, sections[0].points[-1], sections[2].points[-1])):
        sections.append(Section(sid, sid, (a, b), width=4., interpolation=0))
    side = sections[1].points[-1]
    direction = np.array((np.cos(.6+angle), np.sin(.6+angle)))
    finish = tuple(np.r_[np.asarray(side[:2])+20.*direction, 0.])
    sections.append(Section(4, 4, (side, finish), width=4., interpolation=0))
    # Extra terminal incidences belong to actual distinct original Sections.
    for sid, root in ((11, sections[0].points[0]), (12, sections[1].points[0]),
                      (21, sections[0].points[-1]), (22, sections[2].points[-1])):
        end = tuple(np.asarray(root)+np.array((0., 10., 0.)))
        sections.append(Section(sid, sid, (root, end), width=4., interpolation=0))
    sites = (Intersection(100, (End(1, 0), End(10, 0), End(11, 0))),
        Intersection(101, (End(2, 0), End(10, 1), End(12, 0))),
        Intersection(200, (End(1, 1), End(20, 0), End(21, 0))),
        Intersection(201, (End(3, 1), End(20, 1), End(22, 0))),
        Intersection(300, (End(2, 1), End(3, 0), End(4, 0))))
    bundle = Bundle(7, (BundleMember((Traversal(1),)),
                        BundleMember((Traversal(2), Traversal(3)))))
    if reverse:
        sections = [replace(s, points=s.points[::-1]) for s in sections]
        sites = tuple(replace(site, ends=tuple(End(e.section_id, 1-e.end) for e in site.ends)) for site in sites)
        bundle = replace(bundle, members=tuple(BundleMember(tuple(Traversal(v.section_id, False)
            for v in member.path)) for member in bundle.members))
    source = Network('coupled-source-only', tuple(sections), intersections=sites, bundles=(bundle,),
        settings=Settings(curve_tolerance=.005))
    products = {}
    for section in source.sections:
        products[Key('path', section.id)] = Path(section)
        role = 0 if section.id == 1 else 1
        products[Key('profile', section.id)] = Crown(.08, 0., (1 if role == 0 else -1)*(-1 if reverse else 1))
    context = BuildContext(source, products, PatternLibrary())
    members = lambda *ids:tuple(Key('intersection', i) for i in ids)
    head = CompoundGroup(members(100, 101), (SharedClaim(10, *members(100, 101), 0., 0.),), (), ())
    tail = CompoundGroup(members(200, 201), (SharedClaim(20, *members(200, 201), 0., 0.),), (), ())
    center = np.asarray(side[:2]);tangent = np.array((-direction[1], direction[0]))
    polygon = np.array([center-3.*tangent-4.*direction, center+3.*tangent-4.*direction,
                        center+3.*tangent+4.*direction, center-3.*tangent+4.*direction])
    a, b = polygon[1]-polygon[0], polygon[2]-polygon[1]
    if a[0]*b[1]-a[1]*b[0] < 0.:polygon = polygon[::-1]
    key = members(300)[0]
    layout = SimpleNamespace(boundaries=tuple(Line(a, b) for a, b in zip(polygon, np.roll(polygon, -1, axis=0))),
                             approaches=(), mouths={})
    region = Region(key, polygon)
    footprint = SampledFootprint(region, tuple(RimSegment('curb', i, 0., 1.) for i in range(4)), layout)
    local = CompoundGroup((key,), (), (footprint,), RegionUnion().build((region,)))
    groups = (head, tail, local)
    layouts = {site.key:SimpleNamespace(mouths={}) for site in sites}
    layouts[key] = layout
    return context, bundle, groups, layouts


class CoupledSourceOwnersTests(unittest.TestCase):
    def test_real_member_suffix_and_distinct_terminal_owner_are_preserved(self):
        for reverse, angle, shift in ((False, 0., (0., 0.)), (True, .713, (1200., -1800.))):
            context, bundle, groups, layouts = fixture(reverse=reverse, angle=angle, shift=shift)
            saved = context.network.to_dict()
            requests = candidates(context, bundle, _owners(groups))
            self.assertEqual(len(requests), 1)
            request = requests[0]
            self.assertEqual(request.side, Key('intersection', 300))
            self.assertEqual(request.continuation, (Traversal(3, not reverse),))
            self.assertEqual(request.terminal, Key('intersection', 200))
            spans = visible_intervals(context, bundle, groups, layouts)
            body = analytic_body_overlap(context, bundle, spans)
            self.assertGreater(body['overlap_depth'], body['combined_source_allowance'])
            self.assertTrue(body['centerline_separation']['proof_checks'])
            side = side_footprint_overlap(context, request, spans, groups, layouts)
            self.assertEqual(side['opposite_section'], 1)
            self.assertFalse(side['completed_mesh_used'])
            self.assertEqual(context.network.to_dict(), saved)

    def test_pure_separation_owned_overlap_and_nonoverlapping_side_never_admit(self):
        context, bundle, groups, layouts = fixture(separated=True)
        spans = visible_intervals(context, bundle, groups, layouts)
        self.assertIsNone(analytic_body_overlap(context, bundle, spans))
        context, bundle, groups, layouts = fixture()
        request = candidates(context, bundle, _owners(groups))[0]
        # An owner consumes the entire only opposite member: no visible body
        # witness may be borrowed from its query-owned source strip.
        claimed = replace(groups[0], claims=(*groups[0].claims,
            SharedClaim(1, Key('intersection', 100), Key('intersection', 200), 0., 0.)))
        owned = (claimed, *groups[1:])
        spans = visible_intervals(context, bundle, owned, layouts)
        self.assertIsNone(analytic_body_overlap(context, bundle, spans))
        self.assertIsNone(side_footprint_overlap(context, request, spans, owned, layouts))
        context, bundle, groups, layouts = fixture(separated=True)
        request = candidates(context, bundle, _owners(groups))[0]
        self.assertIsNone(side_footprint_overlap(context, request,
            visible_intervals(context, bundle, groups, layouts), groups, layouts))

    def test_scope_selection_wrong_continuation_and_ambiguous_terminal_are_unavailable(self):
        context, bundle, groups, layouts = fixture()
        for changes in (dict(scope_ends=(1,)), dict(selection_ends=(1,))):
            network = replace(context.network, sections=tuple(replace(s, **changes) if s.id == 3 else s
                for s in context.network.sections), intersections=tuple(replace(site,
                    ends=tuple(end for end in site.ends if end != End(3, 1)))
                    if End(3, 1) in site.ends else site for site in context.network.intersections))
            other = BuildContext(network, context.products, context.patterns)
            self.assertFalse(candidates(other, bundle, _owners(groups)))
        changed = replace(bundle, members=(bundle.members[0], BundleMember((Traversal(3), Traversal(2)))))
        self.assertFalse(candidates(context, changed, _owners(groups)))
        self.assertFalse(candidates(context, bundle, _owners((groups[0], groups[2]))))

    def test_old_success_objects_and_explicit_prepared_bundles_are_untouched(self):
        context, bundle, groups, layouts = fixture(separated=True)
        with patch('tested_python_streets.coupled_source_owners._close_steps') as close:
            actual = prepare(context, context.network.intersections, layouts, groups)
            close.assert_not_called()
        self.assertEqual(actual, groups)
        for a, b in zip(actual, groups):self.assertIs(a, b)
        context, bundle, groups, layouts = fixture()
        with patch('tested_python_streets.coupled_source_owners.candidates') as find:
            self.assertEqual(prepare(context, context.network.intersections, layouts, groups,
                excluded_bundle_ids=(bundle.id,)), groups)
            find.assert_not_called()

    def test_only_typed_candidate_chart_unavailability_declines_before_admission(self):
        from tested_python_streets.bundle_axis import AxialChartUnavailable
        context, bundle, groups, layouts = fixture()
        with patch('tested_python_streets.bundle_angular_curves.AngularCrown',
                   side_effect=AxialChartUnavailable('original offset has no finite angular branch')), \
             patch('tested_python_streets.coupled_source_owners._close_steps') as close:
            result = prepare(context, context.network.intersections, layouts, groups)
            close.assert_not_called()
        for actual, old in zip(result, groups):self.assertIs(actual, old)
        failure = ValueError('invalid original analytic curb certificate')
        with patch('tested_python_streets.gapped_runtime_admission.visible_overlap', side_effect=failure):
            with self.assertRaises(ValueError) as raised:
                prepare(context, context.network.intersections, layouts, groups)
        self.assertIs(raised.exception, failure)

    def test_only_typed_expansion_adds_complete_existing_owners(self):
        from tested_python_streets.compound_preparation import CompoundExpansionRequired
        context, bundle, groups, layouts = fixture()
        owners = _owners(groups);calls = []
        result = SimpleNamespace(members=tuple(sorted(owners)))
        def exterior(keys):
            calls.append(keys)
            yield keys
            if len(calls) == 1:raise CompoundExpansionRequired((Key('intersection', 100),))
            return result
        planner = SimpleNamespace(context=context, sites=layouts, maximum_members=8, exterior_steps=exterior)
        steps = _close_steps(planner, {Key('intersection', 200), Key('intersection', 300)}, owners)
        self.assertEqual(next(steps), calls[0])
        self.assertEqual(next(steps), calls[1])
        with self.assertRaises(StopIteration) as completed:
            next(steps)
        got, history = completed.exception.value
        self.assertIs(got, result)
        self.assertEqual(calls[0], tuple(Key('intersection', i) for i in (200, 201, 300)))
        self.assertEqual(calls[1], tuple(sorted(owners)))
        self.assertEqual(len(history), 1)
        failure = ValueError('actual full-profile certificate failed')
        planner.exterior_steps = lambda keys:(_ for _ in ()).throw(failure)
        with self.assertRaises(ValueError) as raised:_close(planner, {Key('intersection', 300)}, owners)
        self.assertIs(raised.exception, failure)

    def test_sample_allocation_and_pair_search_are_bounded_before_work(self):
        context, bundle, groups, layouts = fixture()
        spans = visible_intervals(context, bundle, groups, layouts)
        calls = []
        original = Path.samples
        def samples(path, *args, **kwargs):
            calls.append(kwargs['max_samples'])
            return original(path, *args, **kwargs)
        with patch.object(Path, 'samples', samples):
            self.assertIsNotNone(analytic_body_overlap(context, bundle, spans))
        self.assertGreater(len(calls), 1)
        self.assertEqual(calls[0], 20000)
        self.assertTrue(all(a > b for a, b in zip(calls, calls[1:])))
        with self.assertRaisesRegex(ValueError, 'budget'):
            analytic_body_overlap(context, bundle, spans, maximum_samples=4)
        # A nonoverlapping first pair cannot launch an unbounded second search.
        with patch('tested_python_streets.gapped_runtime_admission.visible_overlap', return_value=None):
            with self.assertRaisesRegex(ValueError, 'original pair budget'):
                analytic_body_overlap(context, bundle, spans, maximum_pairs=1)

    def test_independent_terminal_couplings_rediscover_updated_owners_and_stop(self):
        context, bundle, groups, layouts = fixture()
        # This test isolates the owner-forest scheduler after the geometric
        # admission proofs, which the transformed real-curve test exercises.
        groups = (groups[0], replace(groups[1], physical_layout=SimpleNamespace(
            coupled_source_dependencies=(33,), coupled_source_owner_proofs=({'inherited':True},))), groups[2])
        other = Key('intersection', 301)
        groups = (*groups, CompoundGroup((other,), (), (), ()))
        layouts = dict(layouts)
        layouts[other] = SimpleNamespace(mouths={})
        first = SourceSideTerminal(bundle, Key('intersection', 300), Key('intersection', 200),
            (Traversal(3),), Traversal(1), 1, 32.)
        second = SourceSideTerminal(bundle, other, Key('intersection', 100),
            (Traversal(2, False),), Traversal(1), 0, 32.)
        discoveries, closures = [], []
        def discover(ctx, b, owners):
            discoveries.append(tuple(sorted(tuple(g.members) for g in set_by_identity(owners.values()))))
            return tuple(request for request in (first, second)
                if owners[request.side] is not owners[request.terminal])
        def close(planner, members, owners):
            keys = tuple(sorted(members));closures.append(keys)
            yield
            return CompoundGroup(keys, (), (), (), SimpleNamespace()), []
        with patch('tested_python_streets.coupled_source_owners.candidates', side_effect=discover), \
             patch('tested_python_streets.coupled_source_owners.visible_intervals', return_value=()), \
             patch('tested_python_streets.coupled_source_owners.side_footprint_overlap', return_value={'positive':True}), \
             patch('tested_python_streets.coupled_source_owners.analytic_body_overlap', return_value={'positive':True}), \
             patch('tested_python_streets.coupled_source_owners._close_steps', side_effect=close):
            result = prepare(context, context.network.intersections, layouts, groups)
        self.assertEqual(len(closures), 2)
        self.assertEqual(len(discoveries), 3)
        self.assertEqual([len(row) for row in discoveries], [4, 3, 2])
        self.assertEqual({frozenset(g.members) for g in result}, {
            frozenset(Key('intersection', i) for i in (200, 201, 300)),
            frozenset(Key('intersection', i) for i in (100, 101, 301))})
        self.assertEqual(sum(g.physical_layout.coupled_source_owner_proofs[-1]['original_owner_forest_joins']
                             for g in result), 2)
        old_tail = next(g for g in result if Key('intersection', 200) in g.members)
        self.assertEqual(old_tail.physical_layout.coupled_source_dependencies, (1, 2, 3, 33))
        self.assertEqual(old_tail.physical_layout.coupled_source_owner_proofs[0], {'inherited':True})

    def test_plan_policy_and_declared_sources_include_the_new_source_only_admission(self):
        from tested_python_streets import ownership
        from tested_python_streets.coupled_source_owners import COUPLED_OWNER_POLICY
        context, bundle, _, _ = fixture()
        task = next(task for task in ownership.JunctionOwnership.tasks(context.network)
                    if task.key == ownership.PLAN)
        self.assertEqual(task.version.split(':')[:2], [ownership.SOURCE_SAMPLING_POLICY, COUPLED_OWNER_POLICY])
        for sid in bundle.sections:
            for kind in ('path', 'profile'):self.assertIn(Key(kind, sid), task.dependencies)
        self.assertFalse(any(key.kind in ('road', 'junction_surface', 'bundle_plan', 'separator')
                             for key in task.dependencies))

    def test_remote_proof_sources_refresh_projected_products_and_unrelated_edits_reuse_them(self):
        from _python_street_fixtures import short_junction_pair
        from tested_python_streets import model, ownership
        from tested_python_streets.backend import PythonMeshBackend
        original = short_junction_pair(model)
        remote = Section(99, 99, ((500., 0., 0.), (520., 0., 0.)), width=4., interpolation=0)
        unrelated = replace(remote, id=98, street_id=98, points=((500., 100., 0.), (520., 100., 0.)))
        network = replace(original, sections=(*original.sections, remote, unrelated))
        tasks, build = ownership.JunctionOwnership.tasks, ownership.OwnershipPlan.build
        def source_tasks(source):
            for task in tasks(source):
                # The real coupling's Bundle sources are already declared.
                # This small source-preparer fixture declares its remote
                # witness explicitly while retaining every actual Task body.
                yield replace(task, dependencies=(*task.dependencies, Key('path', 99), Key('profile', 99))) \
                    if task.key == ownership.PLAN else task
        def source_plan(context, junctions):
            plan = build(context, junctions)
            group = plan.groups[0]
            if group.physical_layout is None:
                from tested_python_streets.compound_preparation import DeferredCompoundPlanner
                layouts = {j.key:context[ownership.footprint_key(j.key)] for j in junctions}
                group = DeferredCompoundPlanner(context, junctions, layouts).exterior(group.members)
                groups = (group, *plan.groups[1:])
                plan = ownership.OwnershipPlan(groups, {site:g for g in groups for site in g.members})
            group.physical_layout.coupled_source_dependencies = (99,)
            group.physical_layout.coupled_source_owner_proofs = ({'remote_width':
                context[Key('path', 99)].section.width},)
            return plan
        with patch.object(ownership.JunctionOwnership, 'tasks', side_effect=source_tasks), \
             patch.object(ownership.OwnershipPlan, 'build', side_effect=source_plan):
            engine = PythonMeshBackend(network.settings)
            first = engine.prepare(network);engine.commit(first)
            owner = first.products[ownership.PLAN].groups[0].members[0]
            surface, layout = ownership.surface_key(owner), Key('mouths', owner.id)
            changed = replace(network, sections=tuple(replace(s, width=s.width+.25) if s.id == 99 else s
                                                       for s in network.sections))
            incremental = engine.prepare(changed)
            fresh = PythonMeshBackend(changed.settings).prepare(changed)
            for key in (ownership.PLAN, layout, surface):
                self.assertIn(key, incremental.rebuilt)
                self.assertNotEqual(first.task_fingerprints[key], incremental.task_fingerprints[key])
            actual, expected = incremental.products[surface], fresh.products[surface]
            self.assertEqual(actual.pavement.metadata['coupled_source_owner_proofs'], ({'remote_width':4.25},))
            for name in ('physical_pavement', 'physical_grass'):
                a, b = getattr(actual, name), getattr(expected, name)
                if a is None:self.assertIsNone(b);continue
                np.testing.assert_array_equal(a.xyz, b.xyz)
                np.testing.assert_array_equal(a.quads, b.quads)
            engine.commit(incremental)
            other = replace(changed, sections=tuple(replace(s, width=s.width+.25) if s.id == 98 else s
                                                   for s in changed.sections))
            unchanged = engine.prepare(other)
            for key in (ownership.PLAN, layout, surface):
                self.assertNotIn(key, unchanged.rebuilt)
                self.assertIs(unchanged.products[key], incremental.products[key])


def set_by_identity(values):
    result = {}
    for value in values:result[id(value)] = value
    return result.values()


if __name__ == '__main__':unittest.main()
