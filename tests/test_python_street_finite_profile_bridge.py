"""Finite local profiles meet an independently certified original run field."""
from types import SimpleNamespace
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_terminal_ports import fixture
from tested_python_streets.finite_terminal_entrances import FinitePairedEntrance
from tested_python_streets.bundle_angular_curves import AngularCrown
from tested_python_streets.merged_run_inventory import inventory_frame
from tested_python_streets.bundle_merged_bridge import MergedProfileBridge
from tested_python_streets.mesh import Mesh
from tested_python_streets.uv import apply_mapping,PlanarMapping
from tested_python_streets import model
from tested_python_streets.paths import Path
from tested_python_streets.updates import BuildContext
from tested_python_streets.one_sided_paired_junctions import InteriorPairedEntrance
from tested_python_streets.bundle_visible_runs import BundleRun,MemberRun


def prepared(mode=3,reverse=False,reflection=False,end=0):
    context,pair,_=fixture(mode,reverse,reflection,.713 if reflection else 0.,(1800.,-1200.) if reverse else (0.,0.))
    approach=FinitePairedEntrance(context,pair.run,end,source_tolerance=.00025)
    frame=inventory_frame(context,context.network.bundles[0],approach)
    field=AngularCrown(frame,[ray.points for ray in approach.field.curves],datum=approach.field.datum,rise=approach.field.rise)
    inventory=SimpleNamespace(frame=frame,field=field,source_sections=approach.source_run_sections,entrance=approach)
    station=approach.minimum_station+5.
    approach.profile=approach.cross_section(station)
    u=approach.mouth_parameters(4,1.)[::-1]
    xy=approach.point(station,approach.width*(u-.5))
    mouth=np.c_[xy,approach.field.values(xy)]
    return context,approach,station,mouth,inventory


class FiniteProfileBridgeTests(unittest.TestCase):
    def test_terminal_incident_prefix_or_suffix_retains_the_whole_original_run_progress_proof(self):
        from tested_python_streets.profiles import Crown
        for end in (0, 1):
            with self.subTest(end=end):
                c, pair, _ = fixture()
                extra = tuple(model.Section(sid, sid, ((100., y, 0.), (150., y, 0.)), width=10.8)
                              for sid, y in ((3, 5.), (4, -3.)))
                members = tuple(model.BundleMember((m.visits[0], model.Traversal(3+role)))
                                for role, m in enumerate(pair.run.members))
                bundle = replace(c.network.bundles[0], members=members)
                seams = tuple(model.Connection(300+role, (model.End(1+role, 1), model.End(3+role, 0)))
                              for role in (0, 1))
                n = replace(c.network, sections=(*c.network.sections, *extra), bundles=(bundle,), connections=seams)
                products = dict(c.products)
                for role, section in enumerate(extra):
                    products[model.Key('path', section.id)] = Path(section)
                    products[model.Key('profile', section.id)] = Crown(.08, 0., 1 if role == 0 else -1)
                c = BuildContext(n, products, c.patterns)
                run = BundleRun(tuple(MemberRun(m.path, (model.End(1+role, 0), model.End(3+role, 1)),
                                    pair.run.boundaries) for role, m in enumerate(members)), pair.run.boundaries)
                approach = FinitePairedEntrance(c, run, end, source_tolerance=.00025)
                frame = inventory_frame(c, bundle, approach)
                field = AngularCrown(frame, [ray.points for ray in approach.field.curves], datum=0., rise=.08)
                selected = tuple((m.path[0 if end == 0 else -1],) for m in members)
                inventory = SimpleNamespace(frame=frame, field=field, entrance=approach,
                    source_sections=tuple(chain[0].section_id for chain in selected), source_visits=selected)
                station = approach.minimum_station+5.
                approach.profile = approach.cross_section(station)
                u = approach.mouth_parameters(4, 1.)[::-1]
                xy = approach.point(station, approach.width*(u-.5))
                mouth = np.c_[xy, approach.field.values(xy)]
                bridge = MergedProfileBridge()
                mesh, report = bridge.build(c, approach, station, mouth, target_inventory=inventory)
                np.testing.assert_array_equal(mesh.xyz[bridge.front_ids], mouth)
                proof = report['original_terminal_run_piece']
                self.assertEqual(set(proof['original_run_sections']), {1, 2, 3, 4})
                self.assertEqual(set(proof['piece_sections']), set(inventory.source_sections))
                self.assertTrue(proof['whole_original_run_progress_certified'])
                self.assertTrue(proof['requires_remaining_run_ownership_audit'])
                wrong = SimpleNamespace(**vars(inventory))
                wrong.source_visits = tuple((m.path[-1 if end == 0 else 0],) for m in members)
                wrong.source_sections = tuple(chain[0].section_id for chain in wrong.source_visits)
                with self.assertRaisesRegex(ValueError, 'exact incident original prefix or suffix'):
                    MergedProfileBridge().build(c, approach, station, mouth, target_inventory=wrong)
                wrong = SimpleNamespace(**vars(inventory));del wrong.entrance
                with self.assertRaisesRegex(ValueError, 'complete original source sections'):
                    MergedProfileBridge().build(c, approach, station, mouth, target_inventory=wrong)

    def test_opposite_role_interior_cut_retains_local_subset_of_a_longer_original_piece(self):
        c,pair,_=fixture()
        continuation=model.Section(3,1,((100.,5.,0.),(150.,5.,0.)),width=10.8)
        stem=model.Section(103,103,((100.,-3.,0.),(100.,-30.,0.)))
        bundle=replace(c.network.bundles[0],members=(model.BundleMember((model.Traversal(1),model.Traversal(3))),
            c.network.bundles[0].members[1]))
        n=replace(c.network,sections=(*c.network.sections,continuation,stem),bundles=(bundle,),
            connections=(model.Connection(301,(model.End(1,1),model.End(3,0))),),
            intersections=(*c.network.intersections,model.Intersection(202,(model.End(2,1),model.End(103,0)))))
        products=dict(c.products)
        products[model.Key('path',3)]=Path(continuation);products[model.Key('path',103)]=Path(stem)
        c=BuildContext(n,products,c.patterns)
        a=InteriorPairedEntrance(c,model.Traversal(2),model.Traversal(1),1,1)
        tokens=pair.run.boundaries
        run=BundleRun((MemberRun(bundle.members[0].path,(model.End(1,0),model.End(3,1)),tokens),pair.run.members[1]),tokens)
        frame=inventory_frame(c,bundle,SimpleNamespace(run=run,source_ranges=True))
        field=AngularCrown(frame,[r.points for r in a.field.curves],datum=a.field.datum,rise=a.field.rise)
        inventory=SimpleNamespace(frame=frame,field=field,source_sections=(1,3,2),
            source_visits=(bundle.members[0].path,bundle.members[1].path),entrance=a)
        at=a.minimum_station+5.;a.profile=a.cross_section(at)
        u=a.mouth_parameters(4,1.)[::-1];xy=a.point(at,a.width*(u-.5));mouth=np.c_[xy,a.field.values(xy)]
        bridge=MergedProfileBridge();mesh,report=bridge.build(c,a,at,mouth,target_inventory=inventory)
        np.testing.assert_array_equal(mesh.xyz[bridge.front_ids],mouth)
        self.assertEqual([x.original_incidence for x in a.source_mouths(at)],[False,True])
        self.assertEqual(set(a.source_run_sections),{1,2})
        self.assertEqual(set(inventory.source_sections),{1,2,3})
        self.assertGreater(report['whole_original_run_progress_lower'],0.)

    def test_interior_cut_uses_explicit_original_visits_and_keeps_the_derived_alias_nonincidental(self):
        c,pair,_=fixture()
        sections=[]
        for source in c.network.sections:
            if source.id==2:source=replace(source,points=((-50.,-4.,0.),*source.points))
            elif source.id==101:source=replace(source,points=((-50.,-4.,0.),(-50.,-34.,0.)))
            sections.append(source)
        n=replace(c.network,sections=tuple(sections));products=dict(c.products)
        for source in sections:products[model.Key('path',source.id)]=Path(source)
        c=BuildContext(n,products,c.patterns)
        a=InteriorPairedEntrance(c,model.Traversal(1),model.Traversal(2),0,0)
        frame=inventory_frame(c,n.bundles[0],SimpleNamespace(run=pair.run,source_ranges=True))
        field=AngularCrown(frame,[ray.points for ray in a.field.curves],datum=a.field.datum,rise=a.field.rise)
        inventory=SimpleNamespace(frame=frame,field=field,source_sections=(1,2),
            source_visits=((model.Traversal(1),),(model.Traversal(2),)),entrance=a)
        at=a.minimum_station+5.;a.profile=a.cross_section(at)
        u=a.mouth_parameters(4,1.)[::-1];xy=a.point(at,a.width*(u-.5));mouth=np.c_[xy,a.field.values(xy)]
        bridge=MergedProfileBridge();mesh,report=bridge.build(c,a,at,mouth,target_inventory=inventory)
        np.testing.assert_array_equal(mesh.xyz[bridge.front_ids],mouth)
        aliases=a.source_mouths(at)
        self.assertEqual([alias.original_incidence for alias in aliases],[True,False])
        self.assertGreater(aliases[1].station,0.)
        self.assertFalse(hasattr(a,'run'))
        self.assertGreater(report['whole_original_run_progress_lower'],0.)
        wrong=SimpleNamespace(**vars(inventory));wrong.source_visits=wrong.source_visits[::-1]
        with self.assertRaisesRegex(ValueError,'visit or role'):
            MergedProfileBridge().build(c,a,at,mouth,target_inventory=wrong)
        wrong=SimpleNamespace(**vars(inventory));del wrong.source_visits
        with self.assertRaisesRegex(ValueError,'explicit original visits'):
            MergedProfileBridge().build(c,a,at,mouth,target_inventory=wrong)

    def test_exact_and_smooth_both_source_directions_keep_the_authoritative_shared_crest(self):
        for mode in (0,3):
            for reverse,reflection,end in ((False,False,0),(True,True,1)):
                with self.subTest(mode=mode,reverse=reverse,end=end):
                    c,a,s,mouth,inventory=prepared(mode,reverse,reflection,end)
                    before=c.network.to_dict();bridge=MergedProfileBridge()
                    mesh,report=bridge.build(c,a,s,mouth,target_inventory=inventory)
                    np.testing.assert_array_equal(mesh.xyz[bridge.front_ids],mouth)
                    np.testing.assert_array_equal(bridge.body_cut,mesh.xyz[bridge.back_ids])
                    np.testing.assert_allclose(bridge.body_cut[:,2],inventory.field.values(bridge.body_cut[:,:2]),atol=1e-12)
                    crest=report['shared_crest_column']
                    self.assertAlmostEqual(bridge.body_cut[crest,2],a.field.ceiling,places=10)
                    stored=Mesh(mesh.xyz.astype(np.float32).astype(float),mesh.quads)
                    self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)
                    self.assertGreaterEqual(report['published_quality']['minimum_corner_sine'],.03)
                    self.assertLessEqual(report['physical_slope'],.1+1e-6)
                    self.assertLessEqual(max(report['entrance_bounds']),.03+1e-8)
                    self.assertLessEqual(report['source_containment']['maximum_combined_error'],.1*c.network.settings.curve_tolerance)
                    self.assertGreater(report['source_curve_identity']['checks'],0)
                    self.assertTrue(np.isfinite(apply_mapping(mesh,PlanarMapping()).uv_layers['UVMap']).all())
                    self.assertEqual(c.network.to_dict(),before)

    def test_wrong_run_frame_or_changed_inner_source_curbs_cannot_supply_a_bridge(self):
        c,a,s,mouth,inventory=prepared()
        with self.assertRaisesRegex(ValueError,'independently certified'):MergedProfileBridge().build(c,a,s,mouth)
        wrong=SimpleNamespace(frame=a.frame_chart,field=inventory.field)
        with self.assertRaisesRegex(ValueError,'consistent'):MergedProfileBridge().build(c,a,s,mouth,target_inventory=wrong)
        wrong=SimpleNamespace(frame=inventory.frame,field=inventory.field,source_sections=(999,))
        with self.assertRaisesRegex(ValueError,'source sections'):MergedProfileBridge().build(c,a,s,mouth,target_inventory=wrong)
        curves=[ray.points.copy() for ray in a.field.curves]
        curves[1]+=np.array((0.,.01))
        changed=AngularCrown(inventory.frame,curves,datum=a.field.datum,rise=a.field.rise)
        wrong=SimpleNamespace(frame=inventory.frame,field=changed)
        with self.assertRaisesRegex(ValueError,'source curb'):MergedProfileBridge().build(c,a,s,mouth,target_inventory=wrong)

    def test_actual_finite_end_and_changed_mouth_stay_rejected(self):
        c,a,s,mouth,inventory=prepared()
        altered=mouth.copy();altered[len(mouth)//2,2]+=.001
        with self.assertRaisesRegex(ValueError,'authoritative'):MergedProfileBridge().build(c,a,s,altered,target_inventory=inventory)
        station=a.length-.1;a.profile=a.cross_section(station)
        u=a.mouth_parameters(4,1.)[::-1];xy=a.point(station,a.width*(u-.5));last=np.c_[xy,a.field.values(xy)]
        with self.assertRaisesRegex(ValueError,'finite|extrapolate'):MergedProfileBridge().build(c,a,station,last,target_inventory=inventory)


if __name__=='__main__':unittest.main()
