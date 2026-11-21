"""Preserve original run visits, actual side-junction rims and sampling phase."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_entrances import network, updater, context_for
from tested_python_streets import model
from tested_python_streets.bundle_axis import AxialChartUnavailable
from tested_python_streets.bundle_merged_entrances import PairedAngularEntrance, preparation_key
from tested_python_streets.merged_run_inventory import (MergedRunInventory, MergedRunNeedsPairedJunction,
                                                       source_spans, inventory_frame, sample_sources)
from tested_python_streets.merged_run_inventory import evaluation_spans, hydrate_evaluation_spans
from tested_python_streets.separator_ports import PORTS, SeparatorPortIndex


def side_junction(reverse=False, segment=1):
    source=network(0)
    original=source.section_index[17]
    midpoint=tuple(.5*(np.asarray(original.points[segment])+original.points[segment+1]))
    sections=[]
    for s in source.sections:
        if s.id==17:
            sections.extend((replace(s,id=18,points=(*s.points[:segment+1],midpoint)),
                             replace(s,points=(midpoint,*s.points[segment+1:]))))
        else:sections.append(s)
    radial=np.asarray(midpoint[:2]);radial/=np.linalg.norm(radial)
    finish=np.asarray(midpoint);finish[:2]-=20.*radial
    sections.append(model.Section(401,401,(midpoint,tuple(finish)),width=4.,interpolation=0))
    sites=tuple(replace(j,ends=tuple(model.End(18,0) if e==model.End(17,0) else e for e in j.ends))
                for j in source.intersections)
    sites+= (model.Intersection(600,(model.End(18,1),model.End(17,0),model.End(401,0)),radius=2.),)
    bundle=replace(source.bundles[0],members=(model.BundleMember((model.Traversal(18),model.Traversal(17),model.Traversal(117))),
                                             source.bundles[0].members[1]))
    if reverse:
        sections=[replace(s,points=s.points[::-1],offset=-s.offset) for s in sections]
        sites=tuple(replace(j,ends=tuple(model.End(e.section_id,1-e.end) for e in j.ends)) for j in sites)
        bundle=replace(bundle,members=tuple(model.BundleMember(tuple(model.Traversal(v.section_id,False) for v in m.path))
                                           for m in bundle.members))
    return replace(source,sections=tuple(sections),intersections=sites,bundles=(bundle,))


def prepared(source):
    engine=updater(source);result=engine.prepare(source)
    context=context_for(source,result,engine.patterns)
    context.products[PORTS]=SeparatorPortIndex.build(context)
    spec=context[preparation_key(source.bundles[0])]
    approach=next(a for a in spec.layout.approaches if isinstance(a,PairedAngularEntrance) and a.end_index==1)
    return context,approach


class MergedRunInventoryTests(unittest.TestCase):
    def test_hydrated_paired_leaves_keep_complete_finite_source_reserve_and_separate_registries(self):
        from test_python_street_paired_source_body import fixture
        from tested_python_streets.paired_source_body import PairedSourceBodyPiece
        for mode, reverse, rotation, offset in ((0, False, 0., (0., 0.)), (3, True, .713, (1200., -1800.))):
            with self.subTest(reverse=reverse):
                c, bundle, original, mouths, frame, products, phase = fixture(mode, reverse, rotation, offset)
                for product in products:
                    c.products[next(site.layout_key for site in c.network.intersections
                                    if site.key == product.owner)] = product.group.physical_layout
                c.task_digests = {key: 'accepted-paired-test-product' for key in c.products}
                registry = SeparatorPortIndex.build(c)
                original_lookups, derived_lookups = [], []
                def resolve(endpoint):
                    original_lookups.append(endpoint)
                    return registry.resolve(endpoint)
                def derived(owner, sid, direction):
                    derived_lookups.append((owner, sid, direction))
                    return registry.resolve_derived(owner, sid, direction)
                c.products[PORTS] = SimpleNamespace(resolve=resolve, resolve_derived=derived)
                before = c.network.to_dict()
                spans = tuple((replace(span, curves=()),) for span in original)
                hydrated = hydrate_evaluation_spans(c, bundle, SimpleNamespace(frame_chart=frame), (),
                                                    spans=spans, frame=frame)
                leaves = tuple(role[0] for role in hydrated)
                self.assertEqual(len(original_lookups), 2)
                self.assertEqual(len(derived_lookups), 2)
                self.assertTrue(all(endpoint.section_id == 2 for endpoint in original_lookups))
                for old, leaf in zip(original, leaves):
                    self.assertEqual(old.source_intervals, leaf.source_intervals)
                    for end, mouth in zip((leaf.head, leaf.tail), mouths):
                        self.assertIs(end.prepared_mouth.approach, mouth.approach)
                        np.testing.assert_array_equal(end.prepared_mouth.xyz, mouth.xyz)
                        self.assertEqual(end.approach_index, mouth.approach_index)
                    row, support = leaf.source_intervals[0], leaf.curve_support[0]
                    self.assertLess(support['first'], row['first'])
                    self.assertGreater(support['last'], row['last'])
                    self.assertEqual(support['first'], 0.)
                    self.assertEqual(support['last'], c[model.Key('path', support['section'])].length)
                    self.assertTrue(all(not curve.flags.writeable for curve in leaf.curves))
                    self.assertFalse(old is leaf)
                    for actual, expected in zip(leaf.curves, old.curves):
                        # Source samples away from both prepared prefixes can
                        # differ, but every authoritative prefix knot survives.
                        for mouth in mouths:
                            source = next(i for i, path in enumerate(mouth.approach.paths)
                                          if path.section.id == row['section'])
                            path = mouth.approach.paths[source]
                            station = np.sort(mouth.approach.source_curve_stations[source])
                            for s in station:
                                point = path.point(s, path.section.offset +
                                    (-1. if leaf.visits[0].forward else 1.)*.5*path.section.width*
                                    (1. if actual is leaf.curves[0] else -1.))
                                self.assertLess(np.linalg.norm(actual-point, axis=1).min(), 1e-10)
                piece = PairedSourceBodyPiece()
                pavement, grass, report = piece.build(c, bundle, leaves, mouths, frame=frame,
                                                       namespace=('hydrated-test', reverse))
                self.assertLessEqual(report['stored_checks'][0]['quality']['maximum_edge_ratio'], 8.)
                for end, mouth in enumerate(mouths):
                    np.testing.assert_array_equal(pavement.xyz[piece.front_ids[end]], mouth.xyz)
                self.assertEqual(c.network.to_dict(), before)
                c.products[PORTS] = SimpleNamespace(resolve=resolve, resolve_derived=lambda *args: None)
                with self.assertRaisesRegex(AxialChartUnavailable, 'accepted owner-scoped physical port'):
                    hydrate_evaluation_spans(c, bundle, SimpleNamespace(frame_chart=frame), (),
                                            spans=spans, frame=frame)
                c.products[PORTS] = registry

    def test_derived_opposite_cuts_partition_original_stations_without_binding_original_ends(self):
        from test_python_street_one_sided_paired_junctions import fixture
        from test_python_street_compounds import footprint_context
        from tested_python_streets.bundle_visible_runs import BundleRun,MemberRun
        from tested_python_streets.one_sided_paired_junctions import OneSidedPairedSpecification
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                context,site,_,visits,opposite=fixture(reverse=reverse)
                network=context.network;sections=list(network.sections);sites=list(network.intersections)
                original=tuple(model.End(v.section_id,int(not v.forward)) for v in (visits[0],opposite))
                terminal=(model.End(visits[-1].section_id,int(visits[-1].forward)),
                          model.End(opposite.section_id,int(opposite.forward)))
                for index,endpoint in enumerate((*original,*terminal)):
                    point=network.section_index[endpoint.section_id].points[0 if endpoint.end==0 else -1]
                    ids=(500+2*index,501+2*index)
                    for sid,delta in zip(ids,((0.,20.,0.),(0.,-20.,0.))):
                        sections.append(model.Section(sid,sid,(point,tuple(np.asarray(point)+delta))))
                    sites.append(model.Intersection(200+index,(endpoint,model.End(ids[0],0),model.End(ids[1],0))))
                network=replace(network,sections=tuple(sections),intersections=tuple(sites));context=footprint_context(network)
                from tested_python_streets.intersections import IntersectionProvider
                local=IntersectionProvider.footprint(context,site)
                spec=OneSidedPairedSpecification.build(context,site,local,visits,opposite,0)
                tokens=(('bundle_end',0),('bundle_end',1))
                run=BundleRun(tuple(MemberRun(m.path,(a,b),tokens) for m,a,b in
                    zip(network.bundles[0].members,original,terminal)),tokens)
                approach=SimpleNamespace(run=run,original_ends=original,frame_chart=None,source_ranges=())
                before=network.to_dict();result=evaluation_spans(context,network.bundles[0],approach,(spec,))
                self.assertEqual(len(result[1]),2)
                ends=(result[1][0].tail,result[1][1].head)
                self.assertTrue(all(end.kind=='DERIVED_MERGED' and end.endpoint is None for end in ends))
                self.assertEqual({end.registry_key for end in ends},set(spec.opposite_cut_keys))
                owned=next(row for row in spec.source_intervals if row.section_id==opposite.section_id)
                total=sum(row['last']-row['first'] for span in result[1] for row in span.source_intervals)
                self.assertAlmostEqual(total+owned.last-owned.first,context[model.Key('path',opposite.section_id)].length)
                self.assertTrue(all(v.section_id==opposite.section_id for span in result[1] for v in span.visits))
                self.assertEqual(network.to_dict(),before)

    def test_finite_prefix_has_separately_certified_whole_run_frame_and_exact_local_chord_phase(self):
        from tested_python_streets.finite_terminal_entrances import FinitePairedEntrance
        source=side_junction();context,polar=prepared(source)
        finite=FinitePairedEntrance(context,polar.run,polar.end_index,source_tolerance=.00025)
        self.assertIs(inventory_frame(context,source.bundles[0],polar),polar.frame_chart)
        frame=inventory_frame(context,source.bundles[0],finite)
        self.assertIsNot(frame,finite.frame_chart)
        spans=source_spans(context,source.bundles[0],finite,frame=frame)
        samples=sample_sources(context,spans,finite)
        for role in spans:
            for span in role:
                for visit in span.visits:
                    path=context[model.Key('path',visit.section_id)]
                    self.assertGreater(frame.certify(path,forward=visit.forward),0.)
                    self.assertEqual(samples[visit.section_id][0],0.)
                    self.assertEqual(samples[visit.section_id][-1],path.length)
        for visit,local in zip(finite.visits,finite.source_curve_stations):
            local=np.sort(local);full=samples[visit.section_id]
            np.testing.assert_array_equal(full[(full>=local[0])&(full<=local[-1])],local)

    def test_source_junction_splits_spans_and_field_keeps_exact_local_merged_chords(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                source=side_junction(reverse);before=source.to_dict()
                context,approach=prepared(source)
                spans=source_spans(context,source.bundles[0],approach)
                self.assertEqual([[tuple(v.section_id for v in s.visits) for s in role] for role in spans],
                                 [[(18,),(17,)],[(29,)]])
                inventory=MergedRunInventory.build(context,source.bundles[0],approach)
                self.assertEqual(inventory.source_sections,(18,17,29))
                self.assertEqual(len(inventory.inner_bridges_by_role[0]),1)
                bridge=inventory.inner_bridges_by_role[0][0]
                self.assertEqual(bridge.site,model.Key('intersection',600))
                self.assertTrue(np.all(np.diff(bridge.angular_stations)>0.))
                actual=context[model.Key('junction_surface',600)].pavement.mesh.xyz
                for point in bridge.physical_xyz:
                    self.assertLess(np.linalg.norm(actual-point,axis=1).min(),1e-8)
                for role in inventory.spans_by_role:
                    for span in role:
                        self.assertEqual(len(span.source_intervals),len(span.visits))
                        for end in (span.head,span.tail):
                            if end.kind=='PHYSICAL':
                                self.assertIsNotNone(end.port)
                                np.testing.assert_array_equal(end.port.coordinates(end.port.mouth.parameters)[0],end.physical_xyz[:,:2])
                # The field shares original whole-Path chord phase at the
                # existing merged mouth and its finite source-body bridge.
                stations=approach.chart_station(np.linspace(1.,4.,7))
                for first,last in zip(inventory.field.curves,approach.field.curves):
                    np.testing.assert_array_equal(first(stations),last(stations))
                self.assertTrue(inventory.proofs['one_sided_ports_separated'])
                self.assertTrue(all(r['policy']=='COMPUTATIONAL_ONLY_UNMESHED_JUNCTION_FIELD'
                                    for r in inventory.computational_outer_bridges))
                self.assertEqual(source.to_dict(),before)

    def test_missing_source_incidence_is_never_replaced_by_coincident_control_points(self):
        source=side_junction();context,approach=prepared(source)
        changed=replace(source,intersections=tuple(j for j in source.intersections if j.id!=600))
        context.network=changed
        with self.assertRaisesRegex(ValueError,'explicit original Section incidence'):
            source_spans(context,changed.bundles[0],approach)

    def test_overlapping_one_sided_original_port_requires_paired_preparation(self):
        source=side_junction(segment=6);context,approach=prepared(source)
        from tested_python_streets.bundle_merged_entrances import preparation_key
        specification=context[preparation_key(source.bundles[0])]
        self.assertIn(model.Key('intersection',600),
            [side.owner for side in specification.paired_sources])
        # The normal provider now prepares this real overlap before inventory.
        # The old strict consumer still cannot concatenate its paired mouth.
        with self.assertRaisesRegex(AxialChartUnavailable,'intermediate paired mouth'):
            MergedRunInventory.build(context,source.bundles[0],approach)

    def test_zero_exaggeration_retains_actual_rim_xy_and_original_physical_profile_values(self):
        source=side_junction();source=replace(source,settings=replace(source.settings,exaggeration=0.))
        context,approach=prepared(source)
        inventory=MergedRunInventory.build(context,source.bundles[0],approach)
        bridge=inventory.inner_bridges_by_role[0][0]
        mesh=context[model.Key('junction_surface',600)].pavement.mesh
        for point in bridge.physical_xyz:
            self.assertLess(np.linalg.norm(mesh.xyz[:,:2]-point[:2],axis=1).min(),1e-8)
        np.testing.assert_allclose(bridge.physical_xyz[:,2],source.section_index[17].rise,atol=1e-10,rtol=0.)


if __name__=='__main__':unittest.main()
