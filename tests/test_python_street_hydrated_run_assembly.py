"""Source-owned leaves assemble by physical identity with complete inventories."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_hydrated_external_runs import fixture as external_fixture
from test_python_street_paired_source_body import fixture as paired_fixture
from tested_python_streets.hydrated_run_assembly import HydratedRunAssembly, interval_inventory
from tested_python_streets.model import Key
from tested_python_streets.ownership import surface_key


def owned_scope(c,spans,mouth):
    opposite=next(i for i,end in enumerate((spans[0].head,spans[1].head)) if end.registry_key is not None)
    product=c[surface_key(mouth.owner)]
    first_sid=spans[opposite].visits[0].section_id
    cuts=[alias.station for a,s in zip(product.group.physical_layout.approaches,
                                      product.group.physical_layout.outward_stations)
          if hasattr(a,'source_mouths') for alias in a.source_mouths(s)
          if not alias.original_incidence and alias.endpoint.section_id==first_sid]
    result=[]
    for role,span in enumerate(spans):
        rows=[]
        for visit in span.visits:
            path=c[Key('path',visit.section_id)]
            low=min(cuts) if role==opposite and visit is span.visits[0] else 0.
            rows.append((visit.section_id,visit.forward,low,path.length))
        result.append(tuple(rows))
    return tuple(result)


class HydratedRunAssemblyTests(unittest.TestCase):
    def test_consumed_terminal_section_keeps_original_identity_and_owner_route_scope(self):
        c,bundle,spans,mouth,frame,_,_,_=external_fixture(short_terminal=True)
        self.assertEqual([v.section_id for v in spans[0].visits],[1,2])
        self.assertEqual([v.section_id for v in spans[0].original_visits],[1,2,3])
        self.assertEqual(spans[0].tail.original_terminal.section_id,3)
        self.assertEqual(spans[0].tail.endpoint.section_id,2)
        whole=next(claim for claim in spans[0].owned_source_intervals if claim['section']==3)
        self.assertEqual((whole['first'],whole['last']),(0.,c[Key('path',3)].length))
        self.assertEqual(whole['owner'],spans[0].tail.owner)
        self.assertLessEqual(whole['containment']['maximum_error'],.002)
        leaves=tuple((span,) for span in spans)
        scope=list(owned_scope(c,spans,mouth))
        scope[0]=(*scope[0],(3,True,0.,c[Key('path',3)].length))
        assembly=HydratedRunAssembly()
        _,_,report=assembly.build(c,bundle,leaves,frame=frame,namespace='moved-terminal',
            expected_intervals=interval_inventory(leaves),scope_intervals=tuple(scope))
        self.assertTrue(report['coverage']['complete_selected_scope'])
        self.assertTrue(any(row['section']==3 and row['proof']=='complete_original_approach_route_prefix'
                            for row in report['junction_owned_intervals']))

    def test_actual_external_piece_has_complete_selected_scope_and_query_only_owned_junction(self):
        c,bundle,spans,mouth,frame,_,_,_=external_fixture()
        leaves=tuple((span,) for span in spans);expected=interval_inventory(leaves)
        scope=owned_scope(c,spans,mouth);before=c.network.to_dict()
        assembly=HydratedRunAssembly()
        pavement,grass,report=assembly.build(c,bundle,leaves,frame=frame,namespace='external-assembly',
                                            expected_intervals=expected,scope_intervals=scope)
        self.assertTrue(report['coverage']['complete_selected_scope'])
        self.assertTrue(report['coverage']['complete_supplied_inventory'])
        self.assertEqual(report['original_interval_inventory'],expected)
        self.assertEqual(report['query_only_junction_owners'],[mouth.owner.text()])
        self.assertFalse(report['junction_geometry_republished'])
        self.assertFalse(report['runtime_discovery_enabled'])
        self.assertEqual(len(assembly.body_query_mesh.quads),len(pavement.quads)+len(grass.quads))
        product=c[surface_key(mouth.owner)]
        cells=product.pavement.mesh.xyz[product.pavement.mesh.quads]
        selected=int(np.argmax(cells[:,:,2].mean(axis=1)))
        point=cells[selected,:,:2].mean(axis=0)
        expected_height=cells[selected,:,2].mean()
        self.assertGreater(expected_height,.001)
        self.assertAlmostEqual(float(assembly.height(point)),float(expected_height),places=9)
        self.assertGreater(len(assembly.query_mesh.quads),len(assembly.body_query_mesh.quads))
        self.assertEqual(c.network.to_dict(),before)
        self.assertEqual(len(assembly.pavement_keys),len(pavement.xyz))
        self.assertEqual(len(set(assembly.pavement_keys)),len(pavement.xyz))

    def test_pair_between_two_products_uses_inventory_only_mode_without_republishing_owners(self):
        c,bundle,spans,mouths,frame,products,_=paired_fixture()
        leaves=tuple((span,) for span in spans)
        assembly=HydratedRunAssembly()
        pavement,grass,report=assembly.build(c,bundle,leaves,frame=frame,namespace='paired-assembly',
                                            expected_intervals=interval_inventory(leaves))
        self.assertFalse(report['coverage']['complete_selected_scope'])
        self.assertEqual(set(assembly.paired_owners),{mouth.owner for mouth in mouths})
        count=len(pavement.quads)+(len(grass.quads) if grass is not None else 0)
        self.assertEqual(len(assembly.body_query_mesh.quads),count)
        actual=count+sum(len(product.pavement.mesh.quads)+
                        (len(product.grass.mesh.quads) if product.grass is not None else 0) for product in products)
        self.assertEqual(len(assembly.query_mesh.quads),actual)
        self.assertEqual(report['total_quads'],count)

    def test_missing_duplicate_unmatched_leaves_and_unowned_scope_are_rejected(self):
        c,bundle,spans,mouth,frame,_,_,_=external_fixture()
        leaves=tuple((span,) for span in spans);expected=interval_inventory(leaves)
        for incomplete in (((),()),(leaves[0],())):
            with self.assertRaisesRegex(ValueError,'complete supplied'):
                HydratedRunAssembly().build(c,bundle,incomplete,frame=frame,namespace='missing',expected_intervals=expected)
        duplicate=tuple((span,span) for span in spans)
        with self.assertRaisesRegex(ValueError,'repeats or reorders'):
            HydratedRunAssembly().build(c,bundle,duplicate,frame=frame,namespace='duplicate',expected_intervals=expected)
        damaged=replace(spans[1],head=replace(spans[1].head,approach_index=999))
        with self.assertRaisesRegex(ValueError,'boundary token'):
            HydratedRunAssembly().build(c,bundle,(leaves[0],(damaged,)),frame=frame,
                namespace='unmatched',expected_intervals=expected)
        scope=tuple(tuple((v.section_id,v.forward,0.,c[Key('path',v.section_id)].length)
                          for v in span.visits) for span in spans)
        with self.assertRaisesRegex(ValueError,'partition|unowned'):
            HydratedRunAssembly().build(c,bundle,leaves,frame=frame,namespace='unowned',
                                        expected_intervals=expected,scope_intervals=scope)


if __name__=='__main__':unittest.main()
