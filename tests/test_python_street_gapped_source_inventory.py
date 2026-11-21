"""Separated side nodes retain real profiles, actual curbs and source topology."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_run_inventory import side_junction,prepared
from tested_python_streets.gapped_source_inventory import GappedSourceInventory,GappedSourceUnavailable
from tested_python_streets.merged_run_inventory import physical_source_spans,hydrate_physical_spans,inventory_frame
from tested_python_streets.model import Key
from tested_python_streets.separator_ports import PORTS
from tested_python_streets.bundle_axis import AxialChartUnavailable


class GappedSourceInventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixtures=[]
        for reverse in (False,True):
            source=side_junction(reverse)
            context,entrance=prepared(source)
            cls.fixtures.append((context,entrance))

    def test_original_physical_partition_hydrates_only_actual_registered_mouths(self):
        for context,entrance in self.fixtures:
            with self.subTest(reversed=not entrance.run.members[0].visits[0].forward):
                bundle=context.network.bundles[0];frame=inventory_frame(context,bundle,entrance)
                before=context.network.to_dict();registry=context[PORTS]
                ports=dict(registry.ports);derived=dict(registry.derived_ports)
                spans=physical_source_spans(context,bundle,entrance.run,frame=frame)
                self.assertEqual([len(own) for own in spans],[2,1])
                self.assertTrue(all(end.kind=='PHYSICAL' and end.endpoint is not None
                    and end.prepared_mouth is None for own in spans for span in own for end in (span.head,span.tail)))
                hydrated=hydrate_physical_spans(context,bundle,spans,frame=frame)
                for role in (0,1):
                    tail=hydrated[role][-1].tail
                    self.assertEqual(tail.kind,'MERGED')
                    self.assertIs(tail.prepared_mouth.approach,entrance)
                    self.assertIsNone(tail.registry_key)
                self.assertEqual(registry.ports,ports);self.assertEqual(registry.derived_ports,derived)
                self.assertEqual(context.network.to_dict(),before)

    def test_actual_separated_gap_keeps_published_inner_curb_and_finite_original_paths(self):
        for context,entrance in self.fixtures:
            with self.subTest(reversed=not entrance.run.members[0].visits[0].forward):
                bundle=context.network.bundles[0];before=context.network.to_dict()
                inventory=GappedSourceInventory.build(context,bundle,run=entrance.run,
                    frame=inventory_frame(context,bundle,entrance))
                self.assertEqual(inventory.source_sections,(18,17,29))
                self.assertEqual(len(inventory.gaps),1);self.assertFalse(inventory.owned_spans)
                gap=inventory.gaps[0];self.assertEqual(gap.bridge.owner,Key('intersection',600))
                self.assertGreater(gap.minimum_gap,.001)
                actual=context[Key('junction_surface',600)].pavement.mesh.xyz
                for point in gap.bridge.physical_xyz:
                    self.assertLess(np.linalg.norm(actual-point,axis=1).min(),1e-8)
                for own in inventory.spans_by_role:
                    for span in own:
                        for piece in span.route.pieces:
                            self.assertIs(piece.path,context[Key('path',piece.section.id)])
                self.assertTrue(all(row['policy']=='COMPUTATIONAL_ONLY_UNMESHED_JUNCTION_FIELD'
                                    for row in inventory.computational_outer_bridges))
                self.assertEqual(inventory.proofs['new_source_aliases'],0)
                stations=entrance.chart_station(np.linspace(1.,4.,7))
                for first,last in zip(inventory.field.curves,entrance.field.curves):
                    np.testing.assert_array_equal(first(stations),last(stations))
                self.assertEqual(context.network.to_dict(),before)

    def test_absent_authored_incidence_and_finite_allocations_cannot_create_a_gap(self):
        context,entrance=self.fixtures[0];network=context.network;bundle=network.bundles[0]
        frame=inventory_frame(context,bundle,entrance)
        context.network=replace(network,intersections=tuple(site for site in network.intersections if site.id!=600))
        try:
            with self.assertRaisesRegex(ValueError,'original Section incidence'):
                physical_source_spans(context,bundle,entrance.run,frame=frame)
        finally:context.network=network
        with self.assertRaises(GappedSourceUnavailable):
            GappedSourceInventory.build(context,bundle,run=entrance.run,frame=frame,maximum_gaps=0)
        with self.assertRaisesRegex(ValueError,'budget'):
            GappedSourceInventory.build(context,bundle,run=entrance.run,frame=frame,max_samples=4)


if __name__=='__main__':unittest.main()
