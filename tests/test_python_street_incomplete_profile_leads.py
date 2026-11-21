"""A computed plane can advance only through a genuine incomplete-profile proof."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Network,Intersection,End,Key,Settings
from tested_python_streets.paths import Path
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.routes import Route,RouteTopology
from tested_python_streets.intersections import Approach
from tested_python_streets.updates import BuildContext
from tested_python_streets.compound_cycles import _port_footprint
from tested_python_streets.regions import Region,RegionUnion
from tested_python_streets.incomplete_profile_leads import capture,admission


def fixture():
    sections=tuple(Section(i,i,((0.,0.,0.),(*point,0.)),width=2.,interpolation=0)
                   for i,point in enumerate(((10.,0.),(-10.,0.),(0.,10.)),1))
    site=Intersection(1,tuple(End(s.id,0) for s in sections),1.)
    n=Network('missing-profile-lead',sections,(site,),settings=Settings(cell_size=1.))
    products={Key('path',s.id):Path(s) for s in sections}
    products.update((Key('profile',s.id),Crown(.08,0.)) for s in sections)
    context=BuildContext(n,products,PatternLibrary());topology=RouteTopology(n)
    arms=tuple(Approach(Route.build(context,topology.trace(End(s.id,0))),0,products[Key('profile',s.id)]) for s in sections)
    partitions={s.id:np.linspace(0.,10.,11) for s in sections}
    ports=tuple(_port_footprint(Key('cycle_port',2*a.source.id),a,at,4,1.,.00025,partitions)
                for a,at in zip(arms,(.2,4.,4.)))
    loops=RegionUnion().build(source.region for source in ports)
    record=capture(arms[0],.2,1.2,ports[0],ports,loops,1)
    group=SimpleNamespace(physical_layout=SimpleNamespace(source_cycle_profile_lead_advances=(record,)))
    return context,arms,ports,loops,record,group


class IncompleteProfileLeadTests(unittest.TestCase):
    def test_actual_partial_profile_and_existing_step_admit_real_recorded_movement(self):
        context,arms,ports,loops,record,group=fixture()
        before=context.network.to_dict()
        self.assertIsNotNone(record)
        self.assertAlmostEqual(record.exposed_measure,.5)
        self.assertTrue(record.validate(context))
        proof=admission(context,group,arms[0],.2,arms[0],1.2)
        self.assertIsNotNone(proof)
        self.assertEqual(proof['advance'],1.)
        self.assertTrue(proof['actual_old_profile_incomplete'])
        self.assertEqual(before,context.network.to_dict())

    def test_nonproof_or_changed_range_cannot_authorize_a_cut_or_override_a_root_donor(self):
        context,arms,ports,loops,record,group=fixture()
        group.physical_layout.source_cycle_profile_lead_advances=({'previous':.2,'candidate':1.2},)
        self.assertIsNone(admission(context,group,arms[0],.2,arms[0],1.2))
        group.physical_layout.source_cycle_profile_lead_advances=(record,)
        self.assertIsNone(admission(context,group,arms[0],.2,arms[0],1.3))
        group.physical_layout.source_cycle_root_support_envelopes=({'donor':1},)
        self.assertIsNone(admission(context,group,arms[0],.2,arms[0],1.2))
        for changed in (replace(record,candidate=1.3),replace(record,exposed_measure=1.),
                        replace(record,previous=-1.),replace(record,iteration=32)):
            with self.assertRaises(ValueError):changed.validate(context)

    def test_complete_profile_stale_source_geometry_and_profile_replacement_are_not_partial_proofs(self):
        context,arms,ports,loops,record,group=fixture()
        full=capture(arms[2],4.,5.,ports[2],ports,loops,1)
        self.assertIsNone(full)
        foreign=dict(context.products)
        foreign[Key('path',1)]=Path(context.network.section_index[1])
        with self.assertRaisesRegex(ValueError,'Path identity'):
            record.validate(BuildContext(context.network,foreign,context.patterns))
        foreign=dict(context.products);foreign[Key('profile',1)]=Crown(.12,0.)
        with self.assertRaisesRegex(ValueError,'complete Profile'):
            record.validate(BuildContext(context.network,foreign,context.patterns))
        changed=replace(record,exposed_measure=.25)
        with self.assertRaisesRegex(ValueError,'actual missing profile'):changed.validate(context)
        stale=replace(ports[1],region=Region(ports[1].region.owner,ports[1].region.xy+.01))
        changed=replace(record,previous_sources=(ports[0],stale,ports[2]))
        with self.assertRaisesRegex(ValueError,'original source operands'):changed.validate(context)


if __name__=='__main__':unittest.main()
