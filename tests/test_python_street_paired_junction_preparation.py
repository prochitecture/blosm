"""Original incidence discovery and fresh canonical paired preparation."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_paired_source_sampling import two_sites
from test_python_street_compounds import footprint_context
from test_python_street_merged_run_inventory import side_junction, prepared
from tested_python_streets import model
from tested_python_streets.bundle_merged_entrances import preparation_key
from tested_python_streets.bundle_visible_runs import BundleRun, MemberRun
from tested_python_streets.intersections import IntersectionProvider
from tested_python_streets.ownership import footprint_key
from tested_python_streets.paired_junction_preparation import (discover, prepare,
    PairedJunctionPreparationUnavailable)
from tested_python_streets.paired_source_sampling import verify_shared_source_phase


def discovery_context(reverse=False):
    context, _ = two_sites()
    source = context.network
    if reverse:
        angle=.713
        rotation=np.array(((np.cos(angle),-np.sin(angle)),(np.sin(angle),np.cos(angle))))
        sections=tuple(replace(s, points=tuple(tuple((*xy,z)) for xy,z in zip(
            np.asarray(s.points)[::-1,:2]@rotation.T+(1200.,-1800.),np.asarray(s.points)[::-1,2])),
            offset=-s.offset) for s in source.sections)
        sites=tuple(replace(j, ends=tuple(model.End(e.section_id,1-e.end) for e in j.ends))
                    for j in source.intersections)
        bundle=replace(source.bundles[0],members=tuple(model.BundleMember(tuple(
            model.Traversal(v.section_id,not v.forward) for v in m.path)) for m in source.bundles[0].members))
        source=replace(source,sections=sections,intersections=sites,bundles=(bundle,))
        context=footprint_context(source)
    bundle=source.bundles[0]
    members=tuple(MemberRun(tuple(m.path),(),()) for m in bundle.members)
    # This fixture uses the real original member inventory; no physical End or
    # synthetic persistent junction is introduced by its run descriptor.
    approach=SimpleNamespace(run=BundleRun(members,()),source_ranges=())
    terminal=SimpleNamespace(bundle_id=bundle.id,layout=SimpleNamespace(approaches=(approach,)))
    for site in source.intersections:
        context.products[footprint_key(site.key)]=IntersectionProvider.footprint(context,site)
    return context,bundle,terminal


class PairedJunctionPreparationTests(unittest.TestCase):
    def test_discovery_uses_all_original_incidence_and_explicit_opposite_visits(self):
        for reverse in (False,True):
            context,bundle,terminal=discovery_context(reverse)
            before=context.network.to_dict()
            requests=discover(context,bundle,terminal)
            self.assertEqual([r.site.id for r in requests],[100,101])
            self.assertEqual([tuple(v.section_id for v in r.near_visits) for r in requests],[(1,2),(2,5)])
            self.assertTrue(all(r.opposite_visit.section_id==3 for r in requests))
            self.assertTrue(all(r.opposite_visit.forward != reverse for r in requests))
            self.assertEqual(context.network.to_dict(),before)

    def test_missing_incidence_and_overlapping_opposite_choices_reject(self):
        context,bundle,terminal=discovery_context()
        source=context.network
        context.network=replace(source,intersections=source.intersections[1:])
        with self.assertRaisesRegex(PairedJunctionPreparationUnavailable,'authored incidence'):
            discover(context,bundle,terminal)
        context,bundle,terminal=discovery_context()
        other=replace(context.network.section_index[3],id=7,street_id=7)
        bundle=replace(bundle,members=(bundle.members[0],model.BundleMember((model.Traversal(3),model.Traversal(7)))))
        network=replace(context.network,sections=(*context.network.sections,other),bundles=(bundle,))
        changed=footprint_context(network)
        for site in network.intersections:
            changed.products[footprint_key(site.key)]=IntersectionProvider.footprint(changed,site)
        members=tuple(MemberRun(tuple(m.path),(),()) for m in bundle.members)
        terminal=SimpleNamespace(bundle_id=bundle.id,layout=SimpleNamespace(approaches=(
            SimpleNamespace(run=BundleRun(members,()),source_ranges=()),)))
        with self.assertRaisesRegex(PairedJunctionPreparationUnavailable,'explicit opposite source interval'):
            discover(changed,bundle,terminal)

    def test_fresh_terminal_and_sides_share_exact_phase_without_mutating_accepted_mouths(self):
        source=side_junction(segment=6)
        context,_=prepared(source)
        bundle=source.bundles[0]
        terminal=context[preparation_key(bundle)]
        old=[(a,tuple(x.copy() for x in a.source_curve_stations))
             for a in terminal.layout.approaches if hasattr(a,'source_curve_stations')]
        result,sides,parameters,proof=prepare(context,bundle,terminal)
        self.assertEqual([s.owner.id for s in sides],[600])
        self.assertIsNot(result,terminal)
        self.assertEqual(result.paired_sources,sides)
        self.assertIs(result.source_phase_proof,proof)
        self.assertTrue(proof['original_path_identity'])
        self.assertTrue(proof['overlapping_station_arrays_exact'])
        for a,values in old:
            for first,last in zip(a.source_curve_stations,values):
                np.testing.assert_array_equal(first,last)
        entrances=[a for a in result.layout.approaches if hasattr(a,'source_curve_stations')]
        entrances.extend(a for spec in sides for a in spec.entrances)
        self.assertEqual(verify_shared_source_phase(entrances),proof)
        for a in entrances:
            for path in a.paths:
                self.assertIs(path,context[model.Key('path',path.section.id)])
        self.assertEqual(source.to_dict(),context.network.to_dict())
        self.assertTrue(all(isinstance(values,tuple) and len(values)<=20000 for values in parameters.values()))
        again,again_sides,again_parameters,again_proof=prepare(context,bundle,result)
        self.assertEqual(again_parameters,parameters)
        self.assertEqual(again_proof,proof)
        requests=discover(context,bundle,result)
        wrong=replace(requests[0],opposite_visit=requests[0].near_visits[0])
        with self.assertRaisesRegex(PairedJunctionPreparationUnavailable,'original source incidence'):
            prepare(context,bundle,result,requests=(wrong,))


if __name__=='__main__':unittest.main()
