"""A source contact retry cannot spend source width or profile geometry."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_source_cycles import prepare,triangle
from tested_python_streets.compounds import rim_label
from tested_python_streets.finite_free_caps import PreparedCapCurb
from tested_python_streets.model import Key
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.sampling import AnalyticalBoundarySampler,AnalyticalContactError
from tested_python_streets.source_contact_recovery import recover
from tested_python_streets.source_rim_normalization import normalize


def profiles(group):
    inventory={s.region.owner:s for s in group.footprints}
    return np.asarray([loop.xy[i] for loop in group.loops for i,span in enumerate(loop.spans)
        if rim_label(span,inventory).kind=='mouth'])


class SourceContactRecoveryTests(unittest.TestCase):
    def test_full_source_cycle_keeps_profiles_and_rejects_repeated_or_clipped_retry(self):
        context,group=prepare(triangle(8.))
        original=context.network.to_dict();before=profiles(group)
        parent=CompoundQuadChart(context.patterns,spacing=context.network.settings.spacing)
        chart=recover(parent,group,context,AnalyticalContactError('previous local contact exhausted'))
        np.testing.assert_array_equal(profiles(chart.partition.group),before)
        self.assertEqual(context.network.to_dict(),original)
        report=chart.metadata['source_contact_normalization']
        self.assertTrue(report['contact_tolerance_unchanged'])
        self.assertEqual(report['maximum_added_polyline_projection'],1e-6)
        self.assertTrue(report['original_width_proof'] and report['final_width_proof'])
        with self.assertRaisesRegex(ValueError,'unnormalized'):
            recover(parent,normalize(group),context,AnalyticalContactError('again'))
        loop=group.loops[0];center=loop.xy.mean(axis=0)
        clipped=replace(group,loops=(replace(loop,xy=center+.1*(loop.xy-center)),))
        with self.assertRaises(ValueError):recover(parent,clipped,context,AnalyticalContactError('clipped'))

    def test_only_exhausted_contact_branches_receive_the_recovery_error_type(self):
        sampler=AnalyticalBoundarySampler(tolerance=.005,spacing=.7)
        owners=(Key('intersection',1),Key('intersection',2))
        spans=tuple(SimpleNamespace(owner=k) for k in owners)
        labels=(SimpleNamespace(index=0,last=.5),SimpleNamespace(index=0,first=.5))
        def layouts(a,b):return {k:SimpleNamespace(boundaries=(curve,),rows=8) for k,curve in zip(owners,(a,b))}
        horizontal=PreparedCapCurb(np.array(((0.,0.),(1.,0.))))
        parallel=PreparedCapCurb(np.array(((0.,1.),(1.,1.))))
        with self.assertRaises(AnalyticalContactError):sampler._contact(*spans,labels,layouts(horizontal,parallel))
        vertical=PreparedCapCurb(np.array(((.5,-1.),(.5,1.))))
        horizontal.contact_reference_points=horizontal.points+100.
        with self.assertRaises(ValueError) as failure:
            sampler._contact(*spans,labels,layouts(horizontal,vertical))
        self.assertIs(type(failure.exception),ValueError)
        self.assertIn('source approximation budget',str(failure.exception))


if __name__=='__main__':unittest.main()
