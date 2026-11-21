"""A real moved mouth retains whole and partial original owner-route prefixes."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_hydrated_external_runs import fixture
from tested_python_streets.gapped_source_inventory import source_ledger
from tested_python_streets.gapped_source_plan import _coverage
from tested_python_streets.model import Key


def inventory(spans):
    return SimpleNamespace(run=SimpleNamespace(members=tuple(
        SimpleNamespace(visits=tuple(span.original_visits or span.visits)) for span in spans)),
        spans_by_role=tuple((span,) for span in spans),owned_spans=())


class GappedOwnedPrefixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixtures=[fixture(short_terminal=True),
                      fixture(True,.713,(1200.,-1800.),short_terminal=True)]

    def test_real_routed_terminal_preserves_whole_and_partial_prefixes_in_original_order(self):
        for context,bundle,spans,*_ in self.fixtures:
            before=context.network.to_dict();selected=inventory(spans)
            original=[tuple(v.section_id for v in member.visits) for member in selected.run.members]
            self.assertEqual(original,[(1,2,3),(4,5,6)])
            self.assertEqual(tuple(v.section_id for v in spans[0].visits),(1,2))
            positions=[end.physical_xyz.copy() for span in spans for end in (span.head,span.tail)]
            prefixes,proof=source_ledger(context,selected.run,selected.spans_by_role,())
            self.assertEqual([[row['section'] for row in own] for own in proof['original_ordered_visits']],
                             [[1,2,3],[4,5,6]])
            whole=next(row for row in prefixes if row['section']==3)
            partial=next(row for row in prefixes if row['section']==2)
            self.assertEqual((whole['first'],whole['last']),(0.,context[Key('path',3)].length))
            self.assertLess(partial['last']-partial['first'],context[Key('path',2)].length)
            for row in prefixes:
                self.assertIs(row['source_path'],context[Key('path',row['section'])])
                self.assertEqual(row['owner'],spans[0].tail.owner)
                self.assertLessEqual(row['containment']['maximum_error'],.002)
            for end,saved in zip((end for span in spans for end in (span.head,span.tail)),positions):
                np.testing.assert_array_equal(end.physical_xyz,saved)
            self.assertEqual(context.network.to_dict(),before)

    def test_complete_extent_includes_query_owned_section_without_a_body_visit(self):
        for context,_,spans,*_ in self.fixtures:
            from tested_python_streets.intersections import IntersectionProvider
            context=replace(context,products=dict(context.products))
            for connection in context.network.connections:
                context.products[connection.layout_key]=IntersectionProvider.layout(context,connection)
            # Select the actual first-role source interval; its remote ordinary
            # owner consumes all Section 3 and a finite part of Section 2.
            selected=inventory(spans[:1])
            rows=_coverage(context,selected)
            self.assertEqual([row['section'] for row in rows],[1,2,3])
            whole=rows[-1]['body_or_whole_owner_intervals']
            self.assertEqual(len(whole),1)
            self.assertTrue(whole[0]['complete_original_route_prefix'])
            self.assertEqual((whole[0]['first'],whole[0]['last']),(0.,context[Key('path',3)].length))
            self.assertTrue(all(row['complete_original_extent'] for row in rows))

    def test_duplicate_missing_reordered_and_wrong_owner_prefixes_reject(self):
        context,_,spans,*_=self.fixtures[0];selected=inventory(spans)
        span=spans[0];claims=span.owned_source_intervals
        whole=next(row for row in claims if row['section']==3)
        partial=next(row for row in claims if row['section']==2)
        variants=(
            replace(span,owned_source_intervals=(*claims,whole)),
            replace(span,owned_source_intervals=tuple(row for row in claims if row is not whole)),
            replace(span,owned_source_intervals=tuple(dict(row,owner=Key('intersection',202)) for row in claims)),
            replace(span,owned_source_intervals=tuple(dict(row,source_path=context[Key('path',1)]) for row in claims)),
            replace(span,owned_source_intervals=tuple(dict(row,first=row['first']+.001) if row is partial else row for row in claims)),
            replace(span,owned_source_intervals=tuple(dict(row,forward=not row['forward']) for row in claims)),
            replace(span,owned_source_intervals=tuple(dict(row,section=999999) for row in claims)),
            replace(span,original_visits=span.original_visits[::-1]),
            replace(span,tail=replace(span.tail,owner=Key('intersection',202)),
                    owned_source_intervals=tuple(dict(row,owner=Key('intersection',202)) for row in claims)),
            replace(span,owned_source_intervals=tuple(dict(row,containment=dict(row['containment'],maximum_error=.003)) for row in claims)),
        )
        for index,wrong in enumerate(variants):
            with self.subTest(index=index),self.assertRaises(ValueError):
                source_ledger(context,selected.run,((wrong,),selected.spans_by_role[1]),())
        with self.assertRaisesRegex(ValueError,'reorders or repeats'):
            source_ledger(context,selected.run,((span,span),selected.spans_by_role[1]),())


if __name__=='__main__':unittest.main()
