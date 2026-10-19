"""Face-free joins count only with native evidence and both attached surfaces."""
import importlib.util
from pathlib import Path
import sys
import unittest
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location('_original_coverage_tests',
    Path(__file__).resolve().parents[1]/'way/street_network/original_span_coverage.py')
module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module
spec.loader.exec_module(module)
Owner, Coverage = module.OriginalSurfaceOwner, module.OriginalSurfaceCoverage
Region = module.SeparatorRegionOwner


class NativeRecordComponentTests(unittest.TestCase):
    @staticmethod
    def component(roles, ids=None, site=7):
        ids = ids or [31 if role==2 else 0 for role in roles]
        fields={'_cn_cache_part':('INT',roles), '_cn_cache_site':('INT',[site if role==3 else 0 for role in roles]),
                'valid':('BOOLEAN',[True]*len(roles)), 'st_cj_input_index':('INT',ids),
                'st_cj_seam':('BOOLEAN',[True]*len(roles)), 'st_cj_valid':('BOOLEAN',[True]*len(roles))}
        return SimpleNamespace(attributes={name:SimpleNamespace(domain='POINT',data_type=kind,
            data=[SimpleNamespace(value=v) for v in values]) for name,(kind,values) in fields.items()})

    def matches(self,*components):
        return module.matching_connection_records(components,7,(31,32))

    def test_legacy_native_and_mixed_storage_agree(self):
        full=self.component([1,2,2,3],[0,31,32,0])
        self.assertTrue(self.matches(full,None))
        self.assertTrue(self.matches(self.component([1]),full))
        self.assertTrue(self.matches(self.component([3]),self.component([2,2],[32,31])))

    def test_duplicate_cross_component_status_or_mouth_rejects(self):
        full=self.component([2,2,3],[31,32,0])
        self.assertFalse(self.matches(full,self.component([3])))
        self.assertFalse(self.matches(full,self.component([2],[31])))

    def test_foreign_identity_invalid_or_incomplete_records_reject(self):
        self.assertFalse(self.matches(self.component([3],site=8),self.component([2,2],[31,32])))
        self.assertFalse(self.matches(self.component([3]),self.component([2,2],[31,33])))
        self.assertFalse(self.matches(self.component([2,3],[31,0])))
        for name in ('valid','st_cj_valid','st_cj_seam'):
            row=self.component([2,2,3],[31,32,0])
            for item in row.attributes[name].data:item.value=False
            self.assertFalse(self.matches(row),name)

    def test_malformed_relevant_schema_rejects(self):
        for defect in ('missing','domain','type','length'):
            row=self.component([2,2,3],[31,32,0]);attr=row.attributes['st_cj_valid']
            if defect=='missing':del row.attributes['st_cj_valid']
            elif defect=='domain':attr.domain='FACE'
            elif defect=='type':attr.data_type='INT'
            else:attr.data.pop()
            self.assertFalse(self.matches(row),defect)

    def test_components_need_only_their_role_schema(self):
        status=self.component([3]);mouths=self.component([2,2],[31,32])
        for key in ('st_cj_input_index','st_cj_valid','st_cj_seam'):del status.attributes[key]
        for key in ('_cn_cache_site','valid'):del mouths.attributes[key]
        self.assertTrue(self.matches(status,mouths))
        self.assertFalse(module.matching_connection_records((status,mouths),7,(31,31)))


class SharedRegionCoverageTests(unittest.TestCase):
    def setUp(self):
        self.rows = [(902, 17, 900, 101), (902, 17, 901, 102),
                     (902, 17, 900, 103), (902, 17, 901, 104)]

    def test_region_face_owner_retains_both_original_bundles(self):
        shared, = Coverage.shared_separator_owners(self.rows)
        self.assertEqual(shared, Owner(7, 902, (101, 102, 103, 104), bundles=(900, 901)))
        coverage = Coverage((Owner(0, 17, shared.streets), shared))
        missing = coverage.evaluate([(0, 17)])
        self.assertEqual(missing['missing_bundles'], [900, 901])
        self.assertEqual(missing['failed_streets'], [])
        self.assertFalse(missing['complete'])
        self.assertTrue(coverage.evaluate([(0, 17), (7, 902)])['complete'])
        self.assertEqual(Coverage.shared_separator_owners(self.rows[::-1]), (shared,))

    def test_shared_region_does_not_hide_another_missing_region(self):
        shared, = Coverage.shared_separator_owners(self.rows)
        coverage = Coverage((Owner(0, 17, shared.streets), shared, Owner(4, 999, (200, 201), bundles=(999,))))
        result = coverage.evaluate([(0, 17), (7, 902)])
        self.assertFalse(result['complete']); self.assertEqual(result['missing_bundles'], [999])
        with self.assertRaises(ValueError): coverage.evaluate([(0, 17), (7, 999)])

    def test_incomplete_repeated_conflicting_and_invalid_declarations_reject(self):
        cases = [self.rows[:-1], self.rows+self.rows[:1],
                 [(902, 18, 900, 101)]+self.rows[1:],
                 [(901, p, b, s) for _, p, b, s in self.rows],
                 [(True, p, b, s) for _, p, b, s in self.rows],
                 self.rows+[(903, p, b, s+100) for _, p, b, s in self.rows]]
        for rows in cases:
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                Coverage.shared_separator_owners(rows)


class RegionCoverageTests(unittest.TestCase):
    def setUp(self):
        self.coverage=Coverage((Owner(1,0,(10,)),Owner(4,7,(10,11),bundles=(7,))),
            (Region(7,7,(10,11)),Region(7,20,(10,12))))

    def test_one_visible_region_cannot_complete_the_bundle(self):
        result=self.coverage.evaluate([(1,0),(4,7)],accepted_regions=[(7,7)])
        self.assertFalse(result['complete']);self.assertEqual(result['missing_regions'],[[7,20]])
        self.assertEqual(result['missing_bundles'],[7]);self.assertEqual(result['affected_streets'],[10,12])

    def test_all_region_evidence_is_required_and_original_ids_remain(self):
        result=self.coverage.evaluate([(1,0),(4,7)],accepted_regions=[(7,20),(7,7),(7,20)])
        self.assertTrue(result['complete']);self.assertEqual(result['accepted_regions'],[[7,7],[7,20]])
        self.assertEqual(result['missing_bundles'],[])

    def test_missing_face_attribute_unknown_region_and_wrong_bundle_reject(self):
        self.assertFalse(self.coverage.evaluate([(1,0),(4,7)])['complete'])
        for accepted,regions in (([(1,0),(4,7)],[(7,99)]), ([(1,0)],[(7,20)])):
            with self.assertRaises(ValueError):self.coverage.evaluate(accepted,accepted_regions=regions)

    def test_duplicate_and_orphan_declarations_reject(self):
        for regions in ((Region(7,20),Region(7,20)),(Region(8,20),)):
            with self.assertRaises(ValueError):Coverage((Owner(4,7),),regions)


class MatchingConnectionCoverageTests(unittest.TestCase):
    def setUp(self):
        self.owners = (Owner(0, 0, (10, 11), (100,), connection_mouths=(2, 3),
                             connection_owners=((1, 0), (3, 0))),
                       Owner(1, 0, (10,)), Owner(3, 0, (11,)))
        self.coverage = Coverage(self.owners)

    def test_native_matching_row_and_both_surfaces_complete_the_connection(self):
        result = self.coverage.evaluate(((1, 0), (3, 0)), matching_connections=((0, 0),))
        self.assertTrue(result['complete'])
        self.assertEqual(result['affected_streets'], [])

    def test_faces_alone_do_not_prove_a_face_free_connection(self):
        result = self.coverage.evaluate(((1, 0), (3, 0)))
        self.assertFalse(result['complete'])
        self.assertEqual(result['missing_junctions'], [100])
        self.assertEqual(result['affected_streets'], [10, 11])

    def test_missing_incident_surface_cannot_be_hidden_by_a_valid_row(self):
        for retained in (((1, 0),), ((3, 0),), ()):
            with self.subTest(retained=retained):
                result = self.coverage.evaluate(retained, matching_connections=((0, 0),))
                self.assertFalse(result['complete'])
                self.assertIn([0, 0], result['missing_owners'])
                self.assertEqual(result['missing_junctions'], [100])

    def test_unresolved_incidence_cannot_complete(self):
        incomplete = Owner(0, 0, (10, 11), (100,), connection_mouths=(2, 3),
                           connection_owners=((1, 0),))
        result = Coverage((incomplete, *self.owners[1:])).evaluate(
            ((1, 0), (3, 0)), matching_connections=((0, 0),))
        self.assertFalse(result['complete'])

    def test_unequal_transition_still_requires_its_own_faces(self):
        self.assertTrue(self.coverage.evaluate(((0, 0), (1, 0), (3, 0)))['complete'])
        self.assertFalse(self.coverage.evaluate(((1, 0), (3, 0)))['complete'])

    def test_unknown_and_non_connection_claims_reject(self):
        for key in ((0, 9), (1, 0), (3, 0)):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'matching connection'):
                self.coverage.evaluate(((1, 0), (3, 0)), matching_connections=(key,))

    def test_unknown_surface_and_duplicate_owner_checks_remain(self):
        with self.assertRaisesRegex(ValueError, 'Undeclared native surface'):
            self.coverage.evaluate(((1, 99),))
        with self.assertRaisesRegex(ValueError, 'Repeated original surface'):
            Coverage((*self.owners, self.owners[0]))


if __name__ == '__main__':
    unittest.main()
