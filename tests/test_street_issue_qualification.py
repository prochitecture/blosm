"""Import outcome and current READY reports must qualify together."""
from copy import deepcopy
import unittest
from street_issue_qualification import ground_qualification_passed


def accepted():
    return dict(operator=['FINISHED'], import_result=dict(outcome='success', errors=[]),
        import_reports=[dict(ground_coverage='COMPLETE', surfaces=[dict(status='READY')])])


class StreetIssueQualificationTests(unittest.TestCase):
    def test_complete_import_and_explicit_elevated_exclusion_are_accepted(self):
        row = accepted()
        self.assertTrue(ground_qualification_passed(row))
        report = row['import_reports'][0]
        report.update(coverage='PARTIAL', preparation_issues=[dict(code='OUTSIDE_GROUND_SCOPE')],
            ground_scope=dict(deferred_sections=[7]))
        self.assertTrue(ground_qualification_passed(row))

    def test_ready_reports_reject_late_failure_and_stale_prior_success(self):
        row = accepted()
        row['import_result'].update(outcome='partial', errors=['failure after COMPLETE report'])
        self.assertFalse(ground_qualification_passed(row))
        row = accepted()
        row['error'] = 'failure before second-input cleanup'
        row.pop('operator')
        self.assertFalse(ground_qualification_passed(row))

    def test_missing_cancelled_failed_or_inconsistent_imports_reject(self):
        for outcome in (None, 'partial', 'failed', 'cancelled'):
            row = accepted(); row['import_result']['outcome'] = outcome
            with self.subTest(outcome=outcome):
                self.assertFalse(ground_qualification_passed(row))
        for operator in (None, [], ['CANCELLED'], ['FINISHED', 'CANCELLED']):
            row = accepted(); row['operator'] = operator
            with self.subTest(operator=operator):
                self.assertFalse(ground_qualification_passed(row))
        for result in (None, {}, dict(outcome='success'), dict(outcome='success', errors=['cleanup failed'])):
            row = accepted(); row['import_result'] = result
            self.assertFalse(ground_qualification_passed(row))

    def test_zero_missing_or_unfinished_ground_reports_reject(self):
        for reports in (None, [], [{}], [dict(ground_coverage='COMPLETE')],
                [dict(ground_coverage='COMPLETE', surfaces=[])],
                [dict(ground_coverage='IN_PROGRESS', surfaces=[dict(status='READY')])],
                [dict(ground_coverage='COMPLETE', surfaces=[dict(status='FAILED')])],
                [dict(ground_coverage='COMPLETE', surfaces=[{}])]):
            row = accepted(); row['import_reports'] = reports
            with self.subTest(reports=reports):
                self.assertFalse(ground_qualification_passed(row))
        row = accepted(); row['import_reports'].append(deepcopy(row['import_reports'][0]))
        row['import_reports'][1]['ground_coverage'] = 'PARTIAL'
        self.assertFalse(ground_qualification_passed(row))


if __name__ == '__main__':
    unittest.main()
