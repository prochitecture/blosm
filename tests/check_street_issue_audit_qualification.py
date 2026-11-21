"""Native three-T acceptance tests, including failure after READY publication.

Run in an isolated factory Blender with -- --artifacts DIR. This tests the
audit gate; its controlled faults do not alter street geometry algorithms.
"""
import argparse
from contextlib import nullcontext
import json
from pathlib import Path
import sys
from unittest.mock import patch
from xml.etree import ElementTree as ET

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_street_issue_audit as audit
from street_issue_qualification import ground_qualification_passed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts', required=True, type=Path)
    parser.add_argument('--accepted-audit', action='append', default=[], type=Path)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.artifacts.mkdir(parents=True, exist_ok=True)
    fixture = args.artifacts / 'three_t.osm'
    audit.support.synthetic_osm(fixture)
    excluded = args.artifacts / 'three_t_elevated_exclusion.osm'
    tree = ET.parse(fixture)
    for way in tree.getroot().findall('way'):
        if int(way.get('id')) >= 120:
            ET.SubElement(way, 'tag', k='bridge', v='yes')
            ET.SubElement(way, 'tag', k='layer', v='1')
    tree.write(excluded, encoding='utf-8', xml_declaration=True)

    enable = audit.support.enable
    package = None
    def enable_once():
        nonlocal package
        if package is None:
            package = enable()
        return package
    audit.support.enable = enable_once
    records = []

    def run(name, paths, expected, fault=None):
        target = args.artifacts / name
        sys.argv = [__file__, '--', '--artifacts', str(target), '--full-mesh', '--require-ground-coverage']
        for path in paths:
            sys.argv.extend(('--osm', str(path)))
        error = None
        with fault if fault is not None else nullcontext():
            try:
                audit.main()
            except AssertionError as exc:
                error = str(exc)
        rows = json.loads((target / 'preparation_audit.json').read_text())
        assert bool(error is None) == expected, (name, error, rows)
        assert all(ground_qualification_passed(row) == row['ground_qualification_passed'] for row in rows)
        records.append(dict(case=name, passed=True, audit_accepted=error is None,
            import_outcomes=[row['import_result'].get('outcome') for row in rows],
            ground_accepted=[row['ground_qualification_passed'] for row in rows],
            report=str(target / 'preparation_audit.json')))
        return rows

    run('normal', (fixture,), True)
    from blosm.way.street_network.python_backend.renderer import PythonStreetRenderer
    def fail_finalization(self):
        raise RuntimeError('CONTROLLED acceptance regression: late renderer finalization failure')
    late = run('late_finalize_failure', (fixture,), False,
        patch.object(PythonStreetRenderer, 'finalize', fail_finalization))
    assert late[0]['import_result']['outcome'] == 'partial'
    assert late[0]['import_reports'][0]['ground_coverage'] == 'COMPLETE'
    assert all(s['status'] == 'READY' for s in late[0]['import_reports'][0]['surfaces'])

    def cancel_finalization(self):
        self.app.osm_session.cancel()
    cancelled = run('cancelled_after_ready', (fixture,), False,
        patch.object(PythonStreetRenderer, 'finalize', cancel_finalization))
    assert cancelled[0]['import_result']['outcome'] == 'cancelled'

    original_prepare = audit.support.StreetProgressiveChecks.prepare
    calls = 0
    def fail_second_prepare(self, *positional, **keywords):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise RuntimeError('CONTROLLED acceptance regression: failure before second-input cleanup')
        return original_prepare(self, *positional, **keywords)
    stale = run('stale_prior_reports', (fixture, fixture), False,
        patch.object(audit.support.StreetProgressiveChecks, 'prepare', fail_second_prepare))
    assert stale[1].get('error') and not stale[1].get('operator') and not stale[1]['import_result']

    original_report = PythonStreetRenderer._write_report
    def remove_final_report(self, *positional, **keywords):
        context = original_report(self, *positional, **keywords)
        if keywords['complete']:
            del context.collection['blosm_street_report']
        return context
    missing = run('missing_report', (fixture,), False,
        patch.object(PythonStreetRenderer, '_write_report', remove_final_report))
    assert missing[0]['import_result']['outcome'] == 'success' and not missing[0]['import_reports']

    elevated = run('explicit_elevated_exclusion', (excluded,), True)
    reports = elevated[0]['import_reports']
    assert any(report['ground_scope']['deferred_sections'] for report in reports)
    assert all(report['ground_coverage'] == 'COMPLETE' for report in reports)

    revalidated = []
    for path in args.accepted_audit:
        rows = json.loads(path.read_text())
        assert rows and all(ground_qualification_passed(row) for row in rows), path
        revalidated.append(dict(path=str(path), accepted=True, inputs=len(rows),
            scope='existing actual audit record rechecked; no new native mesh run'))
    report = dict(passed=True, cases=records, revalidated_existing_audits=revalidated)
    (args.artifacts / 'qualification_checks.json').write_text(json.dumps(report, indent=2))
    print('STREET_AUDIT_QUALIFICATION_CHECKS', json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
