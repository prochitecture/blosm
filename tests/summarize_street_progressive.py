"""Summarize GUI import evidence, including partial coverage, without Blender."""
import argparse
import json
from collections import Counter
from pathlib import Path

def percentile(values, fraction):
    values = sorted(values)
    if not values:
        return None
    index = (len(values) - 1) * fraction
    lo = int(index)
    hi = min(lo + 1, len(values) - 1)
    return values[lo] * (1 - (index - lo)) + values[hi] * (index - lo)


def summarize(folder):
    report = json.loads((folder / 'ui_checks.json').read_text(encoding='utf-8'))
    resources = json.loads((folder / 'process/resources.json').read_text(encoding='utf-8'))
    beats = [row for row in report['heartbeats'] if row['phase'] == 'normal completion']
    gaps = [dict(seconds=b['seconds'] - a['seconds'], before=a['session_phase'],
                 after=b['session_phase']) for a, b in zip(beats, beats[1:])]
    periods = [row['seconds'] for row in gaps]
    surfaces = [s for r in report['street_import_reports'] for s in r['surfaces']]
    issues = [issue for r in report['street_import_reports'] for issue in r['preparation_issues']]
    prepared = sum(len(s['sections']) for s in surfaces)
    accepted = sum(len(s['sections']) for s in surfaces if s['status'] == 'READY')
    preview_beats = [b for b in beats if b['visible_previews']]
    visible_beats = [b for b in beats if b['visible_components']]
    return dict(run=folder.name, source=report['source_path'],
        source_sha256=report['source_sha256'], qualification_passed=report['passed'],
        result=report['normal_result'], counts=report['scene_counts'],
        tested_prepared_sections=prepared, accepted_sections=accepted,
        accepted_section_percent=(100 * accepted / prepared if prepared else None),
        heartbeat_interval=report['observer_interval_seconds'],
        first_visible_preview_seconds=(preview_beats[0]['import_seconds'] if preview_beats else None),
        first_visible_accepted_seconds=(visible_beats[0]['import_seconds'] if visible_beats else None),
        observed_preview_counts=sorted(set(b['visible_previews'] for b in preview_beats)),
        heartbeat_gap_p50=percentile(periods, .5), heartbeat_gap_p95=percentile(periods, .95),
        heartbeat_gap_p99=percentile(periods, .99), heartbeat_gap_max=max(periods, default=None),
        gaps_over_1s=sum(x > 1 for x in periods), gaps_over_5s=sum(x > 5 for x in periods),
        longest_heartbeat_gaps=sorted(gaps, key=lambda x:x['seconds'], reverse=True)[:8],
        longest_callbacks=report['wall_timing']['scopes']['normal_import']['controller_tick']['longest_callbacks'],
        failed_components=[dict(sections=len(s['sections']), streets=len(s['streets']),
                                reason=s.get('reason')) for s in surfaces if s['status'] != 'READY'],
        preparation_issue_counts=dict(Counter(i['reason'] for i in issues)),
        peak_resident_gib=resources['peak_resident']/2**30,
        sampled_peak_private_gib=resources['peak_private']/2**30,
        resource_status=resources['status'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifacts', type=Path, required=True)
    args = parser.parse_args()
    summaries = [summarize(p.parent) for p in sorted(args.artifacts.glob('*/ui_checks.json'))]
    (args.artifacts / 'summary.json').write_text(json.dumps(summaries, indent=2), encoding='utf-8')
    for row in summaries:
        print(json.dumps(dict(run=row['run'], result=row['result']['outcome'],
            seconds=row['result']['elapsed_seconds'],
            first_preview_seconds=row['first_visible_preview_seconds'],
            maximum_gap_seconds=row['heartbeat_gap_max'],
            accepted_sections=row['accepted_sections'], tested_sections=row['tested_prepared_sections'],
            counts=row['counts'], peak_private_gib=row['sampled_peak_private_gib'])))


if __name__ == '__main__':
    main()
