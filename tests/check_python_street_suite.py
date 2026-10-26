"""Reproducible pure-Python backend, preparation and asset regression suite.

Run with the host Python plus NumPy, without importing the Blender addon entry
point. Native publication/reopening checks remain separate Blender scripts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BASE_MODULES = (
    'test_aggregate_surface_coverage', 'test_asset_loading_policy', 'test_asset_recovery',
    'test_boundary_junctions', 'test_bounded_bundle_network', 'test_bundle_junctions',
    'test_bundle_topology', 'test_clipped_source_context', 'test_compound_network',
    'test_deferred_dimensions', 'test_intake_programs', 'test_junction_interiors',
    'test_junction_origins', 'test_library_store', 'test_logical_boundary_program',
    'test_longitudinal_separators', 'test_original_boundary_runs', 'test_original_span_coverage',
    'test_original_span_free_bundles', 'test_original_span_plan', 'test_original_span_separator_chains',
    'test_original_span_separator_continuations', 'test_original_span_separator_cycles',
    'test_original_span_separator_regions', 'test_original_span_separators', 'test_patch_program_catalog',
    'test_physical_approaches', 'test_physical_clusters', 'test_physical_connections',
    'test_physical_dimensions', 'test_physical_junctions', 'test_physical_ports', 'test_physical_sources',
    'test_port_clusters', 'test_port_compounds', 'test_port_connections', 'test_port_network',
    'test_port_source_layout', 'test_provider_evaluation', 'test_selection_boundaries',
    'test_shared_port_clusters', 'test_street_coverage', 'test_street_surface_validation',
)


def source_digest():
    digest = hashlib.sha256()
    for path in sorted(ROOT.rglob('*.py')):
        digest.update(path.relative_to(ROOT).as_posix().encode('utf8') + b'\0')
        digest.update(path.read_bytes() + b'\0')
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    temporary = args.output / 'temporary'
    temporary.mkdir(exist_ok=True)
    tempfile.tempdir = str(temporary.resolve())
    sys.path.insert(0, str(ROOT / 'tests'))
    names = sorted(set(BASE_MODULES) | {p.stem for p in (ROOT / 'tests').glob('test_python_street_*.py')})
    before = source_digest()
    with (args.output / 'tests.log').open('w', encoding='utf8') as stream:
        suite = unittest.defaultTestLoader.loadTestsFromNames(names)
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    after = source_digest()
    report = dict(passed=result.wasSuccessful() and before == after, tests=result.testsRun,
                  modules=names, failures=len(result.failures), errors=len(result.errors),
                  unchanged_sources=before == after, source_sha256_before=before,
                  source_sha256_after=after)
    (args.output / 'suite_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps({key: value for key, value in report.items() if key != 'modules'}), flush=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
