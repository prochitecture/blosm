"""Run host import contracts and the native extension discovery/import check.

Launch in factory-startup Blender with this chat's isolated user directories.
--repository contains a blosm link to this exact development checkout.
"""
import argparse
import json
from pathlib import Path
import runpy
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / 'tests'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--artifacts', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.artifacts.mkdir(parents=True, exist_ok=True)
    suite = unittest.defaultTestLoader.loadTestsFromNames(('test_extension_contract', 'test_vendor_imports'))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    (args.artifacts / 'host_checks.json').write_text(json.dumps(dict(
        passed=result.wasSuccessful(), count=result.testsRun,
        failures=len(result.failures), errors=len(result.errors)), indent=2), encoding='utf8')
    if not result.wasSuccessful():
        raise AssertionError('Extension import contracts failed')
    sys.argv = [str(ROOT / 'tests/check_development_extension.py'), '--',
                '--repository', str(args.repository), '--artifacts', str(args.artifacts / 'native')]
    runpy.run_path(str(ROOT / 'tests/check_development_extension.py'), run_name='__main__')


if __name__ == '__main__':
    main()
