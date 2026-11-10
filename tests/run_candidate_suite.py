"""Install the candidate before running an existing feature suite in Blender.

The standalone mode extracts identical ZIP bytes for tests that need ordinary
Python, without importing bpy. Both routes prohibit runtime code from the checkout.
"""
import argparse
import hashlib
import importlib
from importlib.abc import MetaPathFinder
from importlib.machinery import PathFinder
import json
import os
from pathlib import Path
import runpy
import sys
import zipfile

sys.dont_write_bytecode = True
TESTS = Path(__file__).resolve().parent
SOURCE = TESTS.parent
sys.path.insert(0, str(TESTS))
from release_profiles import inspect_candidate


class SourceIsolation(MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        spec = PathFinder.find_spec(fullname, path)
        if spec and spec.origin and spec.origin not in ('built-in', 'frozen'):
            origin = Path(spec.origin).resolve()
            if origin.is_relative_to(SOURCE) and not origin.is_relative_to(TESTS):
                raise ImportError(f'Candidate test attempted to load checkout module: {fullname}')
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate-zip', required=True, type=Path)
    parser.add_argument('--profile', required=True, choices=('base', 'pro'))
    parser.add_argument('--suite', required=True)
    parser.add_argument('--artifacts', required=True, type=Path)
    parser.add_argument('suite_args', nargs=argparse.REMAINDER)
    argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv and 'bpy' in sys.modules else sys.argv[1:]
    args = parser.parse_args(argv)
    args.artifacts.mkdir(parents=True, exist_ok=True)
    manifest, names = inspect_candidate(args.candidate_zip)
    os.environ['BLOSM_TEST_PROFILE'] = args.profile
    # The wrapper needs only test helpers; no addon path from the caller is used.
    sys.path[:] = [p for p in sys.path if Path(p or '.').resolve() not in (SOURCE, SOURCE.parent)]
    sys.meta_path.insert(0, SourceIsolation())
    installed = args.suite != 'standalone'
    if installed:
        import bpy
        import addon_utils
        repo_dir = args.artifacts / 'installed'
        repo_dir.mkdir()
        repo = bpy.context.preferences.extensions.repos.new(
            name='Blosm release tests', module='blosm_release_tests', custom_directory=str(repo_dir))
        package = f'bl_ext.{repo.module}.{manifest["id"]}'
        assert bpy.ops.extensions.package_install_files(
            filepath=str(args.candidate_zip), repo=repo.module, enable_on_install=True) == {'FINISHED'}
        addon = importlib.import_module(package)
        root = Path(addon.__file__).resolve().parent
        assert root == (repo_dir / manifest['id']).resolve(), 'Candidate was not installed in the isolated repository'
        assert addon_utils.check(package)[1]
        assert bool(addon.blenderApp.app.isPro) == (args.profile == 'pro'), 'Candidate edition does not match profile'
        os.environ['BLOSM_TEST_PACKAGE'] = package
    else:
        root = args.artifacts / 'candidate'
        root.mkdir()
        with zipfile.ZipFile(args.candidate_zip) as archive:
            archive.extractall(root)  # inspect_candidate rejected traversal and symlinks.
        package = None
        os.environ.pop('BLOSM_TEST_PACKAGE', None)
    os.environ['BLOSM_TEST_ROOT'] = str(root)
    metadata = dict(profile=args.profile, package=package, installed=installed,
                    sha256=hashlib.sha256(args.candidate_zip.read_bytes()).hexdigest(),
                    root=str(root), manifest=manifest, file_count=len(names))
    (args.artifacts / 'candidate.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    remaining = args.suite_args[1:] if args.suite_args[:1] == ['--'] else args.suite_args
    if installed:
        script = TESTS / args.suite
        assert script.is_file() and script.parent == TESTS
        sys.argv = [str(script), '--', '--artifacts', str(args.artifacts), *remaining]
        runpy.run_path(str(script), run_name='__main__')
    else:
        sys.argv = ['unittest', 'discover', '-s', str(TESTS), '-p', 'test_*.py', '-v']
        runpy.run_module('unittest', run_name='__main__')


if __name__ == '__main__':
    main()
