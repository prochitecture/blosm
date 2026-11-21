"""Copy a coherent working tree into a disposable, immutable test candidate.

This is not a checkout, build, installation or recovery source. Maintained code
stays in the repository. Launch the copied tests in an isolated factory profile
so later development edits cannot mix old loaded modules with new source files.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def git(checkout, *args):
    return subprocess.check_output(['git', '-C', str(checkout), *args])


def hashes(checkout, names):
    return {name: hashlib.sha256((checkout/name).read_bytes()).hexdigest()
            for name in names}


def candidate_files(checkout):
    names = {name for name in git(checkout, 'ls-files', '--cached', '--others',
                                 '--exclude-standard', '-z').decode().split('\0') if name}
    for name in ('assets/prochitecture_street_materials.blend',
                 'assets/prochitecture_street_terrain.blend'):
        if (checkout/name).is_file():
            names.add(name)
    return sorted(name for name in names if name != 'AGENTS.md'
                  and '__pycache__' not in Path(name).parts
                  and (checkout/name).is_file())


def freeze(checkout, output):
    checkout, output = checkout.resolve(), output.resolve()
    artifact_root = Path('C:/dev/workspace/artifacts/streets').resolve()
    if not output.is_relative_to(artifact_root) or output == artifact_root:
        raise ValueError('Test snapshots belong in a dedicated street artifact folder')
    actual_root = Path(git(checkout, 'rev-parse', '--show-toplevel').decode().strip()).resolve()
    if actual_root != checkout:
        raise ValueError('The source must be an explicit repository checkout root')
    if any(entry.startswith(b'160000 ') for entry in git(checkout, 'ls-files', '--stage', '-z').split(b'\0')):
        raise ValueError('Submodule source needs explicit snapshot handling')
    names = candidate_files(checkout)
    for name in names:
        if not (checkout/name).resolve().is_relative_to(checkout):
            raise ValueError('Snapshot input leaves the source checkout: '+name)
    before = hashes(checkout, names)
    output.mkdir(parents=True, exist_ok=False)
    destination = output/'addon'
    destination.mkdir()
    manifest = dict(qualified=False, source_checkout=str(checkout),
                    snapshot=str(destination),
                    head=git(checkout, 'rev-parse', 'HEAD').decode().strip(),
                    branch=git(checkout, 'branch', '--show-current').decode().strip(),
                    files=before, copied_runtime_assets=True,
                    source_of_truth='The original repository, not this disposable candidate')
    target = output/'source_manifest.json'
    target.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    for name in names:
        copied = destination/name
        copied.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(checkout/name, copied)
    after = hashes(checkout, names)
    copied = hashes(destination, names)
    changed = [name for name in names if before[name] != after[name]
               or before[name] != copied[name]]
    changed.extend(sorted(set(names) ^ set(candidate_files(checkout))))
    manifest.update(qualified=not changed, changed_during_copy=changed,
                    file_count=len(names), bytes=sum((destination/name).stat().st_size for name in names))
    target.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    if changed:
        raise RuntimeError('Source changed while copying; do not launch this candidate')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkout', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = freeze(args.checkout, args.output)
    print(json.dumps({key: value for key, value in result.items() if key != 'files'}, indent=2))
