"""Build a release ZIP from Git-visible source, excluding ignored local files."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--blender', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    paths = subprocess.check_output(
        ['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard'], cwd=root
    ).decode('utf-8').split('\0')
    # Stage in the output directory so development outputs stay outside source.
    with tempfile.TemporaryDirectory(prefix='extension-source-', dir=output) as staging:
        staging = Path(staging)
        for name in sorted(set(paths)):
            if not name or not (root / name).is_file():
                continue
            target = staging / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / name, target)
        subprocess.run([str(args.blender), '--background', '--factory-startup',
                        '--command', 'extension', 'build', '--source-dir', str(staging),
                        '--output-dir', str(output)], check=True)


if __name__ == '__main__':
    main()
