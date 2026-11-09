"""Render an already-qualified UV gallery with Cycles CPU, without rebuilding."""
import argparse
import json
from pathlib import Path
import sys
import time

import bpy

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_street_uv import render_cpu


def main():
    default = Path(__file__).resolve().parents[2] / 'artifacts/streets/uv_mapping'
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, default=default / 'street_uv_mapping.blend')
    parser.add_argument('--output', type=Path, default=default)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--width', type=int, default=1000)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    args.output.mkdir(parents=True, exist_ok=True)
    if not 1 <= args.samples <= 256 or not 128 <= args.width <= 4096:
        parser.error('Samples must be 1–256 and width 128–4096 pixels')
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    started = time.perf_counter()
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()), use_scripts=False)
    render_cpu(bpy.context.scene, args.output, samples=args.samples, width=args.width)
    image = args.output / 'uv_mapping.png'
    assert image.is_file() and image.stat().st_size > 0
    report = dict(passed=True, engine='CYCLES', device='CPU', samples=args.samples, width=args.width,
                  scene=str(args.scene.resolve()), image=str(image.resolve()),
                  seconds=time.perf_counter() - started)
    (args.output / 'render_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print('UV_RENDER_OK', json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
