"""Normal OSM import followed by physical UV verification on every mesh.

Run with the arguments of check_python_street_import.py. The UV native helper
selects this checkout explicitly, including when its worktree is not named blosm.
"""
import json
from pathlib import Path
import sys
import time
sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_street_uv as uv_check
import bpy
from check_python_street_import import main as import_main


def main():
    import_main()
    started = time.perf_counter()
    args = sys.argv[sys.argv.index('--') + 1:]
    output = Path(args[args.index('--output') + 1])
    controllers = [o for o in bpy.data.objects if o.get('st_py_network_uuid')]
    assert controllers, 'No generated street networks to verify'
    report = dict(passed=True, networks=[uv_check.verify_mapping(c) for c in controllers])
    report['seconds'] = time.perf_counter() - started
    (output / 'uv_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print('UV_IMPORT_OK', json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
