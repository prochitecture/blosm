"""Idle factory-Blender focus target for the isolated two-instance timing test.

No addon is enabled and no generation is run here. This owned test window exits
once the specified result receipt exists. Its profile and outputs must be in
this chat's artifacts; never use a normal user Blender session as the target.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import bpy


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--result', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    bpy.context.preferences.view.show_splash = False
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'

    def ready():
        windows = list(bpy.context.window_manager.windows)
        if not windows:
            return .25
        with bpy.context.temp_override(window=windows[0]):
            bpy.ops.wm.save_as_mainfile(filepath=str(args.output / 'focus_target.blend'))
        (args.output / 'ready.json').write_text(json.dumps(dict(pid=os.getpid(), role='idle isolated focus target')))
        bpy.app.timers.register(watch_result, first_interval=1.)
        return None

    def watch_result():
        if not args.result.exists():
            return 1.
        bpy.ops.wm.quit_blender()
        return None

    bpy.app.timers.register(ready, first_interval=1.)


if __name__ == '__main__':
    main()
