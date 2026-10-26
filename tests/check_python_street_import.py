"""Normal OSM preparation and Python backend qualification in background Blender."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import bpy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--osm', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--capture-separator-boundaries', action='store_true')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    captures = []
    if args.capture_separator_boundaries:
        # Diagnostic-only observer. It preserves numerical arguments/results
        # and captures a small, reusable chart input before a possible failure.
        import numpy as np
        from blosm.way.street_network.python_backend.separator_surfaces import FixedRimSeparator
        original = FixedRimSeparator.build
        def captured(self, boundary, **kwargs):
            path = args.output / f'separator_boundary_{len(captures):02d}.npz'
            np.savez_compressed(path, boundary=boundary, **kwargs)
            captures.append(str(path))
            print(json.dumps(dict(separator_boundary_capture=str(path))), flush=True)
            return original(self, boundary, **kwargs)
        FixedRimSeparator.build = captured
    from blosm.app import BaseApp
    from blosm.app.blender import BlenderApp
    from blosm.parse.osm import Osm
    from blosm.setup import Setup
    from blosm.action.generate_streets import StreetGenerator
    from blosm.way.renderer_streets import StreetRenderer
    app = BlenderApp()
    app.mode = BaseApp.realistic
    app.loadMissingMembers = False
    app.straightAngleThreshold = 175.5
    app.streetBackend = 'python'
    app.osmFilepath = str(args.osm.resolve())
    app.assetsDir = ''
    app.pmlFilepathStreet = None
    app.trafficSide = 'auto'
    app.assetLibraryMode = 'APPEND'
    BaseApp.initOsm(app)
    osm = Osm(app)
    setup = Setup(app, osm)
    setup.skipWays()
    setup.roadsAndPaths()
    manager = setup.getWayManager()
    generator = StreetGenerator(manager.styleStore, getStyle=manager.getStyle, doPixelCorrection=False)
    manager.addAction(generator)
    start = time.perf_counter()
    app.createLayers(osm)
    osm.parse(str(args.osm.resolve()), forceExtentCalculation=True)
    app.initLayers()
    app.process()
    prepare_seconds = time.perf_counter() - start
    renderer = StreetRenderer(app)
    renderer.prepare()
    context = renderer.render(manager, osm)
    renderer.finalize()
    renderer.cleanup()
    report = context.report
    report.update(input=str(args.osm.resolve()), input_sha256=hashlib.sha256(args.osm.read_bytes()).hexdigest(),
                  preparation_seconds=prepare_seconds, blender=bpy.app.version_string,
                  gn_modifiers=sum(m.type == 'NODES' for o in bpy.data.objects for m in o.modifiers))
    report['separator_boundary_captures'] = captures
    from blosm.way.street_network.python_backend.scene import SourceObjects
    controllers = [o for o in bpy.data.objects if o.get('st_py_source_data')]
    templates = [json.loads(o['st_py_source_data']) for o in controllers]
    numerical = [SourceObjects.read(o).to_dict() for o in controllers]
    (args.output / 'prepared_input_templates.json').write_text(json.dumps(templates, indent=2), encoding='utf8')
    (args.output / 'numerical_inputs.json').write_text(json.dumps(numerical, indent=2), encoding='utf8')
    assert report['gn_modifiers'] == 0
    bpy.context.preferences.filepaths.file_preview_type = 'NONE'
    scene = args.output / 'python_street_import.blend'
    bpy.ops.wm.save_as_mainfile(filepath=str(scene.resolve()))
    (args.output / 'import_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(dict(scene=str(scene), coverage=report['coverage'],
        ready=sum(r['status'] == 'READY' for r in report['surfaces']), components=len(report['surfaces']),
        quads=sum(r.get('quads', 0) for r in report['surfaces']),
        failures=[r.get('reason') for r in report['surfaces'] if r['status'] != 'READY'])), flush=True)


if __name__ == '__main__':
    main()
