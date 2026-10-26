"""Read effective authored inputs from a saved scene, including failed imports.

Stored templates can precede Blender float32 transport or later source edits.
Export both; geometry diagnosis must consume SourceObjects.read() results.
"""
import argparse
import json
from pathlib import Path
import sys
import bpy
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
from blosm.way.street_network.python_backend.model import Network
from blosm.way.street_network.python_backend.scene import SourceObjects


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scene', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    bpy.ops.wm.open_mainfile(filepath=str(args.scene.resolve()))
    args.output.mkdir(parents=True, exist_ok=True)
    inputs, templates, records = [], [], []
    for controller in bpy.data.objects:
        if not controller.get('st_py_source_data'):
            continue
        template = json.loads(controller['st_py_source_data'])
        original = Network.from_dict(template)
        current = SourceObjects.read(controller)
        assert original.section_index.keys() == current.section_index.keys()
        assert current == Network.from_dict(json.loads(json.dumps(current.to_dict())))
        differences = []
        for sid, section in current.section_index.items():
            old = original.section_index[sid]
            differences.append(dict(section=sid,
                maximum_control_difference_m=float(np.max(abs(np.asarray(section.points) - old.points))),
                maximum_dimension_difference_m=max(abs(getattr(section, k) - getattr(old, k))
                                                    for k in ('width', 'offset', 'rise'))))
        record = dict(network=current.uuid, sections=len(current.sections), differences=differences,
                      selection_ends=sum(len(s.selection_ends) for s in current.sections))
        if controller.get('st_py_source_support'):
            record['initial_source_support'] = json.loads(controller['st_py_source_support'])
        inputs.append(current.to_dict())
        templates.append(template)
        records.append(record)
    (args.output / 'numerical_inputs.json').write_text(json.dumps(inputs, indent=2), encoding='utf8')
    (args.output / 'prepared_input_templates.json').write_text(json.dumps(templates, indent=2), encoding='utf8')
    report = dict(passed=True, scene=str(args.scene.resolve()), effective_authored_inputs=True,
                  generated_meshes_changed=False, networks=records)
    (args.output / 'input_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(dict(passed=True, networks=len(records), effective_authored_inputs=True,
                         generated_meshes_changed=False)), flush=True)


if __name__ == '__main__':
    main()
