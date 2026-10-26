"""Native grouped/per-Street/mixed authored inputs produce identical output."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
import bpy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT / 'tests'))
from check_python_street_scene import mesh_buffers
from check_python_street_reopen import equal
from _python_street_fixtures import fixture_networks
from blosm.way.street_network.python_backend import model
from blosm.way.street_network.python_backend.scene import NetworkScene


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
    args.output.mkdir(parents=True, exist_ok=True)
    network = fixture_networks(model)[-1]
    controllers = [NetworkScene.create(replace(network, uuid='grouping/' + mode), grouping=grouping)
                   for mode, grouping in (('grouped', 'grouped'), ('street', 'street'), ('mixed', 'street'))]
    baseline = mesh_buffers(controllers[0])
    for controller in controllers[1:]:
        assert equal(baseline, mesh_buffers(controller))
    mixed = controllers[2]
    sources = list(mixed['st_py_sources'].objects)
    original_count = len(sources)
    bpy.ops.object.select_all(action='DESELECT')
    for obj in sources[:2]:
        obj.select_set(True)
    bpy.context.view_layer.objects.active = sources[0]
    assert bpy.ops.object.join() == {'FINISHED'}
    assert len(mixed['st_py_sources'].objects) == original_count - 1
    update = NetworkScene.update(mixed)
    assert update.rebuilt == (), 'Regrouping changed semantic inputs'
    assert equal(baseline, mesh_buffers(mixed))
    member = next(o for o in mixed['st_py_sources'].objects
                  if 1 in {v.value for v in o.data.attributes['st_section_id'].data})
    # Locate Section 1 by its persistent edge ID and endpoint role, not the
    # merged object's name or its incidental vertex-array offset.
    edge = next(e for e in member.data.edges if member.data.attributes['st_section_id'].data[e.index].value == 1)
    role = member.data.attributes['st_endpoint_role']
    vertex = next(member.data.vertices[i] for i in edge.vertices if role.data[i].value == 0)
    original = vertex.co.copy()
    vertex.co.y += .2
    member.data.update()
    assert equal(baseline, mesh_buffers(mixed)), 'Authored geometry updated output without Update'
    NetworkScene.update(mixed)
    assert not equal(baseline, mesh_buffers(mixed))
    vertex.co = original
    member.data.update()
    NetworkScene.update(mixed)
    assert equal(baseline, mesh_buffers(mixed))
    report = dict(passed=True, grouped_per_street_equivalent=True, mixed_equivalent=True,
                  regrouping_preserves_semantic_identity=True, mixed_edit_update=True,
                  exact_restoration=True)
    (args.output / 'grouping_checks.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
