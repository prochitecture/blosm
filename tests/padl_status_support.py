"""Private GN status resources may live only inside an active import batch."""
import bpy


def assert_status_resources(runtime):
    session = getattr(runtime.app, 'osm_session', None)
    active = session is not None and not session.job.done and not session.finished
    groups = [item for item in bpy.data.node_groups if item.name.startswith('PADL temporary status reader')]
    objects = [item for item in bpy.data.objects if item.name.startswith('PADL temporary status owner')]
    scenes = [item for item in bpy.data.scenes if item.name.startswith('PADL private status')]
    meshes = [item for item in bpy.data.meshes if item.name.startswith('PADL temporary status empty mesh')]
    if not active:
        assert not groups + objects + scenes + meshes, 'Status evaluator survived batch teardown/history'
    elif groups or objects or scenes or meshes:
        assert len(groups) == len(objects) == len(scenes) == len(meshes) == 1
        assert all(item.get('_blosm_transient_import') for item in groups + objects + scenes + meshes)
        assert scenes[0] != session.scene and scenes[0] != bpy.context.scene
        assert objects[0].name not in session.scene.objects, 'Status helper leaked into import scene'
        assert objects[0].name in scenes[0].objects
