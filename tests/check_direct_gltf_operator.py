"""Check the actual Blosm import operator with a deterministic mixed tileset."""
import argparse
import importlib
import json
from math import radians
from pathlib import Path
import sys
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gltf_test_support import (
    bpy, REPO, fixture, make_glb, make_b3dm, reset_scene, no_tile_writes_or_operator,
)
import addon_utils

parser = argparse.ArgumentParser()
parser.add_argument("--artifacts", type=Path, required=True)
args = parser.parse_args(sys.argv[sys.argv.index("--")+1:])
args.artifacts.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO.parent))
package = REPO.name
addon_utils.enable(package, default_set=True)
try:
    assert package in bpy.context.preferences.addons
    bpy.context.preferences.addons[package].preferences.dataDir = str(args.artifacts)
    reset_scene()
    settings = bpy.context.scene.blosm
    settings.dataType = "3d-tiles"
    settings.threedTilesSource = "custom"
    settings.threedTilesUrl = "https://fixture.invalid/tileset.json"
    settings.minLon, settings.maxLon = 24.93, 24.95
    settings.minLat, settings.maxLat = 60.16, 60.18
    settings.join3dTilesObjects = False
    settings.cacheJsonFiles = False
    settings.cache3dFiles = True  # A cache miss must still never save a tile.
    manager_module = importlib.import_module(package + ".threed_tiles.manager")
    center = manager_module.BaseManager.fromGeographic(60.17, 24.94, 0)
    doc, binary, glb = fixture()
    b3dm = make_b3dm(glb, list(center))
    doc["nodes"][0]["translation"] = [center.x, center.z, -center.y]
    glb = make_glb(doc, binary)
    bounds = {"region": [radians(24.93), radians(60.16),
                         radians(24.95), radians(60.18), -100, 100]}
    tileset = {"asset": {"version": "1.0"}, "geometricError": 0, "root": {
        "geometricError": 0, "boundingVolume": bounds, "children": [
            {"geometricError": 0, "boundingVolume": bounds, "content": {"uri": "test.glb"}},
            {"geometricError": 0, "boundingVolume": bounds, "content": {"uri": "test.b3dm"}},
        ]}}
    payloads = {"tileset.json": json.dumps(tileset).encode(), "test.glb": glb, "test.b3dm": b3dm}
    requested = []
    def download(self, uri):
        name = uri.rsplit("/", 1)[-1]
        requested.append(name)
        return payloads[name]
    renderer_module = importlib.import_module(package + ".threed_tiles.blender")
    real_import = renderer_module.import_gltf
    from io_scene_gltf2.blender.imp.blender_gltf import BlenderGlTF
    real_create = BlenderGlTF.create
    forwarded = []
    def import_with_report_check(filepath, *, content=None, report=None):
        assert report is not None, "Blosm did not supply its report callback"
        def forward(level, message):
            forwarded.append((sorted(level), message))
            report(level, message)
        return real_import(filepath, content=content, report=forward)
    def create_with_warning(importer):
        importer.log.warning("fixture glTF warning", popup=True)
        real_create(importer)
    with patch.object(manager_module.BaseManager, "download", download), \
            patch.object(renderer_module, "import_gltf", side_effect=import_with_report_check), \
            patch.object(BlenderGlTF, "create", side_effect=create_with_warning), \
            no_tile_writes_or_operator():
        result = bpy.ops.blosm.import_data()
    assert result == {"FINISHED"}
    meshes = [obj for obj in bpy.context.scene.objects if obj.type == "MESH"]
    assert len(meshes) == 2
    assert all(len(obj.data.vertices) == 3 and len(obj.data.polygons) == 1 for obj in meshes)
    assert "City of Test" in settings.copyright and "Café Example" in settings.copyright
    assert requested == ["tileset.json", "test.glb", "test.b3dm"]
    assert forwarded == [(["WARNING"], "fixture glTF warning")] * 2
    successful_reports = list(forwarded)

    # A failed tile must report its cause and enter Blosm's error accounting,
    # while a later tile is still imported. Blender's Python operator wrapper
    # raises for ERROR reports even when the operator finishes a partial import.
    reset_scene()
    payloads["test.glb"] = b"invalid glTF JSON"
    requested.clear()
    forwarded.clear()
    tile_errors = []
    real_process_error = manager_module.BaseManager.processError
    def process_error(self, error, uri):
        tile_errors.append((str(error), uri))
        real_process_error(self, error, uri)
    with patch.object(manager_module.BaseManager, "download", download), \
            patch.object(manager_module.BaseManager, "processError", process_error), \
            patch.object(renderer_module, "import_gltf", side_effect=import_with_report_check), \
            no_tile_writes_or_operator():
        try:
            bpy.ops.blosm.import_data()
        except RuntimeError as exc:
            assert len(tile_errors) == 1 and tile_errors[0][0] in str(exc)
        else:
            raise AssertionError("Expected the reported glTF error from Blender's operator wrapper")
    assert len(tile_errors) == 1 and tile_errors[0][1].endswith("test.glb")
    assert forwarded == [(["ERROR"], tile_errors[0][0])]
    assert len([obj for obj in bpy.context.scene.objects if obj.type == "MESH"]) == 1
    assert requested == ["tileset.json", "test.glb", "test.b3dm"]
    assert not list(args.artifacts.rglob("*.glb"))
    assert not list(args.artifacts.rglob("*.b3dm"))
    report = {"blender": bpy.app.version_string, "status": "passed",
              "meshes": len(meshes), "requests": requested,
              "tile_files_written": 0, "gltf_operator_calls": 0,
              "forwarded_gltf_reports": successful_reports,
              "failed_tile_reports": forwarded, "continued_after_tile_error": True,
              "copyright": settings.copyright}
    (args.artifacts / "operator_check.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("BLOSM_OPERATOR_CHECK", json.dumps(report))
finally:
    addon_utils.disable(package, default_set=True)
