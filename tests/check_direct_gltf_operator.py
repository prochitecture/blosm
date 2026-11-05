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
    with patch.object(manager_module.BaseManager, "download", download), no_tile_writes_or_operator():
        result = bpy.ops.blosm.import_data()
    assert result == {"FINISHED"}
    meshes = [obj for obj in bpy.context.scene.objects if obj.type == "MESH"]
    assert len(meshes) == 2
    assert all(len(obj.data.vertices) == 3 and len(obj.data.polygons) == 1 for obj in meshes)
    assert "City of Test" in settings.copyright and "Café Example" in settings.copyright
    assert requested == ["tileset.json", "test.glb", "test.b3dm"]
    assert not list(args.artifacts.rglob("*.glb"))
    assert not list(args.artifacts.rglob("*.b3dm"))
    report = {"blender": bpy.app.version_string, "status": "passed",
              "meshes": len(meshes), "requests": requested,
              "tile_files_written": 0, "gltf_operator_calls": 0,
              "copyright": settings.copyright}
    (args.artifacts / "operator_check.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("BLOSM_OPERATOR_CHECK", json.dumps(report))
finally:
    addon_utils.disable(package, default_set=True)
