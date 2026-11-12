"""Probe a real area's poses through the production pipeline without mass cameras.

Use an isolated factory-startup Blender process. The token is read from a local
file, never from a command-line value. Generated reports and bounded metadata
samples belong in the supplied artifacts directory, outside the repository.
"""
import argparse
from collections import Counter
import importlib
import json
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))
import bpy

addon = importlib.import_module(ROOT.name)
addon.register()
adapter = importlib.import_module(ROOT.name + '.mapillary.progressive')
Projection = importlib.import_module(ROOT.name + '.util.transverse_mercator').TransverseMercator
sanitize = importlib.import_module(ROOT.name + '.downloads.diagnostics').redact_text

parser = argparse.ArgumentParser()
parser.add_argument('--osm', type=Path, required=True)
parser.add_argument('--token-file', type=Path, required=True)
parser.add_argument('--artifacts', type=Path, required=True)
args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:])
with args.osm.open('rb') as stream:
    for _, element in ET.iterparse(stream, events=('start',)):
        if element.tag == 'bounds':
            bounds = tuple(float(element.attrib[key]) for key in ('minlon', 'minlat', 'maxlon', 'maxlat'))
            break
    else:
        raise RuntimeError('The OSM fixture has no bounds')
args.artifacts.mkdir(parents=True, exist_ok=True)
token = args.token_file.read_text(encoding='utf-8-sig').strip()
if not token or any(c.isspace() for c in token):
    raise RuntimeError('Invalid local development token')
headers = {'Authorization': 'OAuth ' + token}
origin = ((bounds[1] + bounds[3])/2, (bounds[0] + bounds[2])/2)
projection = Projection(lat=origin[0], lon=origin[1])
bpy.context.scene.blosm.concurrentDownloads = True
bpy.context.scene.blosm.downloadWorkers = 4


class PoseProbe(adapter.MapillaryImport):
    def __init__(self):
        super().__init__(bpy.context, bounds, projection, origin, False, token)
        self.orientations = Counter()
        self.models = Counter()
        self.downward = 0
        self.minimum_up_z = 1.
        self.samples = {'up': [], 'down': []}

    def commit(self, record, metadata, camera_pose):
        # Keep discovery, filtering, metadata batching and pose interpretation
        # identical to production. Replace only the bpy camera allocation.
        up_z = camera_pose.rotation[2][1]
        direction = 'down' if up_z < 0 else 'up'
        self.orientations[(str(metadata.get('exif_orientation')), direction)] += 1
        self.models[camera_pose.properties['camera_model']] += 1
        self.minimum_up_z = min(self.minimum_up_z, up_z)
        self.downward += int(direction == 'down')
        self.count += 1
        self.approximate += int(camera_pose.properties['pose_approximate'])
        # Retain at most 32 affected entries and four ordinary reference entries.
        # No complete decoded dataset or set of camera datablocks accumulates.
        if len(self.samples[direction]) < (32 if direction == 'down' else 4):
            self.samples[direction].append({
                'record': record, 'metadata': metadata,
                'position': camera_pose.position, 'rotation': camera_pose.rotation,
                'properties': camera_pose.properties, 'up_z': up_z,
            })


task = PoseProbe()
started = last_report = time.monotonic()
try:
    while not task.job.done:
        task.job.step()
        now = time.monotonic()
        if now - last_report > 5:
            print('Checked %s image poses; %s downward' % (task.count, task.downward), flush=True)
            last_report = now
        if now - started > 600:
            task.cancel()
            raise RuntimeError('Area pose probe exceeded its ten-minute deadline')
        if task.job.waiting:
            time.sleep(.001)
finally:
    if not task.job.done:
        task.cancel()
    report = {
        'bounds': bounds, 'blender': bpy.app.version_string,
        'checked': task.count, 'downward': task.downward,
        'minimum_up_z': task.minimum_up_z if task.count else None,
        'orientation_counts': {'/'.join(k): v for k, v in task.orientations.items()},
        'camera_models': dict(task.models), 'approximate': task.approximate,
        'errors': dict(task.errors),
        'fatal_error': sanitize(task.job.error, headers=headers) if task.job.error else None,
        'cancelled': task.job.cancelled, 'elapsed_seconds': time.monotonic() - started,
        'metrics': task.metrics, 'samples': task.samples,
    }
    text = sanitize(json.dumps(report, indent=2), headers=headers)
    assert token not in text, 'Credential found in report'
    (args.artifacts/'area_pose.json').write_text(text, encoding='utf-8')
    print('Checked %s image poses; %s downward. Report: %s' % (
        task.count, task.downward, args.artifacts/'area_pose.json'), flush=True)
assert task.job.error is None, 'Area probe stopped; see sanitized report'
assert task.count, 'No image entries found'
addon.unregister()
