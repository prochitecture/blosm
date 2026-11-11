"""Mapillary-specific producer/consumer over the production Job and transport."""
from collections import Counter
from functools import partial
import hashlib
from itertools import islice
import json
import time
from urllib.parse import urlencode

import bpy
from mathutils import Matrix
from ..downloads import Request
from ..downloads.diagnostics import redact_text
from ..downloads.transport import DownloadTransport
from ..progressive.core import Job, Download, Drain
from ..progressive.blender import notify_result, preserve_context
from . import pose
from .tiles import MAX_TILE_BYTES, contains, identifier, image_points, tile_coordinates, validate_bounds

TILE_URL = 'https://tiles.mapillary.com/maps/vtp/mly1_computed_public/2/{z}/{x}/{y}'
GRAPH_URL = 'https://graph.mapillary.com/images'
# GraphQueryCreator in MapillaryJS uses images?image_ids=...&fields=... .
# 32 is our bounded working batch, not a claimed provider maximum.
METADATA_BATCH = 32
MAX_METADATA_BYTES = 256 * 1024
MAX_IMAGES = 100000
FIELDS = ('id', 'geometry', 'computed_geometry', 'altitude', 'computed_altitude',
          'computed_rotation', 'compass_angle', 'computed_compass_angle', 'camera_type',
          'camera_parameters', 'width', 'height', 'exif_orientation', 'sequence', 'captured_at', 'creator')


class MapillaryImport:
    def __init__(self, context, bounds, projection, origin, set_origin, token, *, transport=None,
                 tile_url=TILE_URL, graph_url=GRAPH_URL, max_images=MAX_IMAGES):
        validate_bounds(bounds)
        if not token or any(c.isspace() for c in token):
            raise ValueError('Set a valid Mapillary client access token in Blosm preferences')
        self.scene, self.view_layer = context.scene, context.view_layer
        self.bounds, self.projection, self.origin = tuple(bounds), projection, tuple(origin)
        self.set_origin = set_origin
        self.token = token
        self.auth_scope = hashlib.sha256(token.encode()).hexdigest()
        self.headers = {'Authorization': 'OAuth ' + token}
        self.tile_url, self.graph_url = tile_url, graph_url
        self.max_images = max_images
        self.collection = None
        self.count = self.revision = self.approximate = 0
        self.errors = Counter()
        self.details = []
        self.seen = set()
        self.metrics = {'decode_seconds': 0., 'metadata_seconds': 0., 'camera_seconds': 0.,
                        'max_decode_turn': 0., 'max_metadata_parse': 0., 'max_camera_commit': 0.,
                        'peak_tile_bytes': 0, 'peak_metadata_records': 0}
        self.metrics['peak_dedup_ids'] = 0
        offset = pose.number(self.scene.get('height_offset'))
        self.height_source = 'scene' if offset is not None else 'absolute_provider_altitude'
        if offset is None:
            terrain = self.scene.objects.get(self.scene.blosm.terrainObject)
            offset = pose.number(terrain.get('height_offset')) if terrain else None
            if offset is not None:
                self.height_source = 'terrain'
        self.height_offset = offset if offset is not None else 0.
        settings = self.scene.blosm
        workers = settings.downloadWorkers if settings.concurrentDownloads else 1
        self.transport = transport if transport is not None else DownloadTransport(
            max_requests=workers, max_bytes=MAX_TILE_BYTES,
            max_buffer_bytes=workers*MAX_TILE_BYTES)
        self.job = Job(self.steps(), self.transport, validate=self.validate)

    def validate(self):
        if self.collection is not None:
            try:
                if self.collection in self.scene.collection.children_recursive:
                    return
            except ReferenceError:
                pass
            raise RuntimeError('The Mapillary destination collection was deleted or unlinked from the import scene')

    def request(self, url, *, tile=False):
        if tile:
            url += ('&' if '?' in url else '?') + urlencode({'access_token': self.token})
        return Request(url, headers={} if tile else self.headers, decode_gzip=True,
                       auth_scope=self.auth_scope)

    def failure(self, stage, error):
        # Transport errors are sanitized already, but also sanitize parse errors
        # and any future adapter diagnostics against this import's credentials.
        self.errors[stage] += 1
        if len(self.details) < 20:
            self.details.append(redact_text('%s: %s' % (stage, error), headers=self.headers))
        if getattr(error, 'code', None) == 401 or getattr(error, 'auth_failed', False):
            raise RuntimeError('Mapillary authentication failed; check the client access token')
        if getattr(error, 'code', None) in (400, 403):
            raise RuntimeError('Mapillary rejected the request (HTTP %s); check API access and request fields' % error.code)

    def _coverage(self, output, tile, content, error):
        if error is not None:
            self.failure('coverage', error)
        else:
            output.append((tile, content))

    def steps(self):
        # At most one window of raw tiles, one tile dictionary iterator and one
        # bounded window of metadata replies. Decoded dataset records never form
        # a whole-area list. Global ID deduplication has an explicit ceiling.
        tiles = iter(tile_coordinates(self.bounds))
        raw = []
        try:
            while group := list(islice(tiles, self.transport.max_requests)):
                for tile in group:
                    z, x, y = tile
                    yield Download(self.request(self.tile_url.format(z=z, x=x, y=y), tile=True),
                                   partial(self._coverage, raw, tile))
                yield Drain()
                self.metrics['peak_tile_bytes'] = max(self.metrics['peak_tile_bytes'], sum(len(b) for _, b in raw))
                while raw:
                    tile, content = raw.pop()
                    records = []
                    iterator = image_points(content, tile)
                    processed = 0
                    checkpoint = time.monotonic()
                    try:
                        while True:
                            started = time.monotonic()
                            try:
                                record = next(iterator)
                            except StopIteration:
                                break
                            elapsed = time.monotonic()-started
                            self.metrics['decode_seconds'] += elapsed
                            self.metrics['max_decode_turn'] = max(self.metrics['max_decode_turn'], elapsed)
                            processed += 1
                            if record is not None and contains(self.bounds, *record['coverage_geometry']):
                                identity = record['id']
                                if identity not in self.seen:
                                    if len(self.seen) >= self.max_images:
                                        raise RuntimeError('Mapillary image limit reached; select a smaller area')
                                    self.seen.add(identity)
                                    self.metrics['peak_dedup_ids'] = len(self.seen)
                                    records.append(record)
                                    if len(records) == METADATA_BATCH:
                                        yield self.metadata_download(records)
                                        records = []
                            # A dense city tile contains many points outside a
                            # small ROI. Keep traversal cooperative without an
                            # entire UI timer turn per few dozen rejected points.
                            if record is None or processed % 16 == 0 or time.monotonic()-checkpoint >= .002:
                                yield None
                                checkpoint = time.monotonic()
                    except ValueError as error:
                        self.failure('coverage_decode', error)
                    finally:
                        iterator.close()
                        del iterator, content
                    if records:
                        yield self.metadata_download(records)
                    records = None
                    yield Drain()
            return self.count
        finally:
            raw.clear()
            self.seen.clear()
            # No credentials belong to the persisted Blender scene or reports.
            # Retain credentials only in memory through finish(), so even an
            # unexpected fatal exception can be redacted before reporting it.
            self.token = ''

    def metadata_download(self, records):
        url = self.graph_url + '?' + urlencode({'image_ids': ','.join(r['id'] for r in records),
                                               'fields': ','.join(FIELDS)})
        return Download(self.request(url), partial(self._metadata, tuple(records)))

    def _metadata(self, records, content, error):
        # A cooperative consumer, driven by Job.step(), never by paused poll().
        metadata = {}
        if error is not None:
            self.failure('metadata_request', error)
        else:
            started = time.monotonic()
            try:
                if len(content) > MAX_METADATA_BYTES:
                    raise ValueError('Mapillary metadata reply exceeds 256 KiB per 32-image batch')
                document = json.loads(content)
                if isinstance(document, dict) and 'error' in document:
                    raise RuntimeError('Mapillary returned a Graph API error; check API access and requested fields')
                if not isinstance(document, dict):
                    raise ValueError('Invalid Mapillary metadata document')
                data = document.get('data')
                if not isinstance(data, list) or len(data) > METADATA_BATCH:
                    raise ValueError('Invalid Mapillary metadata batch')
                wanted = {record['id'] for record in records}
                for entry in data:
                    if not isinstance(entry, dict):
                        raise ValueError('Invalid Mapillary image entity')
                    identity = identifier(entry.get('id'))
                    if identity in wanted:
                        metadata[identity] = entry
                self.metrics['peak_metadata_records'] = max(self.metrics['peak_metadata_records'], len(metadata))
            except (ValueError, TypeError, RecursionError) as failure:
                self.failure('metadata_decode', failure)
            elapsed = time.monotonic()-started
            self.metrics['metadata_seconds'] += elapsed
            self.metrics['max_metadata_parse'] = max(self.metrics['max_metadata_parse'], elapsed)
        content = None  # Release wire bytes before committing camera objects.
        yield None
        for record in records:
            entry = metadata.pop(record['id'], {})
            if not entry:
                self.errors['missing_metadata'] += 1
            try:
                camera_pose = pose.interpret(entry, record, self.projection, self.height_offset)
                if contains(self.bounds, camera_pose.properties['source_longitude'], camera_pose.properties['source_latitude']):
                    self.commit(record, entry, camera_pose)
            except (ValueError, TypeError, OverflowError) as failure:
                self.failure('camera_metadata', failure)
            yield None

    def commit(self, record, metadata, camera_pose):
        self.validate()
        started = time.monotonic()
        camera = obj = None
        with preserve_context():
            try:
                if self.collection is None:
                    self.collection = bpy.data.collections.new('Mapillary cameras')
                    self.scene.collection.children.link(self.collection)
                    self.collection['source'] = 'Mapillary'
                    self.collection['source_url'] = 'https://www.mapillary.com/'
                camera = bpy.data.cameras.new('Mapillary ' + record['id'])
                camera.type = camera_pose.camera_type
                camera.lens = camera_pose.lens
                camera.sensor_fit = 'HORIZONTAL'
                camera.sensor_width = 36.
                camera.sensor_height = camera_pose.sensor_height
                camera.display_size = 0.75
                if camera.type == 'PANO' and hasattr(camera, 'panorama_type'):
                    camera.panorama_type = 'EQUIRECTANGULAR'
                obj = bpy.data.objects.new('Mapillary ' + record['id'], camera)
                obj.location = camera_pose.position
                obj.rotation_mode = 'QUATERNION'
                obj.rotation_quaternion = Matrix(camera_pose.rotation).to_quaternion()
                props = dict(camera_pose.properties)
                props.update(image_id=record['id'], sequence_id=str(metadata.get('sequence') or record.get('sequence_id') or ''),
                             captured_at=str(metadata.get('captured_at') or record.get('captured_at') or ''),
                             image_url='https://www.mapillary.com/app/?pKey=' + record['id'],
                             coverage_longitude=record['coverage_geometry'][0], coverage_latitude=record['coverage_geometry'][1],
                             height_offset_source=self.height_source,
                             metadata_available=bool(metadata))
                for key, value in props.items():
                    if isinstance(value, str):
                        value = redact_text(value, headers=self.headers)
                    obj['mapillary_' + key] = json.dumps(value) if isinstance(value, (list, dict)) else value
                # Preserve only requested metadata, redacted even if a provider
                # unexpectedly echoes a request credential. Never retain URLs.
                safe = {key: metadata[key] for key in FIELDS if key in metadata}
                obj['mapillary_metadata'] = redact_text(json.dumps(safe), headers=self.headers)
                self.collection.objects.link(obj)
                if self.set_origin:
                    self.scene['lat'], self.scene['lon'] = self.origin
                self.count += 1
                self.approximate += int(camera_pose.properties['pose_approximate'])
                self.revision += 1
            except Exception:
                if obj is not None:
                    bpy.data.objects.remove(obj, do_unlink=True)
                if camera is not None and camera.users == 0:
                    bpy.data.cameras.remove(camera)
                raise
        elapsed = time.monotonic()-started
        self.metrics['camera_seconds'] += elapsed
        self.metrics['max_camera_commit'] = max(self.metrics['max_camera_commit'], elapsed)

    @property
    def status(self):
        return '%s Mapillary cameras, %s approximate, %s errors%s' % (
            self.count, self.approximate, sum(self.errors.values()), ' — downloading' if self.job.waiting else '')

    def cancel(self):
        self.job.cancel()

    def finish(self, silent=False):
        for detail in self.details:
            print('Mapillary:', detail)
        if silent:
            self.headers.clear()
            return
        summary = '%s cameras imported (%s approximate); %s errors.' % (
            self.count, self.approximate, sum(self.errors.values()))
        if self.job.error:
            notify_result('ERROR', 'Mapillary import stopped: %s. %s' % (
                redact_text(self.job.error, headers=self.headers), summary))
        elif self.job.cancelled:
            notify_result('INFO', 'Mapillary import cancelled. Kept ' + summary)
        else:
            notify_result('INFO' if self.count else 'ERROR', 'Mapillary import complete. ' + summary)
        self.headers.clear()

    def run_blocking(self):
        # Same producer, transport and scene consumer as progressive execution.
        with preserve_context():
            while not self.job.done:
                self.job.step()
                if self.job.waiting:
                    time.sleep(.001)
        self.finish()
        return self.job.error is None
