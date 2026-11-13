"""OSM session shared by blocking and timer-driven imports.

Managers may implement process_steps(session), render_steps(session), and
render_extra_steps(session). A bulk producer chooses its own atomic units with
session.unit(callable); it is never required to accept individual OSM ways.
"""
from contextlib import nullcontext
import hashlib
import os
import random
import tempfile
import time
from urllib.error import HTTPError
from xml.etree import ElementTree
from array import array
from pathlib import Path
from types import SimpleNamespace

import bpy
from .blender import app
from ..downloads import Request, DownloadError
from ..downloads.transport import DownloadTransport
from ..parse.osm import Osm
from ..progressive.core import Job, LazyTransport, YieldToUI
from ..progressive.ownership import Ownership, alive
from ..progressive import blender as driver
from ..renderer import Renderer

last_result = None


def result_message(result):
    """Keep original Blosm diagnostics visible in both import drivers."""
    if result.get('preliminary_failure') and result['errors']:
        return result['errors'][0]
    message = 'OSM import %s: %s completed, %s errors.' % (
        result['outcome'], result['committed'], len(result['errors']))
    if result['errors']:
        message += '\n' + result['errors'][0]
    return message


def _linked(root, collection):
    return root == collection or any(_linked(child, collection) for child in root.children)


class OsmJob(Job):
    def __init__(self, session):
        self.session = session
        super().__init__(session.steps(), LazyTransport(lambda: DownloadTransport(
            max_requests=app.downloadWorkers if app.concurrentDownloads else 1)), session.validate)

    def step(self, *args, **kwargs):
        session = self.session
        state = random.getstate()
        random.setstate(session.random_state)
        started = time.monotonic()
        session.in_step = True
        try:
            with driver.preserve_context():
                for obj in bpy.context.selected_objects:
                    obj.select_set(False)
                # This hook exists on buildings_gn. Scope it to a single turn.
                try:
                    from .import_performance import import_performance
                except ImportError:
                    scope = nullcontext()
                else:
                    scope = import_performance(bpy.context, app.addonName)
                with scope:
                    super().step(*args, **kwargs)
        finally:
            session.random_state = random.getstate()
            random.setstate(state)
            session.in_step = False
            session.max_callback = max(session.max_callback, time.monotonic() - started)


class OsmImport:
    def __init__(self, context):
        self.scene, self.view_layer = context.scene, context.view_layer
        self.settings = {p.identifier: getattr(context.scene.blosm, p.identifier)
                         for p in context.scene.blosm.bl_rna.properties
                         if p.identifier != 'rna_type' and p.type not in ('COLLECTION', 'POINTER')}
        self.requested = bool(self.settings['progressiveImport'])
        # Blender ID inputs (for example a street material override) are not
        # scalar settings, but still belong to the captured import request.
        # Do not retain RNA PropertyGroups or re-read the panel on later turns.
        self.pointer_settings = {}
        for prop in context.scene.blosm.bl_rna.properties:
            if prop.type == 'POINTER' and prop.identifier != 'rna_type':
                value = getattr(context.scene.blosm, prop.identifier)
                if value is None or isinstance(value, bpy.types.ID):
                    self.pointer_settings[prop.identifier] = value
        self.progressive = self.requested and not self.settings['commandLineMode']
        self.effective_single = self.settings['singleObject'] and not self.requested
        self.units = self.scene.unit_settings.scale_length
        self.origin_signature = (self.scene.get('lat'), self.scene.get('lon'))
        self.default_levels = tuple((item.levels, item.weight) for item in context.scene.blosm.defaultLevels)
        if not self.default_levels:
            from ..gui import _defaultLevels
            self.default_levels = tuple(_defaultLevels)
        entry = context.preferences.addons.get(app.addonName)
        defaults = dict(dataDir='', assetsDir='', osmServer=app.devOsmServer,
                        enableExperimentalFeatures=False)
        self.preferences = ({app.addonName: SimpleNamespace(preferences=SimpleNamespace(**{
            name: getattr(entry.preferences, name, default) for name, default in defaults.items()}))}
                            if entry else {})
        self.random_state = random.getstate()
        self.phase = 'Preparing import'
        self.revision = self.committed_units = self.work_done = 0
        self.errors, self.warnings = [], []
        self.rendered_nodes = set()
        self.ownership = Ownership()
        self.destination = None
        self.osm = None
        self.initialized = self.finished = self.in_step = False
        self.preflight_complete = False
        self.owns_app = False
        self.started = time.monotonic()
        self.max_callback = 0.
        self.first_result = None
        self.first_preview = None
        self.preview_count = self.preview_batches = self.peak_preview_count = 0
        self.terrain_source = None
        self.terrain_dirty = False
        self.terrain_signature = None
        self.gn_group = None
        self.origin_written = False
        self.job = OsmJob(self)
        self._terrain_handler = self._terrain_update

    @property
    def status(self):
        message = '%s: %s completed, %s errors' % (self.phase, self.committed_units, len(self.errors))
        if self.preview_count:
            message += ' (%s preview meshes)' % self.preview_count
        return message

    def preview_changed(self, count):
        """Record visible provisional output without accepting an import unit."""
        self.preview_count = count
        self.peak_preview_count = max(self.peak_preview_count, count)
        if count:
            self.preview_batches += 1
            if self.first_preview is None:
                self.first_preview = time.monotonic() - self.started
        self.revision += 1

    def report(self, levels, message):
        if 'ERROR' in levels:
            self.errors.append(str(message))
        elif 'WARNING' in levels:
            self.warnings.append(str(message))

    def _terrain_update(self, scene, graph):
        if self.terrain_source and not self.in_step:
            for update in graph.updates:
                if update.id.original in (self.terrain_source, self.terrain_source.data):
                    self.terrain_dirty = True
                    break

    def _terrain_key(self):
        obj = self.terrain_source
        coordinates = array('f', [0.]) * (len(obj.data.vertices) * 3)
        obj.data.vertices.foreach_get('co', coordinates)
        modifiers = []
        for modifier in obj.modifiers:
            values = []
            for prop in modifier.bl_rna.properties:
                if prop.identifier in ('rna_type', 'name', 'show_expanded') or prop.is_readonly:
                    continue
                value = getattr(modifier, prop.identifier, None)
                if prop.type in ('BOOLEAN', 'INT', 'FLOAT', 'STRING', 'ENUM'):
                    values.append((prop.identifier, tuple(value) if getattr(prop, 'is_array', False) else value))
                elif prop.type == 'POINTER' and isinstance(value, bpy.types.ID):
                    values.append((prop.identifier, value.as_pointer()))
            inputs = tuple((key, str(value)) for key, value in modifier.items()) if modifier.type == 'NODES' else ()
            modifiers.append((modifier.type, tuple(values), inputs))
        return (obj.data.as_pointer(), tuple(v for row in obj.matrix_world for v in row),
                hashlib.sha256(coordinates.tobytes()).digest(), len(obj.data.polygons), tuple(modifiers))

    def validate(self):
        for name, value in self.pointer_settings.items():
            if value is not None and not alive(value):
                raise RuntimeError('An import input was deleted: ' + name)
        if not alive(self.scene):
            raise RuntimeError('The import scene was deleted')
        if self.scene.unit_settings.scale_length != self.units:
            raise RuntimeError('Scene units changed during the import')
        if (self.scene.get('lat'), self.scene.get('lon')) != self.origin_signature:
            raise RuntimeError('Scene geographic origin changed during the import')
        if self.gn_group is not None and not alive(self.gn_group):
            raise RuntimeError('The import Geometry Nodes group was deleted')
        if self.destination is not None:
            if not alive(self.destination) or not _linked(self.scene.collection, self.destination):
                raise RuntimeError('The import collection was deleted or unlinked')
        if self.terrain_source is not None:
            if not alive(self.terrain_source) or self.terrain_source.name not in self.scene.objects:
                raise RuntimeError('The import terrain was deleted or unlinked')
            transform = tuple(v for row in self.terrain_source.matrix_world for v in row)
            if transform != self.terrain_signature[1] or self.terrain_dirty:
                if self._terrain_key() != self.terrain_signature:
                    raise RuntimeError('The import terrain changed; completed objects were retained')
                self.terrain_dirty = False

    def call(self, callback):
        # No yield in allocation capture: user data cannot be created between
        # this snapshot and its matching capture. Explicit own() is available
        # for resumable providers; no global snapshot spans UI turns.
        with self.ownership.capture() as allocations:
            result = callback()
        return result, allocations

    def unit(self, callback):
        self.validate()
        _, allocations = self.call(callback)
        if not self.effective_single:
            self.call(Renderer.finishUnit)
        # Deferred texture work belongs to this same, still pending unit.
        exporters = [r.textureExporter for r in app.renderers if getattr(r, 'exportMaterials', False)]
        pending_bakes = any(exporter.bakeQueue and exporter.bakeQueue.tasks for exporter in exporters)
        staged = [(obj, obj.hide_viewport, obj.hide_render) for kind, obj in allocations
                  if kind == 'objects' and alive(obj)] if pending_bakes else []
        if pending_bakes:
            self.phase = 'Rendering export textures'
            for obj, _, _ in staged:
                obj.hide_viewport = obj.hide_render = True
        for exporter in exporters:
            if hasattr(exporter, 'finish_steps'):
                yield from exporter.finish_steps()
        for obj, viewport, render in staged:
            if alive(obj):
                obj.hide_viewport, obj.hide_render = viewport, render
        self.phase = 'Importing OSM features'
        self.validate()
        outputs = [obj for kind, obj in allocations if kind == 'objects' and alive(obj)
                   and obj.type in ('MESH', 'CURVE')]
        if outputs:
            self.committed_units += 1
            self.revision += 1
            if self.first_result is None:
                self.first_result = time.monotonic() - self.started
            self._write_origin()
            # Session initialization dependencies serve these retained results.
            self.ownership.accept()
        self.work_done += 1
        yield YieldToUI()

    def _write_origin(self):
        if not self.origin_written and self.osm and self.osm.projection and self.set_origin:
            self.scene['lat'], self.scene['lon'] = self.osm.lat, self.osm.lon
            self.origin_signature = (self.osm.lat, self.osm.lon)
            self.origin_written = True

    def node_unit(self, node, callback):
        # Several OSM managers share rNodes. A node's selected renderer is global
        # to that source record, not an instruction to repeat it per manager.
        if node in self.rendered_nodes:
            return
        yield from self.unit(callback)
        self.rendered_nodes.add(node)

    def _hook(self, target, method):
        hook = 'render_extra_steps' if method == 'renderExtra' else method + '_steps'
        producer = getattr(target, hook, None) or getattr(target, method + '_steps', None)
        # A newer backend may override the synchronous method without adapting
        # its inherited producer. Preserve the override as one atomic fallback.
        mro = type(target).__mro__
        method_owner = next((base for base in mro if method in base.__dict__), None)
        hook_owner = next((base for base in mro if hook in base.__dict__), None)
        if producer and method_owner and hook_owner and mro.index(method_owner) < mro.index(hook_owner):
            producer = None
        if producer is not None:
            yield from producer(self)
        else:
            callback = getattr(target, method, None)
            if callback:
                if method in ('render', 'renderExtra'):
                    yield from self.unit(callback)
                    return
                self.call(callback)
            yield None

    def _download(self, request, path):
        content = yield request
        partial = None
        try:
            handle = tempfile.NamedTemporaryFile(dir=Path(path).parent, prefix='.blosm-', delete=False)
            partial = handle.name
            try:
                with handle:
                    for offset in range(0, len(content), 65536):
                        handle.write(content[offset:offset+65536])
                        yield None
            finally:
                del content
            os.replace(partial, path)
            partial = None
        finally:
            if partial and os.path.isfile(partial):
                os.unlink(partial)

    def check_prerequisites(self):
        app.validateOsmDirectories(bpy.context, settings=SimpleNamespace(**self.settings),
                                   preferences=self.preferences)
        setup_path = self.settings['setupScript']
        if setup_path:
            setup_path = os.path.realpath(bpy.path.abspath(setup_path))
            if not os.path.isfile(setup_path):
                raise RuntimeError('The setup script does not exist: '+setup_path)
        if self.settings['osmSource'] != 'server':
            path = os.path.realpath(bpy.path.abspath(self.settings['osmFilepath']))
            with open(path, 'rb'):
                pass
        return setup_path

    def steps(self):
        # The operator runs this first turn before activating the progress UI.
        # Keep Blender allocation, downloading and parsing after this boundary.
        setup_path = self.check_prerequisites()
        yield YieldToUI()
        app.osm_session = self
        self.owns_app = True
        app.projection = None
        for name, value in self.settings.items():
            setattr(app, name, value)
        for name, value in self.pointer_settings.items():
            setattr(app, name, value)
        app.singleObject = self.effective_single
        app.osm_default_levels = self.default_levels
        self.call(lambda: app.initOsm(self, bpy.context, defer_download=True,
                                    settings=SimpleNamespace(**self.settings), preferences=self.preferences))
        self.initialized = True
        if app.mode is app.twoD and app.gnSetup2d != '-':
            self.gn_group = bpy.data.node_groups[app.gnSetup2d]
        self.preflight_complete = True
        yield YieldToUI()
        if app.osmSource == 'server':
            self.phase = 'Downloading OSM'
            directory = Path(app.dataDir) / app.osmDir
            counter = 1
            while True:
                path = directory / (app.osmFileName % ('' if counter == 1 else '_%s' % counter))
                if not path.exists():
                    break
                counter += 1
            app.osmFilepath = str(path)
            yield from self._download(Request(app.osmServer + app.osmUrlPath % (
                app.minLon, app.minLat, app.maxLon, app.maxLat)), path)
        self.osm = osm = Osm(app)
        if setup_path:
            from ..util.setup_script import load_setup_script
            setup = load_setup_script(setup_path)
        elif app.mode is app.realistic:
            if app.enableExperimentalFeatures:
                from ..setup.realistic_dev import setup
            else:
                from ..setup.premium_default import setup
        else:
            from ..setup.base import setup
        self.call(lambda: setup(app, osm))
        self.call(lambda: app.createLayers(osm))
        self.set_origin = not (None not in self.origin_signature and app.relativeToInitialImport)
        if not self.set_origin:
            osm.setProjection(*self.origin_signature)
        elif app.osmSource == 'server':
            osm.setProjection((app.minLat+app.maxLat)/2., (app.minLon+app.maxLon)/2.)
        flat = app.mode is app.realistic and app.forests
        extent = flat and app.osmSource == 'file'
        self.phase = 'Parsing OSM'
        yield from osm.parse_steps(app.osmFilepath, forceExtentCalculation=extent)
        if app.loadMissingMembers and app.incompleteRelations:
            self.phase = 'Resolving missing relation members'
            extra = app.osmFileExtraName % app.osmFilepath[:-4]
            if not os.path.isfile(extra):
                cache = None
                try:
                    # Bound the IPC request frame, while preserving query semantics.
                    ids = sorted(app.missingWays)
                    merged = tempfile.NamedTemporaryFile(dir=Path(extra).parent,
                                                         prefix='.blosm-members-', delete=False)
                    cache = merged.name
                    with merged:
                        merged.write(b'<osm version="0.6">')
                        for start in range(0, len(ids), 1000):
                            query = app.overpassWays % ');way('.join(ids[start:start+1000])
                            part = cache + '.part'
                            try:
                                yield from self._download(Request(app.osmServer + app.osmUrlPath2,
                                    method='POST', body=query, retry_safe=True), part)
                                app.loadMissingMembers = False
                                yield from osm.parse_steps(part)
                                for record in Osm._elements(part):
                                    if record is not None:
                                        merged.write(ElementTree.tostring(record))
                                    else:
                                        yield None
                            finally:
                                if os.path.exists(part):
                                    os.unlink(part)
                        merged.write(b'</osm>')
                    os.replace(cache, extra)
                    cache = None
                    app.loadMissingMembers = False
                except (HTTPError, DownloadError) as exc:
                    self.errors.append('Missing relation members: %s' % exc)
                    app.loadMissingMembers = False
                finally:
                    if cache and os.path.isfile(cache):
                        os.unlink(cache)
            else:
                app.loadMissingMembers = False
                yield from osm.parse_steps(extra)
            # Repair uses the same callbacks/order as the synchronous importer.
            for relation, identity, members, tags, condition in app.incompleteRelations:
                relation.process(members, tags, osm)
                if relation.valid:
                    skip = osm.processCondition(condition, relation, identity, osm.parseRelation)
                    if identity not in osm.relations and not skip:
                        osm.relations[identity] = relation
                yield None
            app.incompleteRelations = app.missingWays = None
            if not osm.projection:
                osm.setProjection((osm.minLat+osm.maxLat)/2., (osm.minLon+osm.maxLon)/2.)
        if extent:
            app.minLat, app.maxLat, app.minLon, app.maxLon = osm.minLat, osm.maxLat, osm.minLon, osm.maxLon
            if osm.rNodes and not osm.ways:
                # Point-only extracts need a small border: a ray exactly on a
                # float32 terrain boundary can otherwise miss an extreme point.
                app.minLat -= .00001
                app.maxLat += .00001
                app.minLon -= .00001
                app.maxLon += .00001
        self.phase = 'Preparing terrain and layers'
        self.call(lambda: app.setTerrain(bpy.context, createFlatTerrain=flat, createBvhTree=True,
                                        settings=SimpleNamespace(**self.settings)))
        if app.terrain:
            self.terrain_source = app.terrain.terrain
            self.terrain_signature = self._terrain_key()
            bpy.app.handlers.depsgraph_update_post.append(self._terrain_handler)
        self.call(app.initLayers)
        if self.requested and any(layer.singleObject for layer in app.layers):
            # AreaManager's progressive adapter uses separate area geometry.
            raise RuntimeError('A setup enabled single-object layers during progressive import')
        self.phase = 'Preparing OSM features'
        # The street/facade managers on buildings_gn consume projected bounds.
        if osm.projection:
            app.minX, app.minY = app.projection.fromGeographic(osm.minLat, osm.minLon)[:2]
            app.maxX, app.maxY = app.projection.fromGeographic(osm.maxLat, osm.maxLon)[:2]
        if app.logger:
            app.logger.processStart()
        for manager in app.managers:
            yield from self._hook(manager, 'process')
        if app.logger:
            app.logger.processEnd()
            app.logger.renderStart()
        self.phase = 'Importing OSM features'
        # Establish the destination before prepare(), matching newer building
        # and street backends that allocate their setup resources there.
        self.call(lambda: Renderer.begin(app))
        self.destination = Renderer.collection
        for renderer in app.renderers:
            self.call(renderer.prepare)
            yield None
        yield YieldToUI()
        for manager in app.managers:
            yield from self._hook(manager, 'render')
        self.phase = 'Finalizing OSM import'
        self.call(lambda: Renderer.end(app))
        for manager in app.managers:
            yield from self._hook(manager, 'renderExtra')
        if self.effective_single:
            for layer in app.layers:
                finalize = getattr(layer, 'finalize', None)
                if finalize:
                    self.call(finalize)
                    yield None
        # Optional hooks preserve the newer branch's lifecycle without importing
        # any building or street backend into the generic session.
        for renderer in app.renderers:
            finalize = getattr(renderer, 'finalize', None)
            if finalize:
                self.call(finalize)
                yield None
        if app.logger:
            app.logger.renderEnd()
        self._write_origin()
        self.ownership.accept()
        self.revision += 1

    def cancel(self):
        self.job.cancel()

    def finish(self, silent=False):
        global last_result
        if self.finished:
            return
        self.finished = True
        if self._terrain_handler in bpy.app.handlers.depsgraph_update_post:
            bpy.app.handlers.depsgraph_update_post.remove(self._terrain_handler)
        if self.job.error:
            self.errors.append(str(self.job.error))
        # Clean initialized resources even if setup/start failed before yielding.
        cleanup_errors = []
        with driver.preserve_context():
            for target in (list(app.renderers or []) + list(app.managers or [])) if self.owns_app else []:
                cleanup = getattr(target, 'cleanup', None)
                if cleanup:
                    try:
                        cleanup()
                    except Exception as exc:
                        cleanup_errors.append(str(exc))
            if self.owns_app:
                for target in list(app.renderers or []) + list(app.layers or []):
                    bm = getattr(target, 'bm', None)
                    if bm is not None and bm.is_valid:
                        bm.free()
                        target.bm = None
                terrain = getattr(app, 'terrain', None)
                if terrain and terrain.terrain:
                    terrain.cleanup()
            if self.ownership.pending:
                self.ownership.rollback()
                self.revision += 1
            if self.owns_app:
                app.clean()
        self.errors.extend(cleanup_errors)
        if getattr(app, 'osm_session', None) is self:
            app.osm_session = None
        outcome = ('cancelled' if self.job.cancelled else
                   'partial' if self.errors and self.committed_units else
                   'failed' if self.errors else 'success')
        last_result = dict(outcome=outcome, committed=self.committed_units, errors=list(self.errors),
                           preliminary_failure=outcome == 'failed' and not self.preflight_complete,
                           warnings=list(self.warnings), first_result_seconds=self.first_result,
                           first_preview_seconds=self.first_preview, preview_batches=self.preview_batches,
                           peak_preview_meshes=self.peak_preview_count,
                           max_callback_seconds=self.max_callback, elapsed_seconds=time.monotonic()-self.started)
        for error in self.errors:
            print('OSM import:', error)
        last_result['message'] = result_message(last_result)
        if not silent:
            driver.notify_result('ERROR' if outcome == 'failed' else 'INFO', last_result['message'])
        return last_result


def import_osm(operator, context):
    requested = context.scene.blosm.progressiveImport and not context.scene.blosm.commandLineMode
    if requested and (bpy.app.background or context.window is None):
        operator.report({'ERROR'}, 'Progressive import requires an interactive Blender window. Disable Progressive import for background jobs.')
        return {'CANCELLED'}
    if requested and context.mode != 'OBJECT':
        operator.report({'ERROR'}, 'Start progressive import in Object Mode')
        return {'CANCELLED'}
    if not requested and context.mode != 'OBJECT':
        operator.setObjectMode(context)
    session = OsmImport(context)
    if requested:
        session.job.step(max_steps=1)
        if session.job.done:
            result = session.finish(silent=True)
            operator.report({'ERROR'}, result['message'])
            return {'CANCELLED'}
        try:
            driver.start(session, context)
        except Exception as exc:
            session.cancel()
            session.finish(silent=True)
            operator.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        return {'FINISHED'}
    try:
        while not session.job.done:
            session.job.step()
            if session.job.waiting:
                time.sleep(.005)
    finally:
        result = session.finish(silent=True)
    operator.report({'ERROR'} if result['outcome'] == 'failed' else {'INFO'},
                    result['message'])
    return {'CANCELLED'} if result['outcome'] == 'failed' else {'FINISHED'}
