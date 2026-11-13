"""Simple 3D import manager; future modes can select another PADL program role."""
import json
from contextlib import ExitStack

from ...manager import Manager
from ...parse.osm.relation.multipolygon import Multipolygon
from ...parse.osm.relation.building import Building as BuildingRelation
from ...renderer import Renderer
from ...progressive.ownership import alive
from .assets import BuildingProgramLibrary
from .backends import BuildingGenerationService
from .mapping import OsmFaceAttributeMapper
from .source import ImportReport, OsmFootprintSource, group_buildings, tag_enabled


class PadlBuildingManager(Manager):
    def __init__(self, osm, app, *, mapper=None, service=None):
        super().__init__(osm)
        self.id, self.layerClass, self.app = 'padl_buildings', False, app
        self.ways, self.relations = {}, {}
        self.report = ImportReport()
        self.results, self.report_text = [], None
        self._publish_results()
        self.mapper = mapper or OsmFaceAttributeMapper(projection=osm.projection or app.projection)
        self.service = service or BuildingGenerationService(
            library=BuildingProgramLibrary(app.assetPath, app=app),
            backend=getattr(app, 'buildingBackend', 'geometry_nodes'))
        # Validate the document before parsing an entire extract.
        self.service.library.program(self.service.role)
        app.addManager(self)

    @staticmethod
    def accepts(tags, element):
        return isinstance(element, BuildingRelation) or tag_enabled(tags, 'building') or tag_enabled(tags, 'building:part')

    def parseWay(self, element, elementId):
        source = 'osm:way:' + str(elementId)
        if tag_enabled(element.tags, 'building'):
            element.b = self  # Allows the parser's building-relation ownership pass.
        if element.closed:
            self.ways[source] = element
        else:
            self.report.issue(source, 'invalid_geometry', 'Building way is not closed')

    def parseRelation(self, element, elementId):
        source = 'osm:relation:' + str(elementId)
        if isinstance(element, BuildingRelation):
            return
        if isinstance(element, Multipolygon):
            if tag_enabled(element.tags, 'building'):
                element.b = self
            self.relations[source] = element
        else:
            self.report.issue(source, 'unsupported_relation', 'Only building multipolygons are supported in this profile')

    def render(self):
        for _ in self.render_steps():
            pass

    def render_steps(self, session=None):
        """Publish a complete outline/parts group, including editable inputs.

        Relationship resolution precedes generation so parts never appear as
        unrelated buildings. Both drivers use the same backend and grouping.
        """
        if session is not None and not session.requested:
            # Blocking imports keep one batch/reader and one ownership scope.
            yield from self._unit(session, self.render)
            return
        footprints = OsmFootprintSource(self.osm, self.mapper, self.report).collect(self.ways, self.relations)
        terrain = self.app.terrain.terrain if self.app.terrain else None
        groups = group_buildings(footprints, self.app.singleObject)
        # Reuse the reader for the batch. Progressive status helpers live in a
        # private scene; transient ownership removes them on teardown/history.
        # Each read releases its borrowed building mesh before the next yield.
        progressive = session is not None
        with ExitStack() as batch:
            if groups:
                def begin():
                    batch.enter_context(self.service.batch(progressive=progressive))
                if progressive:
                    session.call(begin)
                else:
                    begin()
            for group in groups:
                def generate():
                    result = self.service.generate(group, collection=Renderer.collection, terrain=terrain,
                        seed=getattr(self.app, 'buildingSeed', 0),
                        name='Buildings simple 3D' if self.app.singleObject else group[0].identity)
                    self.results.append(result)
                    self.report.statuses.extend(result.statuses)
                    if progressive:
                        self._write_report()

                if session is None:
                    generate()
                    yield None
                else:
                    yield from self._unit(session, generate)
        if not progressive or not groups:
            if session is None:
                self._write_report()
            else:
                session.call(self._write_report)
        self._report_warnings()

    def _unit(self, session, callback):
        """Publish metadata with the unit, restoring it if acceptance fails."""
        committed = session.committed_units
        result_count, status_count = len(self.results), len(self.report.statuses)
        text = self.report_text
        content = text.as_string() if text is not None and alive(text) else None
        bindings = [(result.owner, result.owner.get('blosm_building_report'))
                    for result in self.results if alive(result.owner)]
        try:
            yield from session.unit(callback)
        except BaseException:
            # Closing a generator after its accepted UI boundary must retain
            # that building. Before acceptance, the session reclaims new IDs;
            # restore mutations to the already accepted report separately.
            if session.committed_units == committed:
                del self.results[result_count:]
                del self.report.statuses[status_count:]
                if content is not None and alive(text):
                    text.clear()
                    text.write(content)
                self.report_text = text
                for owner, previous in bindings:
                    if alive(owner):
                        if previous is not None and alive(previous):
                            owner['blosm_building_report'] = previous
                        elif 'blosm_building_report' in owner:
                            del owner['blosm_building_report']
                self._publish_results()
            raise

    def _publish_results(self):
        # Previously accepted objects may be renamed or deleted between turns.
        results = tuple(result for result in self.results
                        if alive(result.owner) and alive(result.output))
        self.app.padlBuildingResults = results
        self.app.padlBuildingResult = results[0] if len(results) == 1 else None
        self.app.padlBuildingReport = self.report.as_dict()

    def _write_report(self):
        import bpy
        if self.report_text is None or not alive(self.report_text):
            self.report_text = bpy.data.texts.new('PADL building import report')
        text = self.report_text
        text.clear()
        text.write(json.dumps(self.report.as_dict(), indent=2) + '\n')
        for result in self.results:
            # Assigning even the same ID property tags an accepted GN object
            # for another update. Attach new/recreated reports without touching
            # all previous owners on every progressive turn.
            if alive(result.owner) and result.owner.get('blosm_building_report') != text:
                result.owner['blosm_building_report'] = text
        self._publish_results()

    def _report_warnings(self):
        failures = sum(bool(row.get('error')) for row in self.report.statuses)
        warnings = len(self.report.issues) + failures
        if warnings:
            message = f'PADL: {self.report.input_count} footprints, {failures} generation failures, {len(self.report.issues)} mapping/source notices.'
            if self.report_text is not None and alive(self.report_text):
                message += f' See {self.report_text.name}'
            self.app.op.report({'WARNING'}, message)

    def cleanup(self):
        self.ways.clear()
        self.relations.clear()
        self._publish_results()
        # Legacy terrain initialization may temporarily move its origin. Simple
        # PADL imports have no legacy building renderer to restore that origin.
        terrain = self.app.terrain
        if terrain and terrain.terrain is not None:
            terrain.cleanup()
