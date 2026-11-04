"""Simple 3D import manager; future modes can select another PADL program role."""
import json

from ...manager import Manager
from ...parse.osm.relation.multipolygon import Multipolygon
from ...parse.osm.relation.building import Building as BuildingRelation
from ...renderer import Renderer
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
        import bpy
        footprints = OsmFootprintSource(self.osm, self.mapper, self.report).collect(self.ways, self.relations)
        results = []
        if footprints:
            terrain = self.app.terrain.terrain if self.app.terrain else None
            groups = group_buildings(footprints, self.app.singleObject)
            with self.service.batch():
                for group in groups:
                    result = self.service.generate(group, collection=Renderer.collection, terrain=terrain,
                        seed=getattr(self.app, 'buildingSeed', 0),
                        name='Buildings simple 3D' if self.app.singleObject else group[0].identity)
                    results.append(result)
                    self.report.statuses.extend(result.statuses)
        text = bpy.data.texts.new('PADL building import report')
        text.write(json.dumps(self.report.as_dict(), indent=2) + '\n')
        for result in results:
            result.owner['blosm_building_report'] = text
        self.app.padlBuildingResults = tuple(results)
        self.app.padlBuildingResult = results[0] if len(results) == 1 else None
        self.app.padlBuildingReport = self.report.as_dict()
        failures = sum(bool(row.get('error')) for row in self.report.statuses)
        warnings = len(self.report.issues) + failures
        if warnings:
            self.app.op.report({'WARNING'}, f'PADL: {self.report.input_count} footprints, {failures} generation failures, {len(self.report.issues)} mapping/source notices. See {text.name}')

    def cleanup(self):
        self.ways.clear()
        self.relations.clear()
        # Legacy terrain initialization may temporarily move its origin. Simple
        # PADL imports have no legacy building renderer to restore that origin.
        terrain = self.app.terrain
        if terrain and terrain.terrain is not None:
            terrain.cleanup()
