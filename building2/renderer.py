import bpy

from .. import defs
from ..renderer import Renderer
from .layer import BuildingLayer
from .item_store import ItemStore
from .texture_exporter import TextureExporter

from ..item.building import Building
from ..item.footprint import Footprint
from ..item.facade import Facade
from ..item.level import Level
from ..item.div import Div
from ..item.bottom import Bottom
from ..item.window import Window
from ..item.entrance import Entrance
from ..item.corner import Corner
from ..item.balcony import Balcony
from ..item.chimney import Chimney

from ..item.roof_flat import RoofFlat
from ..item.roof_flat_multi import RoofFlatMulti
from ..item.roof_profile import RoofProfile
from ..item.roof_generatrix import RoofGeneratrix
from ..item.roof_hipped import RoofHipped
from ..item.roof_hipped_multi import RoofHippedMulti
from ..item.roof_side import RoofSide

from .content_generator import MeshGenByIndices, MeshGenByCoords

from ..util.osm import assignTags


_itemClasses = (
        Building,
        Footprint,
        Facade,
        Level,
        Div,
        Bottom,
        Window,
        Entrance,
        Corner,
        Balcony,
        Chimney,
        RoofFlat,
        RoofFlatMulti,
        RoofProfile,
        RoofGeneratrix,
        RoofHipped,
        RoofHippedMulti,
        RoofSide
    )


class BuildingRendererNew(Renderer):
    
    def __init__(self, app, styleStore, itemRenderers, getStyle=None):
        self.app = app
        app.addRenderer(self)
        
        self.styleStore = styleStore
        
        self.assetsDir = app.assetsDir
        self.assetPackageDir = app.assetPackageDir
        
        # do wee need to apply a cladding color for facade textures?
        self.useCladdingColor = True
        
        self.itemRenderers = itemRenderers
        self.footprintRenderer = itemRenderers["Footprint"]
        self.facadeRenderer = itemRenderers["Facade"]
        
        self.exportMaterials = app.enableExperimentalFeatures and app.importForExport
        
        if self.exportMaterials:
            self.textureExporter = TextureExporter(self.assetsDir, self.assetPackageDir)
            # Do we need to cache <claddingTextureInfo> for each cladding material?
            self.cacheCladdingTextureInfo = False
        
        # initialize item renderers
        for itemRenderer in itemRenderers.values():
            itemRenderer.init(itemRenderers, self)
        
        self.getStyle = getStyle
        self.itemStore = ItemStore(_itemClasses)
        
        self._cache = {}
        
        self.actions = []
        # "rev" stands for "render extruded volumes"
        self.revActions = []
    
    def prepare(self):
        self.meshAssets = {} # FIXME: this line is defined in the commented code below
        """ FIXME: uncomment the code below or remove it if you don't need Geometry Nodes setup for buildings
        if self.app.preferMesh:
            # A key to the dictionary below is an object name as it's defined in <self.app.assetStore>.
            # A related value is tuple that consists of 
            # (1) A link to a Blender object in the file defined in the related asset info.
            # (2) A list of custom attributes defined for this object. The addon will
            # create a separate object for each unique set of the attributes. All those created objects will
            # share the same Blender mesh.
            # (3) A dictionary:
            #     A key is a combination of custom attributes described above.
            #     A value is the name of the Blender object created for the unique set of attributes
            #     used in the key.
            self.meshAssets = {}
            # A Python dictionary to store processed information about Blender objects that belong to
            # a Blender collection. Each collection is defined in the asset store.
            self.blenderCollections = {}
            
            _collectionName = "blosm_building_assets"
            if not _collectionName in bpy.data.collections:
                bpy.data.collections.new(_collectionName)
            # a Blender collection for instances on the points cloud
            self.buildingAssetsCollection = bpy.data.collections[_collectionName]
            
            
            # check if the Geometry Nodes setup with the name "blosm_gn_building" is already available
            _gnName, _gnCollection = "blosm_gn_building", "Collection Info"
            node_groups = bpy.data.node_groups
            self.gnBuilding = node_groups[_gnName]\
                if _gnName in node_groups and _gnCollection in node_groups[_gnName].nodes else\
                appendNodeGroupFromFile(self.app.baseAssetPath, _gnName)
            # set the input of <self.gnBuilding.nodes[_gnCollection]> to <self.buildingAssetsCollection>
            self.gnBuilding.nodes[_gnCollection].inputs['Collection'].default_value = self.buildingAssetsCollection
            
            # a TEMPORARY code below to load the Geometry Nodes setup for flat roofs
            self.gnFlatRoof = appendNodeGroupFromFile(
                os.path.join(os.path.dirname(self.app.baseAssetPath), "flat_roof_objects.blend"),
                "blosm_flat_roof_objects"
            )
        """
        
        if self.app.singleObject:
            self.genVolumes = MeshGenByIndices(
                self.createBlenderObject(
                    "Building Volumes",
                    (0., 0., 0.),
                    collection = self.collection,
                    parent = None
                )
            )
            if self.app.preferableResult == defs.Result.FootprintWithGn:
                self.genHoles = MeshGenByIndices(
                    self.createBlenderObject(
                        "Multipolygon Holes",
                        (0., 0., 0.),
                        collection = self.collection,
                        parent = None
                    )
                )

                self.genPartFootprints = MeshGenByIndices(
                    self.createBlenderObject(
                        "Building Part Footprints",
                        (0., 0., 0.),
                        collection = self.collection,
                        parent = None
                    )
                )
                self.genPartFootprints.initFootprintAttributes()
            
            for layer in self.app.layers:
                if isinstance(layer, BuildingLayer):
                    # <self.genVolumes> must be accessible from <layer>
                    layer.genVolumes = self.genVolumes
                    if self.app.preferableResult == defs.Result.FootprintWithGn:
                        layer.genFootprints = MeshGenByIndices(
                            self.createBlenderObject(
                                "Building Footprints",
                                layer.location,
                                collection = self.collection,
                                parent = None
                            )
                        )
                        layer.genFootprints.initFootprintAttributes()
                        # <self.genHoles> must be accessible from <layer>
                        layer.genHoles = self.genHoles
                        # <self.genPartFootprints> must be accessible from <layer>
                        layer.genPartFootprints = self.genPartFootprints
                    layer.prepare()
    
    def finalize(self):
        if self.app.singleObject:
            for layer in self.app.layers:
                if isinstance(layer, BuildingLayer):
                    layer.finalize(self)
            
            self.genVolumes.finalize()
            if self.app.preferableResult == defs.Result.FootprintWithGn:
                self.genPartFootprints.finalize()
                self.genHoles.finalize()
    
    def cleanup(self):
        for action in self.actions:
            action.cleanup()
        
        if self.exportMaterials:
            self.textureExporter.cleanup()
        
        self._cache.clear()
    
    def render(self, building, data):
        parts = building.parts
        itemStore = self.itemStore
        
        #if "id" in outline.tags: print(outline.tags["id"]) #DEBUG OSM id
        
        building.renderInfo = Building(data)
        
        # get the style of the building
        buildingStyle = self.styleStore.get(self.getStyle(building, self.app))
        if not buildingStyle:
            # skip the building
            return
        building.renderInfo.setStyleMeta(buildingStyle)
        
        if self.app.preferableResult == defs.Result.FootprintWithGn or not parts or building.alsoPart:
            # the building has no parts
            footprint = Footprint(building, building)
            # The attribute <footprint> below may be used:
            # * to generate a footprint of the whole building for subsequent
            #   builing generation with Geometry Nodes
            # * to calculate the area of the building footprint
            # * in <action.terrain.Terrain>
            building.footprint = footprint
            itemStore.add(footprint)
        if parts:
            itemStore.add((Footprint(part, building) for part in parts), Footprint, len(parts))
        
        if not self.app.renderAfterExtrude and not self.app.singleObject:
            self.preRender(building)

        for action in self.actions:
            action.do(building, buildingStyle, self)
            if itemStore.skip:
                # <building.polygon> equal to <None> means that <building> was skipped.
                # It can be used later in the code if two building footprints share an edge
                building.polygon = None
                break
        itemStore.clear()
        
        if itemStore.skip:
            itemStore.skip = False
        elif self.app.renderAfterExtrude:
            self.postRender(building)
    
    def renderExtrudedVolumes(self, building, data):
        if not self.app.singleObject:
            l = building.element.l
            gen = building.renderInfo.gen
            gen.init(
                self.createBlenderObject(
                    self.getName(building.element),
                    building.renderInfo.offsetBlenderObject,
                    collection = l.getCollection(self.collection),
                    parent = l.getParent(l.getCollection(self.collection))
                )
            )
            gen.initFootprintAttributes()
            if self.app.preferableResult == defs.Result.FootprintWithGn:
                l.genVolumes = l.genFootprints = l.genPartFootprints = l.genHoles = gen
            else:
                l.genVolumes = gen
            l.prepare()
        
        # render building footprint
        if self.app.preferableResult == defs.Result.FootprintWithGn:
            self.footprintRenderer.render(building.footprint, data)
        if (not building.parts or building.alsoPart) and not building.footprint.doFootprintOnly:
            self.renderExtrudedVolume(building.footprint)
        # render building parts
        for part in building.parts:
            if self.app.preferableResult == defs.Result.FootprintWithGn:
                self.footprintRenderer.renderPart(part.footprint)
            if not part.footprint.doFootprintOnly:
                self.renderExtrudedVolume(part.footprint)
        
        if not self.app.singleObject:
            self.postRender(building)
    
    def renderExtrudedVolume(self, footprint):
        if not footprint:
            return
        
        if not footprint.noWalls:
            for action in self.revActions:
                action.do(footprint)
            
            self.facadeRenderer.render(footprint)
        
        footprint.roofRenderer.render(footprint.roofItem)
    
    def preRender(self, building):
        layer = building.element.l

        gen = building.renderInfo.gen = MeshGenByIndices()
        if self.app.preferableResult == defs.Result.FootprintWithGn:
            layer.genVolumes = layer.genFootprints = layer.genPartFootprints = layer.genHoles = gen
        else:
            layer.genVolumes = gen
    
    def postRender(self, building):
        # assign OSM tags to the blender object
        assignTags(building.renderInfo.gen.obj, building.element.tags)
        building.element.l.finalize()
        building.renderInfo.gen.finalize()
    
    def setUvs(self, face, uvs, layer, uvLayerName):
        # assign uv coordinates
        uvLayer = layer.genVolumes.bm.loops.layers.uv[uvLayerName]
        loops = face.loops
        for loop,uv in zip(loops, uvs):
            loop[uvLayer].uv = uv
    
    def setVertexColor(self, face, color, layer, layerName):
        vertexColorLayer = layer.genVolumes.bm.loops.layers.color[layerName]
        for loop in face.loops:
            loop[vertexColorLayer] = color

    def setMaterial(self, layer, face, materialName):
        """
        Set material (actually material index) for the given <face>.
        """
        materialIndices = layer.materialIndices
        materials = layer.genVolumes.obj.data.materials
        
        if not materialName in materialIndices:
            materialIndices[materialName] = len(materials)
            materials.append(bpy.data.materials[materialName] if materialName else None)
        face.material_index = materialIndices[materialName]