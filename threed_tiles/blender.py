from os.path import basename, join as joinStrings, exists as pathExists
from os import remove as removeFile
from math import radians, atan, sqrt, pi
import re
from operator import itemgetter
import sys

import bpy
from mathutils import Vector, Matrix

from util.blender import createCollection

from io_scene_gltf2.io.imp.gltf2_io_gltf import glTFImporter


def find_layer_collection(layer_col: bpy.types.LayerCollection, target_col: bpy.types.Collection):
    """
    Recursively searches for the LayerCollection that corresponds to the given <target_col> of type <bpy.types.Collection>
    within the hierarchy starting from <layer_col> of type <bpy.types.LayerCollection>
    Returns the found LayerCollection or None if not found.
    """
    if layer_col.collection == target_col:
        return layer_col
    for child in layer_col.children:
        found = find_layer_collection(child, target_col)
        if found:
            return found
    return None


class BlenderRenderer:
    
    def __init__(self, threedTilesName, join3dTilesObjects, instanceName):
        self.threedTilesName = threedTilesName
        self.join3dTilesObjects = join3dTilesObjects
        self.instanceName = instanceName
        
        self.calculateHeightOffset = False
        self.heightOffset = 0.
        
        self.licenseRePattern = re.compile(b'\"copyright\":\s*\"(\w+)\"')
        self.copyrightHolders = {}
        
        # <self._gltfImporterPatched> is used to store some Python objects related to Blender's glTF importer if its patching is performed
        self._gltfImporterPatched = None
        # the original static function <BlenderGlTF.set_convert_functions> of the glTF importer will be stored in the attribute below
        self._set_convert_functions = None
        # the original static function <BlenderScene.select_imported_objects> of the glTF importer will be stored in the attribute below
        self._select_imported_objects = None
    
    def prepare(self, manager):
        self.num_imported_tiles = 0

        self.collection = createCollection(self.threedTilesName)
        # <self.collection_import> is used to import the objects first, then they are moved to <self.collection>.
        # The reason for this is to isolatede the imported objects from other objects in the scene during the import process.
        self.collection_import = createCollection(self.threedTilesName + "_import")
        # LayerCollection corresponding to <self.collection_import>
        self.layer_collection_import = find_layer_collection(
            bpy.context.view_layer.layer_collection,
            self.collection_import
        )
        # set the active LayerCollection to <self.layer_collection_import>
        bpy.context.view_layer.active_layer_collection = self.layer_collection_import
        
        self.centerCoords = manager.fromGeographic(manager.centerLat, manager.centerLon, 0.)
        
        # Blender's glTF importer will be patched by default. A position offset will be applied during the import.
        glTFImporter._patch_convert_functions = True
        # the check <if "io_scene_gltf2" in sys.modules> is performed in <blosm.app>
        self.patchGltfImporter()
    
    def finalize(self, manager):
        bpy.data.collections.remove(self.collection_import)
        self.collection_import = None
        self.layer_collection_import = None

        if not self.collection.objects:
            bpy.data.collections.remove(self.collection)
            self.collection = None

            if self._gltfImporterPatched:
                self.cleanupGltfImporterPatching()
            return
        
        bpy.ops.object.select_all(action='DESELECT')

        # check if there are any parent-child relationships among the imported objects
        if sum(1 for obj in self.collection.objects if obj.parent):
            # Out task is to remove all parent-child relationships among the imported objects,
            # i.e. flatten the hierarchy.

            # Build a dictionary Object Name --> Object for the imported Mesh objects
            mesh_objects = {obj.name: obj for obj in self.collection.objects if obj.type == 'MESH'}
            while True:
                num_processed_objects = 0
                for obj in list(mesh_objects.values()):
                    # Only bottom-level Mesh objects with a parent (i.e at the very bottom of a parent-child hierarchy)
                    # are processed in each iteration
                    if obj.parent and not obj.children:
                        # get final world transform (includes all parents)
                        mw = obj.matrix_world.copy()
                        # clear parent
                        obj.parent = None
                        # apply the final world transform to the object
                        obj.matrix_world = mw
                        # remove <obj> from the dictionary <mesh_objects>
                        del mesh_objects[obj.name]
                        num_processed_objects += 1
                if num_processed_objects == 0 or not mesh_objects:
                    # no more Mesh objects with a parent and without children are left
                    break
            # remove imported Blender Empty objects
            for obj in [obj for obj in self.collection.objects if obj.type == 'EMPTY']:
                bpy.data.objects.remove(obj)
        
        # select all imported Mesh objects:
        for obj in self.collection.objects:
            obj.select_set(True)
        
        # apply possible rotation
        bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)

        #
        # tranformation matrix
        #
        centerCoords = self.centerCoords
        # lat = radians(manager.centerLat - 90.) # gives incorrect result for the expression below
        lat = atan(centerCoords[2]/sqrt(centerCoords[0]*centerCoords[0] + centerCoords[1]*centerCoords[1]))
        # Rotate the mesh, so it will point to the north pole. The rotations are around Z and X axes
        matrix = Matrix.Rotation(lat-pi/2., 4, 'X') @ Matrix.Rotation(radians(-90. - manager.centerLon), 4, 'Z')
        
        locationsAfterRotation = [(matrix @ obj.location) for obj in self.collection.objects]
        
        # find the lowest Z-coordinate if <self.calculateHeightOffset>
        heightOffset = min(location[2] for location in locationsAfterRotation)\
            if self.calculateHeightOffset else\
            self.heightOffset
        if self.calculateHeightOffset:
            self.heightOffset = heightOffset
        

        if self.join3dTilesObjects:
            self.joinObjects()
            
            # set the origin of the resulting Blender object at <centerCoords>
            _cursorLocation = bpy.context.scene.cursor.location.copy()
            bpy.context.scene.cursor.location = (0., 0., 0.) if self._gltfImporterPatched else centerCoords
            bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
            bpy.context.scene.cursor.location = _cursorLocation
            
            #location = locationsAfterRotation[-1]
            #location[2] -= heightOffset
            # set the matrix_local of the resulting Blender object
            self.collection.objects[0].matrix_local = matrix#Matrix.Translation(location) @ matrix
        else:
            if not self._gltfImporterPatched:
                # rotate the vector <centerCoords>
                centerCoords = matrix @ centerCoords
            for obj, location in zip(self.collection.objects, locationsAfterRotation):
                if not self._gltfImporterPatched:
                    location[2] -= centerCoords[2]
                obj.matrix_local = Matrix.Translation(location) @ matrix
                obj.select_set(True)
        
        bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)
        bpy.ops.object.select_all(action='DESELECT')
        
        bpy.context.scene.blosm.copyright = "; ".join(
            entry[0] for entry in reversed(
                sorted(
                    self.copyrightHolders.items(), key=itemgetter(1)
                )
            )
        )

        if self._gltfImporterPatched:
            self.cleanupGltfImporterPatching()
        
        self.collection = None
    
    def renderGlb(self, manager, uri, path, transformMatrix, cacheContent):
        context = bpy.context
        
        filepath = joinStrings(
            manager.tilesDir,
            basename(path) if cacheContent else ("current_file_" + self.instanceName + ".glb" if self.instanceName else "current_file.glb")
        )
        
        if cacheContent:
            if not pathExists(filepath):
                fileContent = manager.download(uri)
                with open(filepath, 'wb') as f:
                    f.write(fileContent)
            bpy.ops.import_scene.gltf(filepath=filepath, import_scene_as_collection=True)
        else:
            fileContent = manager.download(uri)

            # check if <fileContent> contains copyright information
            match = re.search(self.licenseRePattern, fileContent)
            if match:
                self.processCopyrightInfo(match.group(1).decode('utf-8'))
            
            with open(filepath, 'wb') as f:
                f.write(fileContent)
            
            bpy.ops.import_scene.gltf(filepath=filepath, import_scene_as_collection=True)

            if self.collection_import.objects:
                self.finalize_glb_import(filepath, transformMatrix)

    def renderB3dm(self, manager, uri, path, transformMatrix, cacheContent):
        import numpy
        from .py3dtiles.tileset.content.tile_content_reader import read_array
        
        filepath = joinStrings(
            manager.tilesDir,
            basename(path)[:-4] + "glb" if cacheContent else ("current_file_" + self.instanceName + ".glb" if self.instanceName else "current_file.glb")
        )
        
        if cacheContent:
            if not pathExists(filepath):
                fileContent = manager.download(uri)
                with open(filepath, 'wb') as f:
                    f.write(fileContent)
            bpy.ops.import_scene.gltf(filepath=filepath, import_scene_as_collection=True)
        else:
            fileContent = manager.download(uri)
            # check if <fileContent> contains copyright information
            match = re.search(self.licenseRePattern, fileContent)
            if match:
                self.processCopyrightInfo(match.group(1).decode('utf-8'))
            
            b3dmContent = read_array( numpy.frombuffer(fileContent, dtype=numpy.uint8) )
            if b3dmContent is None or b3dmContent.header is None:
                raise Exception("The file doesn't contain a valid data.")
            
            gltfContent = b3dmContent.body.gltf
            # RTC means "Relative To Center"
            rtc_center = b3dmContent.body.feature_table.header.data.get("RTC_CENTER")
            if rtc_center:
                from io_scene_gltf2.io.imp.gltf2_io_gltf import glTFImporter
                # No need to patch Blender's glTF importer if the property <RTC_CENTER> is provided.
                glTFImporter._patch_convert_functions = False
            
            with open(filepath, 'wb') as f:
                f.write(gltfContent.to_array())
            
            bpy.ops.import_scene.gltf(filepath=filepath, import_scene_as_collection=True)
            
            if self.collection_import.objects or self.collection_import.children:
                if rtc_center:
                    self.process_rtc(rtc_center)
                self.finalize_glb_import(filepath, transformMatrix)
    
    def processCopyrightInfo(self, info):
        for copyrightHolder in info.split(';'):
            copyrightHolder = copyrightHolder.strip()
            if not copyrightHolder in self.copyrightHolders:
                self.copyrightHolders[copyrightHolder] = 0
            self.copyrightHolders[copyrightHolder] += 1
    
    def process_rtc(self, rtc_center):
        rtc_center = Vector(rtc_center)

        # Find Blender collection where the imported objects are located
        collection = self.collection_import
        if self.collection_import.children:
            # Find the first non-excluded child collection using the property <exclude> of <LayerCollection>
            # Why only the first? The other appear to be a phatom collection with the property <exclude> set to True.
            collection = next((c for c in self.layer_collection_import.children if not c.exclude)).collection
        
        for obj in collection.objects:
            # only top level objects are adjusted
            if not obj.parent:
                obj.location += rtc_center - self.centerCoords
    
    def finalize_glb_import(self, filepath, transformMatrix):
        removeFile(filepath)

        # Find Blender collection where the imported objects are located
        collection = self.collection_import
        if self.collection_import.children:
            # Find the first non-excluded child collection using the property <exclude> of <LayerCollection>
            # Why only the first? The other appear to be a phatom collection with the property <exclude> set to True.
            collection = next((c for c in self.layer_collection_import.children if not c.exclude)).collection

        # Move all objects from <collection> to the main collection
        for obj in collection.objects:
            collection.objects.unlink(obj)
            self.collection.objects.link(obj)
            if transformMatrix:
                obj.matrix_world = transformMatrix @ obj.matrix_world
        
        if self.collection_import.children:
            # Remove all child collections from <self.collection_import>
            for collection in self.collection_import.children:
                bpy.data.collections.remove(collection)
        
        self.num_imported_tiles += 1

    def joinObjects(self):
        # If a Blender Empty object was imported and it became the active object in the course of the processing,
        # after deleting it, Blender may not have any active object, which is required by <bpy.ops.object.join()>.
        # Therefore, we set the active object to the first object in the collection.
        bpy.context.view_layer.objects.active = self.collection.objects[0]

        if len(self.collection.objects) > 1:
            bpy.ops.object.join()
        
        joinedObject = self.collection.objects[0]
        joinedObject.name = self.threedTilesName

        # Removing doubles produces odd results in some cases, so it is disabled for now.
        #bpy.ops.object.mode_set(mode='EDIT')
        #bpy.ops.mesh.select_all(action='SELECT')
        #bpy.ops.mesh.remove_doubles()
        #bpy.ops.object.mode_set(mode='OBJECT')
    
    def patchGltfImporter(self):
        bv = bpy.app.version

        if (bv[0] == 3 and bv[1] == 6) or (bv[0] == 4 and bv[1] <= 2):
            from .gltf_patch import set_convert_functions_4_5, select_imported_objects_4_1
            from io_scene_gltf2.blender.imp.gltf2_blender_gltf import BlenderGlTF
            from io_scene_gltf2.blender.imp.gltf2_blender_scene import BlenderScene

            glTFImporter._offset = self.centerCoords
            
            self._set_convert_functions = BlenderGlTF.set_convert_functions
            BlenderGlTF.set_convert_functions = set_convert_functions_4_5
            
            self._select_imported_objects = BlenderScene.select_imported_objects
            BlenderScene.select_imported_objects = select_imported_objects_4_1
            
            self._gltfImporterPatched = (glTFImporter, BlenderGlTF, BlenderScene)
        elif (bv[0] == 4 and 3 <= bv[1]) or bv[0] > 4:
            from .gltf_patch import set_convert_functions_4_5, select_imported_objects_4_1
            from io_scene_gltf2.blender.imp.blender_gltf import BlenderGlTF
            from io_scene_gltf2.blender.imp.scene import BlenderScene
            
            glTFImporter._offset = self.centerCoords

            self._set_convert_functions = BlenderGlTF.set_convert_functions
            BlenderGlTF.set_convert_functions = set_convert_functions_4_5
            
            self._select_imported_objects = BlenderScene.select_imported_objects
            BlenderScene.select_imported_objects = select_imported_objects_4_1
            
            self._gltfImporterPatched = (glTFImporter, BlenderGlTF, BlenderScene)
    
    def cleanupGltfImporterPatching(self):
        glTFImporter, BlenderGlTF, BlenderScene = self._gltfImporterPatched
        
        # clean everything up after patching
        
        BlenderGlTF.set_convert_functions = self._set_convert_functions
        self._set_convert_functions = None

        delattr(glTFImporter, "_patch_convert_functions")
        delattr(glTFImporter, "_offset")
        
        BlenderScene.select_imported_objects = self._select_imported_objects
        self._select_imported_objects = None
        
        self._gltfImporterPatched = None