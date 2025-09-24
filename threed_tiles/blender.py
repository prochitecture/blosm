from os.path import basename, join as joinStrings, exists as pathExists
from os import remove as removeFile
from math import radians, atan, sqrt, pi
import re
from operator import itemgetter
import sys

import bpy
from mathutils import Matrix

from util.blender import createEmptyObject, createCollection


class BlenderRenderer:
    
    def __init__(self, threedTilesName, join3dTilesObjects, instanceName):
        self.threedTilesName = threedTilesName
        self.join3dTilesObjects = join3dTilesObjects
        self.instanceName = instanceName
        self.importedObjects = []
        
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
        self.collection = createCollection(self.threedTilesName)
        
        self.centerCoords = manager.fromGeographic(manager.centerLat, manager.centerLon, 0.)
        
        if "io_scene_gltf2" in sys.modules:
            self.patchGltfImporter()
    
    def finalize(self, manager):
        if not self.importedObjects:
            self.collection = None
            return
        
        #
        # tranformation matrix
        #
        centerCoords = self.centerCoords
        # lat = radians(manager.centerLat - 90.) # gives incorrect result for the expression below
        lat = atan(centerCoords[2]/sqrt(centerCoords[0]*centerCoords[0] + centerCoords[1]*centerCoords[1]))
        # Rotate the mesh, so it will point to the north pole. The rotations are around Z and X axes
        matrix = Matrix.Rotation(lat-pi/2., 4, 'X') @ Matrix.Rotation(radians(-90. - manager.centerLon), 4, 'Z')
        
        locationsAfterRotation = [(matrix @ obj.location) for obj in self.importedObjects]
        
        # find the lowest Z-coordinate if <self.calculateHeightOffset>
        heightOffset = min(location[2] for location in locationsAfterRotation)\
            if self.calculateHeightOffset else\
            self.heightOffset
        if self.calculateHeightOffset:
            self.heightOffset = heightOffset
        
        # select the imported objects
        bpy.ops.object.select_all(action='DESELECT')
        for obj in self.importedObjects:
            obj.select_set(True)
        
        # apply possible rotation after Blender's glTF importer
        bpy.ops.object.transform_apply(location=False, rotation=True, scale=False)
        
        if self.join3dTilesObjects:
            self.joinObjects()
            
            # set the origin of the resulting Blender object at <centerCoords>
            _cursorLocation = bpy.context.scene.cursor.location.copy()
            bpy.context.scene.cursor.location = (0., 0., 0.) if self._gltfImporterPatched else centerCoords
            bpy.ops.object.origin_set(type='ORIGIN_CURSOR')
            bpy.context.scene.cursor.location = _cursorLocation
            
            joinedObject = self.importedObjects[-1]
            #location = locationsAfterRotation[-1]
            #location[2] -= heightOffset
            joinedObject.matrix_local = matrix#Matrix.Translation(location) @ matrix
        else:
            if not self._gltfImporterPatched:
                # rotate the vector <centerCoords>
                centerCoords = matrix @ centerCoords
            for obj, location in zip(self.importedObjects, locationsAfterRotation):
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
        
        numImportedTiles = len(self.importedObjects)
        
        self.importedObjects.clear()
        self.collection = None
        
        if self._gltfImporterPatched:
            vnode, BlenderScene = self._gltfImporterPatched
            
            # clean everything up after patching
            
            vnode.get_node_trs = self._get_node_trs
            self._get_node_trs = None
            
            BlenderScene.select_imported_objects = self._select_imported_objects
            
            self._gltfImporterPatched = None
        
        return numImportedTiles
    
    def renderGlb(self, manager, uri, path, cacheContent):
        context = bpy.context
        
        filePath = joinStrings(
            manager.tilesDir,
            basename(path) if cacheContent else ("current_file_" + self.instanceName + ".glb" if self.instanceName else "current_file.glb")
        )
        
        if cacheContent:
            if not pathExists(filePath):
                fileContent = manager.download(uri)
                with open(filePath, 'wb') as f:
                    f.write(fileContent)
            bpy.ops.import_scene.gltf(filepath=filePath)
        else:
            fileContent = manager.download(uri)
            # check if <fileContent> contains copyright information
            match = re.search(self.licenseRePattern, fileContent)
            if match:
                self.processCopyrightInfo(match.group(1).decode('utf-8'))
            with open(filePath, 'wb') as f:
                f.write(fileContent)
            bpy.ops.import_scene.gltf(filepath=filePath)
            removeFile(filePath)
        
        importedObject = context.object
        # unlink <importedObject> from its collection, there can be more than one colection
        for collection in importedObject.users_collection:
            collection.objects.unlink(importedObject)
        # link <importedObject> to <self.collection>
        self.collection.objects.link(importedObject)
        self.importedObjects.append(importedObject)

    def renderB3dm(self, manager, uri, path, cacheContent):
        import numpy
        from .py3dtiles.tileset.content.tile_content_reader import read_array
        
        context = bpy.context
        
        filePath = joinStrings(
            manager.tilesDir,
            basename(path)[:-4] + "glb" if cacheContent else ("current_file_" + self.instanceName + ".glb" if self.instanceName else "current_file.glb")
        )
        
        if cacheContent:
            if not pathExists(filePath):
                fileContent = manager.download(uri)
                with open(filePath, 'wb') as f:
                    f.write(fileContent)
            bpy.ops.import_scene.gltf(filepath=filePath)
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
            rtc_center = b3dmContent.body.feature_table.header.data.get("RTC_CENTER")
            
            # set position
            if rtc_center:
                gltfContent.header["nodes"] = [
                    dict(
                        mesh = 0,
                        translation = (rtc_center[0], rtc_center[2], -rtc_center[1])
                    )
                ]
            
            with open(filePath, 'wb') as f:
                f.write(gltfContent.to_array())
            
            bpy.ops.import_scene.gltf(filepath=filePath)
            removeFile(filePath)
        
        importedObject = context.object
        self.collection.objects.link(importedObject)
        self.importedObjects.append(importedObject)
        context.scene.collection.objects.unlink(importedObject)
    
    def processCopyrightInfo(self, info):
        for copyrightHolder in info.split(';'):
            copyrightHolder = copyrightHolder.strip()
            if not copyrightHolder in self.copyrightHolders:
                self.copyrightHolders[copyrightHolder] = 0
            self.copyrightHolders[copyrightHolder] += 1
    
    def joinObjects(self):
        if len(self.importedObjects) > 1:
            bpy.ops.object.join()
        joinedObject = self.importedObjects[-1]
        joinedObject.name = self.threedTilesName
        bpy.context.view_layer.objects.active = joinedObject
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.remove_doubles()
        bpy.ops.object.mode_set(mode='OBJECT')
    
    def patchGltfImporter(self):
        centerCoords = self.centerCoords
        # <_get_node_trs> will be set later in the code
        _get_node_trs = None
        
        def get_node_trs(gltf, pynode):
            t, r, s = _get_node_trs(gltf, pynode)
            
            # Unit conversion factor in (Blender units) per meter
            u = 1.0 / bpy.context.scene.unit_settings.scale_length
            
            t = t/u
            t = u * (t - centerCoords)
            return t, r, s
        
        
        bv = bpy.app.version
        if (bv[0] == 3 and bv[1] == 6) or (bv[0] == 4 and bv[1] <= 2):
            from .gltf_patch import select_imported_objects_4_1
            import io_scene_gltf2.blender.imp.gltf2_blender_vnode as vnode
            from io_scene_gltf2.blender.imp.gltf2_blender_scene import BlenderScene
            
            _get_node_trs = self._get_node_trs = vnode.get_node_trs
            vnode.get_node_trs = get_node_trs
            
            self._select_imported_objects = BlenderScene.select_imported_objects
            BlenderScene.select_imported_objects = select_imported_objects_4_1
            
            self._gltfImporterPatched = (vnode, BlenderScene)
        elif bv[0] == 4 and 3 <= bv[1]:
            from .gltf_patch import select_imported_objects_4_1
            import io_scene_gltf2.blender.imp.vnode as vnode
            from io_scene_gltf2.blender.imp.scene import BlenderScene
            
            _get_node_trs = self._get_node_trs = vnode.get_node_trs
            vnode.get_node_trs = get_node_trs
            
            self._select_imported_objects = BlenderScene.select_imported_objects
            BlenderScene.select_imported_objects = select_imported_objects_4_1
            
            self._gltfImporterPatched = (vnode, BlenderScene)