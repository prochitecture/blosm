"""
This file is a part of Blosm addon for Blender.
Copyright (C) 2014-2018 Vladimir Elistratov
prokitektura+support@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import bpy, bmesh


def makeActive(obj, context=None):
    if not context:
        context = bpy.context
    obj.select_set(True)
    context.view_layer.objects.active = obj


def createMeshObject(name, location=(0., 0., 0.), mesh=None, collection=None):
    if not mesh:
        mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    obj.location = location
    if not collection:
        collection = bpy.context.scene.collection
    collection.objects.link(obj)
    return obj


def createEmptyObject(name, location, hide=False, collection=None, **kwargs):
    obj = bpy.data.objects.new(name, None)
    obj.location = location
    obj.hide_viewport = hide
    obj.hide_select = hide
    obj.hide_render = True
    if kwargs:
        for key in kwargs:
            setattr(obj, key, kwargs[key])
    if not collection:
        collection = bpy.context.scene.collection
    collection.objects.link(obj)
    return obj


def createCollection(name, parent=None, hide_viewport=False, hide_select=False, hide_render=False):
    collection = bpy.data.collections.new(name)
    if not parent:
        parent = bpy.context.scene.collection
    parent.children.link(collection)
    collection.hide_viewport = hide_viewport
    collection.hide_select = hide_select
    collection.hide_render = hide_render
    return collection


def getBmesh(obj):
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    return bm


def setBmesh(obj, bm):
    bm.to_mesh(obj.data)
    bm.free()
    

def pointNormalUpward(face):
    if face.normal.z < 0.:
        face.normal_flip()


def createDiffuseMaterial(name, color):
    material = bpy.data.materials.new(name)
    material.diffuse_color = (color[0], color[1], color[2], 1.)
    return material


def loadMeshFromFile(filepath, name):
    """
    Loads a Blender mesh with the given <name> from the .blend file with the given <filepath>
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.meshes>
        data_to.meshes = [name]
    return data_to.meshes[0]


def loadParticlesFromFile(filepath, name):
    """
    Loads Blender particles settings with the given <name> from the .blend file
    with the given <filepath>
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.particles>
        data_to.particles = [name]
    return data_to.particles[0]


def loadCollectionFromFile(filepath, name):
    """
    Loads a Blender collection with the given <name> from the .blend file with the given <filepath>
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.meshes>
        data_to.collections = [name]
    return data_to.collections[0]


def linkCollectionFromFile(filepath, name):
    """
    Links a Blender collection with the given <name> from the .blend file with the given <filepath>
    """
    with bpy.data.libraries.load(filepath, link=True) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.meshes>
        data_to.collections = [name]
    return data_to.collections[0]


def loadSceneFromFile(filepath, name):
    """
    Loads a Blender scene with the given <name> from the .blend file with the given <filepath>
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.meshes>
        data_to.scenes = [name]
    return data_to.scenes[0]


def loadGroupFromFile(filepath, name):
    """
    Loads a Blender group with the given <name> from the .blend file with the given <filepath>
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.meshes>
        data_to.groups = [name]
    return data_to.groups[0]


def loadImagesFromFile(filepath, *names):
    """
    Loads images with <names> from the .blend file with the given <filepath>.
    If an image name is available at <bpy.data.images>, the image won't be loaded
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.images>
        data_to.images = [
            name for name in names if not name in bpy.data.images and name in data_from.images
        ]


def loadNodeGroupsFromFile(filepath, *names):
    """
    Loads node groups with <names> from the .blend file with the given <filepath>.
    If a node group is available at <bpy.data.node_groups>, the node groups won't be loaded
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.images>
        data_to.node_groups = [
            name for name in names if not name in bpy.data.node_groups and name in data_from.node_groups
        ]


def appendNodeGroupFromFile(filepath, name):
    """
    Appends node groups with <name> from the .blend file with the given <filepath>.
    """
    with bpy.data.libraries.load(filepath) as (_, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.node_groups>
        data_to.node_groups = [name]
    return data_to.node_groups[0]


def loadTextFromFile(filepath, name):
    """
    Loads a Blender text with the given <name> from the .blend file with the given <filepath>
    """
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.meshes>
        data_to.texts = [name]
    return data_to.texts[0]


def appendObjectsFromFile(filepath, collection, *names):
    with bpy.data.libraries.load(filepath) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.objects>
        data_to.objects = list(names)
    if collection:
        # append the objects to the Blender scene
        for obj in data_to.objects:
            if obj:
                    collection.objects.link(obj)
                    obj.select_set(False)
    # return the appended Blender objects
    return data_to.objects


def linkObjectFromFile(filepath, collection, name):
    with bpy.data.libraries.load(filepath, link=True) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.objects>
        data_to.objects = [name]
    obj = data_to.objects[0]
    if collection:
        # link <obj> to <collection>
        collection.objects.link(obj)
        obj.select_set(False)
    # return the appended Blender objects
    return obj


def getMaterialIndexByName(obj, name, filepath):
    """
    Check if Blender material with the <name> is already set for <obj>,
    if not, check if the material is available in bpy.data.material
    (if yes, append it to <obj>),
    if not, load the material with the <name> from the .blend with the given <filepath>
    and append it to <obj>.
    """
    if name in obj.data.materials:
        material = obj.data.materials[name]
        # find index of the material
        for materialIndex,m in enumerate(obj.data.materials):
            if material == m:
                break
    elif name in bpy.data.materials:
        materialIndex = len(obj.data.materials)
        obj.data.materials.append( bpy.data.materials[name] )
    else:
        with bpy.data.libraries.load(filepath) as (data_from, data_to):
            data_to.materials = [name]
        material = data_to.materials[0]
        materialIndex = len(obj.data.materials)
        obj.data.materials.append(material)
    return materialIndex


def getMaterialByName(obj, name, filepath=None):
    """
    Check if Blender material with the <name> is already set for <obj>,
    if not, check if the material is available in bpy.data.material,
    if not and if <filepath> is provided, load the material with the <name>
    from the .blend with the given <filepath>.
    The material is NOT appended to <obj>.
    """
    material = None
    if name in obj.data.materials:
        material = obj.data.materials[name]
    elif name in bpy.data.materials:
        material = bpy.data.materials[name]
    elif filepath:
        with bpy.data.libraries.load(filepath) as (data_from, data_to):
            data_to.materials = [name]
        material = data_to.materials[0]
    return material


def loadMaterialsFromFile(filepath, link, *names):
    with bpy.data.libraries.load(filepath, link=link) as (data_from, data_to):
        # a Python list (not a Python tuple!) must be set to <data_to.objects>
        data_to.materials = list(names)
    return data_to.materials


def getModifier(obj, modifierType):
    for m in obj.modifiers:
        if m.type == modifierType:
            break
    else:
        m = None
    return m


def addShrinkwrapModifier(obj, target, offset):
    m = obj.modifiers.new(name="Shrinkwrap", type='SHRINKWRAP')
    m.wrap_method = "PROJECT"
    m.use_positive_direction = False
    m.use_negative_direction = True
    m.use_project_z = True
    m.target = target
    m.offset = offset


def addGeometryNodesModifier(obj, nodeGroup, modifierName=''):
    m = obj.modifiers.new(modifierName, "NODES")
    m.node_group = nodeGroup
    return m


def loadImage(fileName, directory):
    image = bpy.data.images.get(fileName if directory else os.path.basename(fileName))
    if not image:
        # absolute path!
        imagePath = os.path.join(directory, fileName) if directory else fileName
        try:
            image = bpy.data.images.load(imagePath)
        except Exception:
            print("Unable to load the image %s" % imagePath)
    return image


def useAttributeForGnInput(modifier, inputId, attributeName):
    # Set "_use_attribute" to 1 to use geometry attributes instead of
    # using manually entered input values
    modifier[inputId + "_use_attribute"] = 1
    # set "_attribute_name" to the related mesh attribute of the Blender object
    modifier[inputId + "_attribute_name"] = attributeName


def createPolylineMesh(obj, bm, polyline, prevVert):
    if obj:
        bm = getBmesh(obj)
    
    if not prevVert:
        prevVert = bm.verts.new((polyline[0][0], polyline[0][1], 0.))
    for i in range(1, len(polyline)):
        vert = bm.verts.new((polyline[i][0], polyline[i][1], 0.))
        bm.edges.new((prevVert, vert))
        prevVert = vert
    
    if obj:
        setBmesh(obj, bm)
    else:
        return prevVert