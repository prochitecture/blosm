import os

from mathutils import Vector

import bpy
from item_renderer.util import getFilepath
from mathutils.geometry import intersect_line_line_2d
from renderer import Renderer
from util.blender import createMeshObject, createCollection, getBmesh, setBmesh, loadMaterialsFromFile, \
    addGeometryNodesModifier, useAttributeForGnInput, createPolylineMesh

from .asset_store import AssetStore, AssetType, AssetPart


_location = Vector((0., 0., 0.))

#from renderer.curve_renderer import CurveRenderer
sidewalkWidth = 5.
crosswalkWidth = 3.6
stopLineWidth = 0.#1.5
# <c>rosswalk and <s>top <l>ine width
cslWidth = crosswalkWidth + stopLineWidth


class StreetRenderer:
    
    def __init__(self, app, itemRenderers):
        self.app = app
        self.itemRenderers = itemRenderers
        
        self.assetsDir = app.assetsDir
        self.assetPackageDir = app.assetPackageDir
        
        self.assetStore = AssetStore(app.assetInfoFilepathStreet)
        
        # initialize item renderers
        for itemRenderer in itemRenderers.values():
            itemRenderer.init(self)
        
        self.sectionRenderer = itemRenderers["Section"]
    
    def initTerrain(self):
        self.terrainObj = None
        
        terrain = self.app.terrain
        if terrain:
            terrain = terrain.terrain
            if terrain:
                # hide Blender terrain object
                terrain.hide_viewport = True
                terrain.hide_render = True
                self.terrainObj = terrain
    
    def prepare(self):
        self.collectionCenterlines = createCollection("Street centerlines", Renderer.collection)
        self.collectionIntersections = createCollection("Intersections", Renderer.collection)
        self.collection3dStreets = createCollection("Streets in 3D", Renderer.collection)
        
        # prepare item renderers
        for itemRenderer in self.itemRenderers.values():
            itemRenderer.prepare()
        
        nodeGroupNames = set((
            "Blosm Init Centerline",
            "Blosm Init Centerline for 3D",
            "blosm_terrain_area",
            "blosm_terrain_street_full"
            #"blosm_terrain_street_centerline"
        ))
        for itemRenderer in self.itemRenderers.values():
            itemRenderer.requestNodeGroups(nodeGroupNames)
        
        nodeGroups = dict(
            (nodeGroupName, bpy.data.node_groups[nodeGroupName]) for nodeGroupName in nodeGroupNames if nodeGroupName in bpy.data.node_groups
        )
        
        nodeGroupsToLoad = [nodeGroupName for nodeGroupName in nodeGroupNames if not nodeGroupName in bpy.data.node_groups]
        
        if nodeGroupsToLoad:
            with bpy.data.libraries.load(os.path.join(os.path.dirname(self.app.baseAssetPath), "prochitecture_streets.blend")) as (_, data_to):
                data_to.node_groups = nodeGroupsToLoad
        
        nodeGroups.update(
            (nodeGroup.name, nodeGroup) for nodeGroup in nodeGroupsToLoad
        )
        
        self.gnInitCenterline = nodeGroups["Blosm Init Centerline"]
        self.gnInitCenterline3d = nodeGroups["Blosm Init Centerline for 3D"]
        self.gnTerrainStreetFull = nodeGroups["blosm_terrain_street_full"]
        #self.gnTerrainStreetCenterline = nodeGroups["blosm_terrain_street_centerline"]
        
        for itemRenderer in self.itemRenderers.values():
            itemRenderer.setNodeGroups(nodeGroups)
        
        self.initTerrain()
        
        gnSidewalk = "blosm_sidewalk"
        gnLineItem = "blosm_line_item"
        gnSeparator = "blosm_roadway_separator"
        gnLamps = "blosm_street_lamps"
        gnTerrainPatches = "blosm_terrain_patches"
        gnProjectOnTerrain = "blosm_project_on_terrain"
        gnProjectTerrainPatches = "blosm_project_terrain_patches"
        gnMeshToCurve = "blosm_mesh_to_curve"
    
    def render(self, manager, data):
        # split neighbor street sections for the side lane transitions
        for sideLane in manager.transitionSideLanes:
            sideLane.splitAffectedSection()
        
        # render instances of the class <Bundle>
        for bundle in manager.iterBundles():
            for street in bundle.streetsHead:
                self.initStreet(street)
        
        # render instances of the class <Street>
        for street in manager.iterStreets():
            self.initStreet(street)

        # render instances of class <Intersection>
        intersectionRenderer = self.itemRenderers["Intersection"]
        for intersection in manager.majorIntersections.values():
            intersectionRenderer.renderItem(intersection)
        
        #
        # Finalize items: apply modifiers
        #
        for street in manager.iterStreets():
            street.obj3d = self.get3dStreetObj(street, _location)
            
            if street.head is street.tail:
                self.finalizeItem(street.head, 0)
            else:
                itemIndex = 1
                for item in street.iterItems():
                    self.finalizeItem(item, itemIndex)
                    itemIndex += 1
              
            if self.terrainObj:
                m = addGeometryNodesModifier(street.obj3d, self.gnTerrainStreetFull, "Streets on terrain")
                m["Input_2"] = self.terrainObj
        
        for itemRenderer in self.itemRenderers.values():
            itemRenderer.finalize()
    
    def initStreet(self, street):
        street.edgeIndexOffset = 0
        
        # Create a Blender object and BMesh for the <street> centeline. An instance of <Street> contains
        # at least one instance of <Section>, so <street.obj> and <street.bm> will be needed anyway.
        street.obj = self.getStreetCenterlineObj(street, _location)
        street.bm = getBmesh(street.obj)
        
        #if self.terrainObj:
        #    m = addGeometryNodesModifier(street.obj, self.gnTerrainStreetCenterline, "Streets on terrain")
        #    m["Input_2"] = self.terrainObj
        
        # Initialization is performed in two passes:
        # In the first pass we create a Blender mesh
        # in the second pass we set attributes
        
        #
        # (1) the first pass
        #
        if street.head is street.tail:
            self.initItemCenterline1(street.head, True)
        else:
            # the first BMesh vertex for the current street section
            street.bmVert = None
            for item in street.iterItems():
                self.initItemCenterline1(item, False)
        
        setBmesh(street.obj, street.bm)
        
        #
        # (2) the second pass
        #
        if street.head is street.tail:
            self.initItemCenterline2(street.head, 0)
        else:
            itemIndex = 1
            for item in street.iterItems():
                self.initItemCenterline2(item, itemIndex)
                itemIndex += 1
        
        addGeometryNodesModifier(street.obj, self.gnInitCenterline, "Init Centerline")
    
    def initItemCenterline1(self, item, singleItem):
        itemRenderer = self.itemRenderers.get(item.__class__.__name__)
        if itemRenderer:
            itemRenderer.initItemCenterline1(item, singleItem)

    def initItemCenterline2(self, item, itemIndex):
        itemRenderer = self.itemRenderers.get(item.__class__.__name__)
        if itemRenderer:
            itemRenderer.initItemCenterline2(item, itemIndex)
    
    def finalizeItem(self, item, itemIndex):
        itemRenderer = self.itemRenderers.get(item.__class__.__name__)
        if itemRenderer:
            itemRenderer.finalizeItem(item, itemIndex)
    
    def renderIntersections(self, manager):
        bm = getBmesh(self.intersectionAreasObj)
        
        for intersectionArea in manager.intersectionAreas:
            polygon = intersectionArea.polygon
            bm.faces.new(
                bm.verts.new(Vector((vert[0], vert[1], 0.))) for vert in polygon
            )
            
        for transition in manager.transitionSymLanes:
            polygon = transition.polygon
            bm.faces.new(
                bm.verts.new(Vector((vert[0], vert[1], 0.))) for vert in polygon
            )
            
            #self.processIntersectionSidewalks(intersectionArea, manager) FIXME
            
            #self.terrainRenderer.processIntersection(intersectionArea) FIXME
        
        
        setBmesh(self.intersectionAreasObj, bm)
        # apply the modifier <self.gnPolygons>
        m = addGeometryNodesModifier(self.intersectionAreasObj, self.gnPolygons, "Intersections")
        self.setMaterial(m, "Input_2", AssetType.material, None, AssetPart.pavement, "asphalt")
        
        setBmesh(self.intersectionSidewalksObj, self.intersectionSidewalksBm)
        # apply the modifier <self.gnPolygons>
        m = addGeometryNodesModifier(self.intersectionSidewalksObj, self.gnPolygons, "Intersection Sidewalks")
        self.setMaterial(m, "Input_2", AssetType.material, None, AssetPart.pavement, None)
        
        
        self.debugIntersectionArea(manager) # FIXME
        
        self.projectOnTerrain(self.intersectionAreasObj, self.gnProjectOnTerrain)
    
    def projectOnTerrain(self, obj, gnModifier):
        terrainObj = self.getTerrainObj()
        if terrainObj:
            m = addGeometryNodesModifier(obj, gnModifier, "Project on terrain")
            m["Input_2"] = terrainObj
        return terrainObj
    
    def getTerrainObj(self):
        terrain = self.app.terrain
        if terrain:
            terrain = terrain.terrain
            if terrain:
                terrain.hide_viewport = True
                terrain.hide_render = True
                return terrain
    
    def processIntersectionSidewalks(self, intersection, manager):
        """
        Process sidewalks around <intersection>
        """
        connectorsInfo = intersection.getConnectorsInfo()
        
        # iterate through all but last connectors in <connectorsInfo>
        for i in range(len(connectorsInfo)-1):
            self.processIntersectionSidewalk(intersection, connectorsInfo[i], connectorsInfo[i+1], manager)
        
        self.processIntersectionSidewalk(intersection, connectorsInfo[-1], connectorsInfo[0], manager)
    
    def processIntersectionSidewalk(self, intersection, connectorInfoL, connectorInfoR, manager):
        """
        Create a mesh for the part of the sidewalk between the street sections attached to the connectors
        described by <connectorInfoL> and <connectorInfoR>.
        
        Important note. The street section attached to the connector described by <connectorInfoR>
        is to the right from the one described by <connectorInfoL>.
        """
        
        bm = self.intersectionSidewalksBm
        polygon = intersection.polygon
        verts = None
        
        # get an instance for the street section attached to the connector described by <connectorInfo1>
        streetSectionL = manager.wayClusters[connectorInfoL[2]]\
            if connectorInfoL[1] else\
            manager.waySectionLines[connectorInfoL[2]]
        streetSectionR = manager.wayClusters[connectorInfoR[2]]\
            if connectorInfoR[1] else\
            manager.waySectionLines[connectorInfoR[2]]
        
        # index of the right point of the left connector
        indexL = (connectorInfoL[0] + 1) % intersection.numPoints
        # index of the left point of the right connector
        indexR = connectorInfoR[0]
        
        offsetToLeftL = not connectorInfoL[3]
        index1L = connectorInfoL[3]
        index2L = -2 if index1L else 1
        offsetDistance = (
            streetSectionL.getLeftBorderDistance()
            if offsetToLeftL else
            streetSectionL.getRightBorderDistance()
        ) + sidewalkWidth
        # Offset points for streetSectionL
        normal1L = streetSectionL.getNormal(index1L, offsetToLeftL)
        point1L = streetSectionL.offsetPoint(
            index1L,
            offsetToLeftL,
            offsetDistance,
            normal1L
        )
        normal2L = streetSectionL.getNormal(index2L, offsetToLeftL)
        point2L = streetSectionL.offsetPoint(
            index2L,
            offsetToLeftL,
            offsetDistance * self.getOffsetWeight(streetSectionL, index2L),
            normal2L
        )
        vectorL = point2L - point1L
        lengthL = vectorL.length
        
        offsetToLeftR = bool(connectorInfoR[3])
        index1R = connectorInfoR[3]
        index2R = -2 if index1R else 1
        offsetDistance = (
            streetSectionR.getLeftBorderDistance()
            if offsetToLeftR else
            streetSectionR.getRightBorderDistance()
        ) + sidewalkWidth
        normal1R = streetSectionR.getNormal(index1R, offsetToLeftR)
        point1R = streetSectionR.offsetPoint(
            index1R,
            offsetToLeftR,
            offsetDistance,
            normal1R
        )
        normal2R = streetSectionR.getNormal(index2R, offsetToLeftR)
        point2R = streetSectionR.offsetPoint(
            index2R,
            offsetToLeftR,
            offsetDistance * self.getOffsetWeight(streetSectionR, index2R),
            normal2R
        )
        vectorR = point2R - point1R
        lengthR = vectorR.length
        
        # Check if the ray with the origin at <point1L> and the direction towards <point2L> intersects
        # with the ray with the origin at <point1R> and the direction towards <point2R>.
        # Actually instead of calculating if the above rays intersect we simply enlarge the segments
        # <point1L>-<point2L> and <point1R>-<point2R> to the length <rayLength> and check if those
        # enlarged segments intersect.
        rayLength = 100.
        intersectionPoint = intersect_line_line_2d(
            point1L,
            point1L + rayLength/lengthL*vectorL,
            point1R,
            point1R + rayLength/lengthR*vectorR
        )
        if intersectionPoint:
            # The rays do intersect.
            
            # Check if the intersection point on the rays belongs to
            # the segment <point1L>-<point2L>
            if (intersectionPoint - point1L).length_squared <= lengthL*lengthL:
                # Calculate the foot of the perpendicular from <intersectionPoint> to
                # the street segment with the indices <index1L> and <index2L> and
                # offset at the border of the roadway in direction defined by
                # normals <normal1L> and <normal2L>
                point1L -= sidewalkWidth*normal1L
                point2L -= sidewalkWidth*normal2L
                # Unit vector along the segment <point1L>-<point2L>
                unitVector = vectorL/lengthL
                # Distance from <point1L> to the foot of the perpendicular
                trimLength = (intersectionPoint - point1L).dot(unitVector)
                if trimLength*trimLength <= (point2L-point1L).length_squared:
                    self.setTrimLength(
                        streetSectionL,
                        trimLength,
                        offsetToLeftL,
                        offsetToLeftL
                    )
                    # Vector to the foot of the perpendicular
                    foot = point1L + trimLength*unitVector
                    verts = [
                        bm.verts.new((polygon[indexL][0], polygon[indexL][1], 0.)),
                        bm.verts.new((foot[0], foot[1], 0.)),
                        bm.verts.new((intersectionPoint[0], intersectionPoint[1], 0.)),
                    ]
            # Check if the intersection point on the rays belongs to
            # the segment <point1R>-<point2R>
            if verts and (intersectionPoint - point1R).length_squared <= lengthR*lengthR:
                # Calculate the foot of the perpendicular from <intersectionPoint> to
                # the street segment with the indices <index1R> and <index2R> and
                # offset at the border of the roadway in direction defined by
                # normals <normal1R> and <normal2R>
                point1R -= sidewalkWidth*normal1R
                point2R -= sidewalkWidth*normal2R
                # Unit vector along the segment <point1R>-<point2R>
                unitVector = vectorR/lengthR
                # Distance from <point1R> to the foot of the perpendicular
                trimLength = (intersectionPoint - point1R).dot(unitVector)
                if trimLength*trimLength <= (point2R-point1R).length_squared:
                    self.setTrimLength(
                        streetSectionR,
                        trimLength,
                        offsetToLeftR,
                        not offsetToLeftR
                    )
                    # Vector to the foot of the perpendicular
                    foot = point1R + trimLength*unitVector
                    verts.append(bm.verts.new((foot[0], foot[1], 0.)))
                
        else:
            verts = [
                bm.verts.new((polygon[indexL][0], polygon[indexL][1], 0.)),
                bm.verts.new((point1L[0], point1L[1], 0.)),
                bm.verts.new((point1R[0], point1R[1], 0.))
            ]
        
        if verts:
            if indexL != indexR:
                verts.append(bm.verts.new((polygon[indexR][0], polygon[indexR][1], 0.)))
                
                if indexL < indexR:
                    if indexR - indexL > 1:
                        verts.extend(
                            bm.verts.new((polygon[i][0], polygon[i][1], 0.)) for i in range(indexR-1, indexL, -1)
                        )
                elif not (indexL==intersection.numPoints-1 and indexR==0):
                    indices = [] if indexL == intersection.numPoints-1 else\
                        list( range(indexL+1, intersection.numPoints) )
                    if indexR == 1:
                        indices.append(0)
                    elif indexR > 1:
                        indices.extend(i for i in range(indexR-1))
                    verts.extend(
                        bm.verts.new((polygon[i][0], polygon[i][1], 0.)) for i in reversed(indices)
                    )
            bm.faces.new(verts)
    
    def setAttributeOffsetWeight(self, inputObj, streetSection, startIndex):
        if self.streetSectionObjNames:
            streetSectionObjData = bpy.data.objects[ self.streetSectionObjNames[streetSection.index] ].data.attributes["offset_weight"].data
            inputObjData = inputObj.data.attributes["offset_weight"].data
            # "ss" stands for "street section"
            for inputObjIndex, ssObjIndex in zip(
                    range(startIndex, startIndex+len(streetSection.centerline)),
                    range(len(streetSection.centerline))
                ):
                # Copy the values of the attribute "offset_weight" from the street section object to <inputObj>
                inputObjData[inputObjIndex].value = streetSectionObjData[ssObjIndex].value
    
    def getOffsetWeight(self, streetSection, index):
        return bpy.data.objects[ self.streetSectionObjNames[streetSection.index] ].data.attributes["offset_weight"].data[index].value
    
    def setTrimLength(self, streetSection, trimLength, left, start):
        if left:
            streetSection.sidewalkL["Input_6" if start else "Input_7"] = trimLength
        else:
            streetSection.sidewalkR["Input_6" if start else "Input_7"] = trimLength
    
    def setOffsetWeights(self, obj, streetSection, pointIndexOffset):
        # Set offset weights. An offset weight is equal to
        # 1/sin(angle/2), where <angle> is the angle between <vec1> and <vec2> (see below the code)
        attributes = obj.data.attributes["offset_weight"].data
        centerline = streetSection.centerline
        numPoints = len(centerline)
        attributes[pointIndexOffset].value = attributes[pointIndexOffset+numPoints-1].value = 1.
        
        if numPoints > 2:
            vec1 = centerline[0] - centerline[1]
            vec1.normalize()
            for centerlineIndex, pointIndex in zip(range(1, numPoints-1), range(pointIndexOffset+1, pointIndexOffset+numPoints-1)):
                vec2 = centerline[centerlineIndex+1] - centerline[centerlineIndex]
                vec2.normalize()
                vec = vec1 + vec2
                vec.normalize()
                attributes[pointIndex].value = abs(1/vec.cross(vec2))
                vec1 = -vec2
    
    def setOffsetWeightsChunked(self, obj, streetSection, pointIndexOffset, transitionStart, transitionEnd):
        # Set offset weights. An offset weight is equal to
        # 1/sin(angle/2), where <angle> is the angle between <vec1> and <vec2> (see below the code)
        attributes = obj.data.attributes["offset_weight"].data
        centerline = streetSection.centerline
        numPoints = streetSection.numPoints
        pointIndex = pointIndexOffset
        attributes[pointIndex].value = attributes[pointIndex+numPoints-1].value = 1.
        
        if numPoints <=2:
            return
        
        # initialize variables
        index1 = 0
        vec1 = None
        
        if transitionStart and not transitionStart.offsetData[0]:
            # there is a side-lane transition and offset in the direction of <transitionStart.incoming>
            centerlinePrev = transitionStart.incoming.centerline
            numPointsPrev = len(centerlinePrev)
            index1 = transitionStart.offsetData[1]
            vec1 = centerlinePrev[index1-1] - centerlinePrev[index1]
            vec1.normalize()
            for centerlineIndex, pointIndex in zip(
                    range(index1, numPointsPrev-1),
                    range(pointIndex+1, pointIndex+numPointsPrev-index1)
                ):
                vec2 = centerlinePrev[centerlineIndex+1] - centerlinePrev[centerlineIndex]
                vec2.normalize()
                vec = vec1 + vec2
                vec.normalize()
                attributes[pointIndex].value = abs(1/vec.cross(vec2))
                vec1 = -vec2
            index1 = 0
        else:
            # The condition in the expession below means:
            # there is a side-lane transition and offset in the direction of <streetSection>
            index1 = transitionStart.offsetData[1]\
                if transitionStart and transitionStart.offsetData[0] else 1
            vec1 = centerline[index1-1] - centerline[index1]
            vec1.normalize()
        
        # The condition in the expession below means:
        # there is a side-lane transition and offset in the direction of <streetSection>
        index2 = transitionEnd.offsetData[1]\
            if transitionEnd and not transitionEnd.offsetData[0] else len(streetSection.centerline)-1
        
        for centerlineIndex, pointIndex in zip(
                range(index1, index2),
                range(pointIndex+1, pointIndex+index2-index1+1)
            ):
            vec2 = centerline[centerlineIndex+1] - centerline[centerlineIndex]
            vec2.normalize()
            vec = vec1 + vec2
            vec.normalize()
            attributes[pointIndex].value = abs(1/vec.cross(vec2))
            vec1 = -vec2
        
        if transitionEnd and transitionEnd.offsetData[0]:
            # there is a side-lane transition and offset in the direction of <transitionEnd.outgoing>
            centerlineNext = transitionEnd.outgoing.centerline
            index2 = transitionEnd.offsetData[1]
            for centerlineIndex, pointIndex in zip(
                    range(0, index2),
                    range(pointIndex+1, pointIndex+index2+1)
                ):
                vec2 = centerlineNext[centerlineIndex+1] - centerlineNext[centerlineIndex]
                vec2.normalize()
                vec = vec1 + vec2
                vec.normalize()
                attributes[pointIndex].value = abs(1/vec.cross(vec2))
                vec1 = -vec2
    
    def finalize(self):
        return
    
    def cleanup(self):
        self.terrainObj = None
    
    def getRoadwayClass(self, streetSection):
        return str(streetSection.totalLanes) + "_lanes"
    
    def setMaterial(self, modifier, modifierAttr, assetType, group, streetPart, cl):
        # get asset info for the material
        assetInfo = self.assetStore.getAssetInfo(
            assetType, group, streetPart, cl
        )
        if assetInfo:
            # set material
            material = self.getMaterial(assetInfo)
            if material:
                modifier[modifierAttr] = material
    
    def getMaterial(self, assetInfo):
        materialName = assetInfo["material"]
        material = bpy.data.materials.get(materialName)
        
        if not material:
            material = loadMaterialsFromFile(
                getFilepath(self, assetInfo),
                False,
                materialName
            )
            material = material[0] if material else None
            
        return material
    
    def getStreetCenterlineObj(self, street, location):
        obj = createMeshObject(
            "Polyline " + (street.getName() or "Street"),
            location = location,
            collection = self.collectionCenterlines
        )
        # create an attribute for the index of the street section
        obj.data.attributes.new("section_index", 'INT', 'EDGE')
        return obj
    
    def get3dStreetObj(self, street, location):
        obj = createMeshObject(
            street.obj.name[9:], # without <Polyline > in the beginning
            location = location,
            collection = self.collection3dStreets
        )
        m = addGeometryNodesModifier(obj, self.gnInitCenterline3d)
        m["Socket_3"] = street.obj
        # setting offset for points where two street sections share a position of their end points
        m["Socket_2"] = 4.
        # intersection at the start of the street
        if street.pred and street.pred.intersection.order == 4 and not street.pred.intersection.isMinor:
            m["Socket_4"] = street.pred.intersection.obj
        # intersection at the end of the street
        if street.succ and street.succ.intersection.order == 4 and not street.succ.intersection.isMinor:
            m["Socket_5"] = street.succ.intersection.obj
        return obj
    
    def debugIntersectionArea(self, manager):
        self.intersectionAreasObj.data.attributes.new("idx", 'INT', 'FACE')
        for idx,intersection in enumerate(manager.intersectionAreas):
            self.intersectionAreasObj.data.attributes["idx"].data[idx].value = idx
    
    def createModifiedPolylineMesh(self, bm, streetSection, transitionStart, transitionEnd):
        # <transitionStart> and <transitionEnd> are defined as follows:
        # transitionStart = streetSection.start if isinstance(streetSection.start, TransitionSideLane) else None
        # transitionEnd = streetSection.end if isinstance(streetSection.end, TransitionSideLane) else None
        # <transitionStart>, <transitionEnd> may therefore seem excessive.
        # They are calculated anyway before the call of this method and provided to prevent
        # an additional calculation of them
        
        # The condition in the expession below means:
        # there is a side-lane transition and offset in the direction of <streetSection>
        index1 = transitionStart.offsetData[1]\
            if transitionStart and transitionStart.offsetData[0] else 1
        
        # The condition in the expession below means:
        # there is a side-lane transition and offset in the direction of <streetSection>
        index2 = transitionEnd.offsetData[1]\
            if transitionEnd and not transitionEnd.offsetData[0] else len(streetSection.centerline)
        
        
        point = transitionStart.offsetData[2] if transitionStart else streetSection.centerline[0]
        prevVert = bm.verts.new((point[0], point[1], 0.))
        
        if transitionStart:
            if transitionStart.offsetData[0]:
                # there is a side-lane transition and offset in the direction of <streetSection>
                pass
            else:
                # there is a side-lane transition and offset in the direction of <transitionStart.incoming>
                streetSectionPrev = transitionStart.incoming
                for i in range(transitionStart.offsetData[1], len(streetSectionPrev.centerline)):
                    vert = bm.verts.new((streetSectionPrev.centerline[i][0], streetSectionPrev.centerline[i][1], 0.))
                    bm.edges.new((prevVert, vert))
                    prevVert = vert
        
        for i in range(index1, index2):
            vert = bm.verts.new((streetSection.centerline[i][0], streetSection.centerline[i][1], 0.))
            bm.edges.new((prevVert, vert))
            prevVert = vert
        
        if transitionEnd:
            if transitionEnd.offsetData[0]:
                # there is a side-lane transition and offset in the direction of <transitionEnd.outgoing>
                streetSectionNext = transitionEnd.outgoing
                for i in range(1, transitionEnd.offsetData[1]):
                    vert = bm.verts.new((streetSectionNext.centerline[i][0], streetSectionNext.centerline[i][1], 0.))
                    bm.edges.new((prevVert, vert))
                    prevVert = vert
            else:
                # there is a side-lane transition and offset in the direction of <streetSection>
                pass
            point = transitionEnd.offsetData[2]
            bm.edges.new(( prevVert, bm.verts.new((point[0], point[1], 0.)) ))
    
    def prepareModificationsForStreetSections(self, transition):
        """
        Prepare data for modification of the street sections due to the side-lane transition.
        Lengthen the street section with the smaller number of lanes and
        shorten the street section with  the larger number of lanes.
        """
        streetSection1 = transition.incoming
        
        # <streetSection1> precedes <streetSection2>
        streetSection2 = transition.outgoing
        
        if transition.totalLanesIncreased:
            centerline2 = streetSection2.centerline
            accumulatedLength = 0.
            for i in range(1, streetSection2.numPoints):
                vec = centerline2[i] - centerline2[i-1]
                # length of the current segment
                l = vec.length
                if accumulatedLength+l > transition.lengthOffset:
                    transition.offsetData = (
                        # offset in the direction of <centerline2> (i.e. in the outgoing direction)
                        True,
                        i,
                        centerline2[i-1] + (transition.lengthOffset - accumulatedLength)/l * vec
                    )
                    # update the number of points in <streetSection1> and <streetSection2>
                    streetSection1.numPoints += i
                    if i > 1:
                        streetSection2.numPoints -= i-1
                    break
                accumulatedLength += l
        else:
            centerline1 = streetSection1.centerline
            accumulatedLength = 0.
            for i in range(streetSection1.numPoints-2, -1, -1):
                vec = centerline1[i+1] - centerline1[i]
                # length of the current segment
                l = vec.length
                if accumulatedLength+l > transition.lengthOffset:
                    transition.offsetData = (
                        # offset in the direction of <centerline1> (i.e. in the incoming direction)
                        False,
                        i+1,
                        centerline1[i] + (l-transition.lengthOffset + accumulatedLength)/l * vec
                    )
                    # update the number of points in <streetSection1> and <streetSection2>
                    streetSection2.numPoints += streetSection1.numPoints-i-1
                    streetSection1.numPoints = i+2
                    break
                accumulatedLength += l
    
    def generateSideLaneTransition(self, obj, transition, laneOnRight):
        # get the street section <streetSection> for which the side-lane transition will be generated
        streetSection, streetSectionMoreLanes =\
            (transition.incoming, transition.outgoing)\
            if transition.totalLanesIncreased else\
            (transition.outgoing, transition.incoming)
        
        m = addGeometryNodesModifier(obj, self.gnSideLaneTransition, "Side-Lane Transition")
        m["Input_2"] = transition.length
        m["Input_3"] = streetSectionMoreLanes.width - streetSection.width
        m["Input_4"] = streetSection.width
        m["Input_5"] = streetSection.offset
        # Reverse the position of the transition lane from right to left or from left to right
        # if the total number of lanes is decreased
        m["Input_6"] = laneOnRight if transition.totalLanesIncreased else not laneOnRight
        self.setMaterial(m, "Input_7", AssetType.material, "demo", AssetPart.side_lane_transition, "default")
        m["Input_8"] = streetSection.streetSectionIndex
        useAttributeForGnInput(m, "Input_9", "offset_weight")
        m["Input_10"] = transition.totalLanesIncreased


class TerrainPatchesRenderer:
    # A renderer for the terrain patches made between street graph cycles
    
    def __init__(self, streetRenderer):
        self.obj = obj = createMeshObject("Terrain Patches", collection = Renderer.collection)
        obj.data.attributes.new("offset_l", 'FLOAT', 'POINT')
        obj.data.attributes.new("offset_r", 'FLOAT', 'POINT')
        obj.data.attributes.new("offset_weight", 'FLOAT', 'POINT')
        
        m = addGeometryNodesModifier(obj, streetRenderer.gnTerrainPatches, "Terrain patches")
        useAttributeForGnInput(m, "Input_2", "offset_l")
        useAttributeForGnInput(m, "Input_3", "offset_r")
        useAttributeForGnInput(m, "Input_4", "offset_weight")
        
        # apply the modifier <self.gnPolygons>
        m = addGeometryNodesModifier(obj, streetRenderer.gnPolygons, "Material for terrain patches")
        streetRenderer.setMaterial(m, "Input_2", AssetType.material, None, AssetPart.vegetation, "grass")
        
        terrainObj = streetRenderer.getTerrainObj()
        if terrainObj:
            m = addGeometryNodesModifier(obj, streetRenderer.gnProjectTerrainPatches, "Project terrain patches")
            m["Input_2"] = terrainObj
        
        self.bm = getBmesh(obj)
    
    def processStreetCenterline(self, streetSection):
        createPolylineMesh(None, self.bm, streetSection.centerline)
    
    def setAttributes(self, manager, renderer):
        setBmesh(self.obj, self.bm)
        
        pointIndex = 0
        
        for streetSection in manager.waySectionLines.values():
            offsetL = -streetSection.getLeftBorderDistance()
            offsetR = streetSection.getLeftBorderDistance()
            renderer.setAttributeOffsetWeight(self.obj, streetSection, pointIndex)
            for _ in streetSection.centerline:
                self.obj.data.attributes["offset_l"].data[pointIndex].value = offsetL
                self.obj.data.attributes["offset_r"].data[pointIndex].value = offsetR
                pointIndex += 1
        
        for streetSection in manager.wayClusters.values():
            offsetL = -streetSection.getLeftBorderDistance()
            offsetR = streetSection.getRightBorderDistance()
            renderer.setAttributeOffsetWeight(self.obj, streetSection, pointIndex)
            for _ in streetSection.centerline:
                self.obj.data.attributes["offset_l"].data[pointIndex].value = offsetL
                self.obj.data.attributes["offset_r"].data[pointIndex].value = offsetR
                pointIndex += 1
    
    def processDeadEnds(self, manager):
        
        for streetSection in manager.waySectionLines.values():
            centerline = streetSection.centerline
            directionVector = None
            if not streetSection.startConnected:
                directionVector = _getDirectionVector(centerline[0], centerline[1])
                
                self.offsetEdgePoint(
                    centerline[0],
                    directionVector,
                    streetSection.getLeftBorderDistance()
                )
            
            if not streetSection.endConnected:
                if len(centerline)>2 or not directionVector:
                    directionVector = _getDirectionVector(centerline[-2], centerline[-1])
                self.offsetEdgePoint(
                    centerline[-1],
                    directionVector,
                    streetSection.getLeftBorderDistance()
                )
        
        for streetSection in manager.wayClusters.values():
            centerline = streetSection.centerline
            directionVector = None
            
            if not streetSection.startConnected:
                directionVector = _getDirectionVector(centerline[0], centerline[1])
                self.offsetEdgePoint(
                    centerline[0],
                    directionVector,
                    streetSection.getLeftBorderDistance(),
                    streetSection.getRightBorderDistance()
                )
            
            if not streetSection.endConnected:
                if len(centerline)>2 or not directionVector:
                    directionVector = _getDirectionVector(centerline[-2], centerline[-1])
                self.offsetEdgePoint(
                    centerline[-1],
                    directionVector,
                    streetSection.getLeftBorderDistance(),
                    streetSection.getRightBorderDistance()
                )
    
    def processIntersection(self, intersection):
        polygon = intersection.polygon
        
        connectorsInfo = intersection.getConnectorsInfo()
        
        #
        # all but the final segments
        #
        for i in range(len(connectorsInfo)-1):
            polylineStartIndex, polylineEndIndex = connectorsInfo[i][0]+1, connectorsInfo[i+1][0]
            # The condition below is used to exclude the case when a connector is directly
            # followed by another connector
            if polylineStartIndex != polylineEndIndex:
                # Create a polyline out of <polygon> that starts at the polygon's point
                # with the index <polylineStartIndex> and ends at the polygon's point
                # with the index <polylineEndIndex>
                prevVert = self.bm.verts.new((polygon[polylineStartIndex][0], polygon[polylineStartIndex][1], 0.))
                for i in range(polylineStartIndex+1, polylineEndIndex+1):
                    vert = self.bm.verts.new((polygon[i][0], polygon[i][1], 0.))
                    self.bm.edges.new((prevVert, vert))
                    prevVert = vert
        
        #
        # the final segment
        #
        # The condition below is used to exclude the case when a connector is directly
        # followed by another connector
        if connectorsInfo[-1][0] != intersection.numPoints-1 or connectorsInfo[0][0]:
            indices = list( range(connectorsInfo[-1][0]+1, intersection.numPoints) )
            indices.extend(i for i in range(connectorsInfo[0][0]+1))
            prevVert = self.bm.verts.new((polygon[indices[0]][0], polygon[indices[0]][1], 0.))
            for i in range(1, len(indices)):
                vert = self.bm.verts.new((
                    polygon[ indices[i] ][0],
                    polygon[ indices[i] ][1],
                    0.
                ))
                self.bm.edges.new((prevVert, vert))
                prevVert = vert
    
    def offsetEdgePoint(self, point, directionVector, offsetLeft, offsetRight=None):
        if offsetRight is None:
            offsetRight = offsetLeft
        
        # point to left from the centerline
        v1 = point - offsetLeft * directionVector
        # point to right from the centerline
        v2 = point + offsetRight * directionVector
        self.bm.edges.new((
            self.bm.verts.new((v1[0], v1[1], 0.)),
            self.bm.verts.new((v2[0], v2[1], 0.))
        ))
    
    def addExtent(self, app):
        # margin in meter
        margin = 2.
        
        minX = app.minX - margin
        minY = app.minY - margin
        maxX = app.maxX + margin
        maxY = app.maxY + margin
        
        # verts
        v1 = self.bm.verts.new((minX, minY, 0.))
        v2 = self.bm.verts.new((maxX, minY, 0.))
        v3 = self.bm.verts.new((maxX, maxY, 0.))
        v4 = self.bm.verts.new((minX, maxY, 0.))
        
        # edges
        self.bm.edges.new((v1, v2))
        self.bm.edges.new((v2, v3))
        self.bm.edges.new((v3, v4))
        self.bm.edges.new((v4, v1))


def _getDirectionVector(point1, point2):
    directionVector = point2 - point1
    directionVector.normalize()
    # vector perpendicular to <directionVector>
    return Vector((directionVector[1], -directionVector[0]))