#-----------------------------------------------------------------------------
# Polygon triangulation by partitioning into y-monotone pieces and subsequent 
# triangulation of these monotone polygons. The algorithms are based on the
# description in:
# 
# Computational Geometry: Algorithms and Applications
# Mark de Berg, Otfried Cheong, Marc van Kreveld, Mark Overmars
# Third Edition, Springer 2008, ISBN 978-3-540-77973-5, Pages 49-58.
#-----------------------------------------------------------------------------

import heapq
from mathutils import Vector
from .MonotoneStates import MonotoneStates
from .EdgeTree import EdgeData
from .Vertex import VertexType, ChainType

class PolygonTriangulation():
    def __init__(self):
        self.sweepStates = None
        self.triangles = []

    def triangulate(self,vPolygon, vHoles=[]):
        """
        Triangulates the polygon with holes.
        polygon: List of vertices of class Vertex in counter-clockwise order
        holes:   List of holes, each hole a list of vertices of class Vertex in 
                 clockwse order. The holes may not intersect each other and
                 have to be completely inside the polygon.
        return:  List of tuples, each of them contains three coordinates of class 
                 mathutils.Vector of a triangle.
        """
        # clear members, when already used once.
        # if self.sweepStates: self.sweepStates = None
        if self.triangles: self.triangles = []

        # Initialize doubly linked lists for polygon and holes and
        # compute type of vertices. Count vertices of type SPLIT or
        # MERGE to decide, whether the polygon is already y-monotone.
        self.sweepStates = MonotoneStates()
        splitOrMergeCount = self.sweepStates.initialize(vPolygon,vHoles)

        # If there are no vertices of type SPLIT or MERGE we can skip
        # the partitioning into y-monotone pieces.
        if splitOrMergeCount:
            # Partitioning into y-monotone pieces
            sweepEvents = self.sweepStates.sweepEvents
            while sweepEvents:
                event = heapq.heappop(sweepEvents)
                if event.type == VertexType.START:
                    self.handleStartVertex(event)
                elif event.type == VertexType.REGULAR:
                    self.handleRegularVertex(event)
                elif event.type == VertexType.MERGE:
                    self.handleMergeVertex(event)
                elif event.type == VertexType.SPLIT:
                    self.handleSplitVertex(event)
                elif event.type == VertexType.END:
                    self.handleEndVertex(event)
                else:
                    raise Exception("should not happen")
        else:
            self.sweepStates.diagList.diagonals.append(self.sweepStates.sweepEvents[0].edge)

        # Triangulate the y-monotone polygons
        for startEdge in self.sweepStates.diagList.diagonals:
            if not startEdge.triangulated:
                self.triangulateMonotone(startEdge)

        # Return the list of triangles, each as tuple of three vertices of class Vertex
        # Sometimes duplicate triangles are created, don't know why. Remove them here.
        self.triangles = list(set(self.triangles) )
        return self.triangles

    def handleStartVertex(self,event):
        # Insert event-edge ei into edge tree and set helper(ei) to event.
        EdgeData.currentY = event[1]
        event.edge.helper = event
        self.sweepStates.edgeTree.insertEdge(event.edge)

    def handleRegularVertex(self,event):
        predV = event.pred
        polyIsLeft = predV[1] < event[1] or (predV[1] == event[1] and predV[0] > event[0])
        # If the interior of the polygon lies on the left of the event
        if polyIsLeft:
            # Search in edge tree to find the edge ei directly left of event
            EdgeData.currentY = event.y
            toTheLeft = self.sweepStates.edgeTree.findToLeft(event)
            # if helper(ej) is a merge vertex
            if toTheLeft.helper.type == VertexType.MERGE:
                # Insert the diagonal connecting the event to helper(ej)
                # into the list of diagonals.
                self.sweepStates.diagList.insertDiagonal(event, toTheLeft)
            # helper(ej) <- event
            toTheLeft.helper = event
        else:   # the interior of the polygon lies on the right of the event
            # if helper(ei−1) is a merge vertex
            predEdge = event.pred.edge
            if predEdge.helper.type == VertexType.MERGE:
                # Insert the diagonal connecting the event to helper(ei−1)
                # into the list of diagonals.
                 self.sweepStates.diagList.insertDiagonal(event, predEdge)
            # Delete ei−1 from edge tree 
            EdgeData.currentY = event.y
            self.sweepStates.edgeTree.removeEdge(predEdge)
            # Insert ei in edge tree and set helper(ei) to event
            event.edge.helper = event
            self.sweepStates.edgeTree.insertEdge(event.edge)

    def handleMergeVertex(self,event):
        # if helper(ei−1) is a merge vertex
        predEdge = event.pred.edge
        if predEdge.helper.type == VertexType.MERGE:
            # Insert the diagonal connecting the event to helper(ei−1)
            # into the list of diagonals.
            self.sweepStates.diagList.insertDiagonal(event, predEdge)
        # Delete ei−1 from edge tree
        EdgeData.currentY = event.y
        self.sweepStates.edgeTree.removeEdge(event.edge)
        # Search in edge tree to find the edge ei directly left of event
        toTheLeft = self.sweepStates.edgeTree.findToLeft(event)
        # if helper(ej) is a merge vertex
        if toTheLeft.helper.type == VertexType.MERGE:
            # Insert the diagonal connecting the event to helper(ei)
            # into the list of diagonals.
            self.sweepStates.diagList.insertDiagonal(event, toTheLeft)
        # helper(ej) <- event
        toTheLeft.helper = event

    def handleSplitVertex(self,event):
        # Search in edge tree to find the edge ei directly left of event
        EdgeData.currentY = event.y
        edgeToTheLeft = self.sweepStates.edgeTree.findToLeft(event)
        # Insert the diagonal connecting the event to helper(ei)
        # into the list of diagonals.
        self.sweepStates.diagList.insertDiagonal(event, edgeToTheLeft)
        # helper(ej) <- event
        edgeToTheLeft.helper = event
        # Insert ei in edge tree and set helper(ei) to event
        event.edge.helper = event
        self.sweepStates.edgeTree.insertEdge(event.edge)

    def handleEndVertex(self,event):
        # if helper(ei−1) is a merge vertex
        predEdge = event.pred.edge
        if predEdge.helper.type == VertexType.MERGE:
            # Insert the diagonal connecting the event to helper(ei−1)
            # into the list of diagonals.
            self.sweepStates.diagList.insertDiagonal(event, predEdge)
        # Delete ei−1 from edge tree
        EdgeData.currentY = event.y
        self.sweepStates.edgeTree.removeEdge(predEdge)


    def triangulateMonotone(self, startEdge):
        vertices = self.getYMonotoneVertices(startEdge)
        if len(vertices) == 3:
            self.triangles.append( tuple(vertices) )
            return

        # Initialize an empty stack and push the first two vertices
        triangulationStack = vertices[:2]
        for i in range(2,len(vertices)):
            curr = vertices[i]
            # Ff the current vertex and the vertex on top of the stack are on different chains
            if triangulationStack[-1].chainType != curr.chainType:
                prevPopped = None
                # Pop all vertices from stack and create the triangles
                while True:
                    popped = triangulationStack.pop()
                    if prevPopped:
                        self.triangles.append( (curr, popped, prevPopped) ) 
                    prevPopped = popped
                    if len(triangulationStack) <= 1:
                        break
                v1 = triangulationStack.pop()
                if prevPopped:
                    self.triangles.append( (curr, v1, prevPopped) )
                triangulationStack.append(vertices[i - 1])
                triangulationStack.append(curr)
            else:  # Similar, when they are on the same chains
                prevPopped = triangulationStack.pop()
                popped = None
                while triangulationStack:
                    popped = triangulationStack.pop()
                    if self.connectionLiesInside(curr, prevPopped, popped):
                        self.triangles.append( (curr, popped, prevPopped) )
                        prevPopped = popped
                    else:
                        triangulationStack.append(popped)
                        popped = prevPopped
                        break
                triangulationStack.append(popped)
                triangulationStack.append(curr)


    def connectionLiesInside(self, pred, this, succ):
        cross = (this[0] - pred[0]) * (succ[1] - this[1]) - \
                (succ[0] - this[0]) * (this[1] - pred[1])
        return cross<0. if pred.chainType == ChainType.LEFT else cross>0.


    def getYMonotoneVertices(self,startEdge):
        # Get the edge with the highest vertex (note that the comparison
        # operator of Edge returns the maximum y as smallest)
        maxEdge = min( startEdge.iter(), key = lambda x: x.orig )
        maxEdge.triangulated = True
        maxEdge.orig.chainType = ChainType.START

        # Traversing left and right chain edges from there automatically
        # sorts their vertices. Classify the edges by ChainType.
        vertices = [maxEdge.orig]
        currSucc = maxEdge.succ
        currPred = maxEdge.pred
        while currSucc != currPred:
            if currSucc.orig < currPred.orig: # currSucc is higher
                currSucc.triangulated = True
                currSucc.orig.chainType = ChainType.LEFT
                vertices.append(currSucc.orig)
                currSucc = currSucc.succ
            else:                             # currPred is higher
                currPred.triangulated = True
                currPred.orig.chainType = ChainType.RIGHT
                vertices.append(currPred.orig)
                currPred = currPred.pred
        currSucc.orig.triangulated = True
        vertices.append(currSucc.orig)

        return vertices

