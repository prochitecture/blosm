import heapq
from itertools import tee, cycle, islice
from .EdgeTree import EdgeTree
from .Vertex import VertexType
from .HalfEdge import HalfEdge
from .DiagonalsList import DiagonalsList

def iterCircular(iterable):
    preds, this, succs = tee(iterable, 3)
    preds = islice(cycle(preds), len(iterable) - 1, None)
    succs = islice(cycle(succs), 1, None)
    return zip(preds, this, succs)

# Prepares the doubly linked lists of polygon and holes and computes
# the sweep events and their types.
# Holds the status data of the monotonization of a polygon with holes.
class MonotoneStates():
    def __init__(self):
        self.edgeTree = EdgeTree()
        self.diagList = DiagonalsList()
        self.sweepEvents = []
        self.splitOrMergeCount = 0

    def initialize(self, polygon, holes):
        """
            polygon: List of vertices of class Vertex in counter-clockwise order
            holes:   List of holes, each hole a list of vertices of class in clockwise order.
                     The holes must not intersect each other and have to be completely
                     inside the polygon.
            return:  The number of vertices classified as type MERGE or SPLIT.
        """
        for predV, thisV, succV in iterCircular(polygon):
            self.splitOrMergeCount += self.getVertexType(thisV,succV,predV)
            heapq.heappush(self.sweepEvents,thisV)
        HalfEdge.createCircularLinkedList(polygon)

        for hole in holes:
            for predV, thisV, succV in iterCircular(hole):
                self.splitOrMergeCount += self.getVertexType(thisV,succV,predV)
                heapq.heappush(self.sweepEvents,thisV)
            HalfEdge.createCircularLinkedList(hole)

        return self.splitOrMergeCount

    def getVertexType(self,thisV,succV,predV):
        isSuccVbelow = thisV[1] > succV[1] or (thisV[1] == succV[1] and thisV[0] < succV[0])
        isPredVbelow = thisV[1] > predV[1] or (thisV[1] == predV[1] and thisV[0] < predV[0])
        # isReflex is computed for counter-clockwise order
        isReflex = (thisV[0] - predV[0]) * (succV[1] - thisV[1]) - \
                   (succV[0] - thisV[0]) * (thisV[1] - predV[1]) < 0.
        if isSuccVbelow and isPredVbelow:   # start or split
            if isReflex:
                thisV.type = VertexType.SPLIT
                return 1
            else:
                thisV.type = VertexType.START
                return 0
        elif not isSuccVbelow and not isPredVbelow:   # end oder merge
            if isReflex:
                thisV.type = VertexType.MERGE
                return 1
            else:
                thisV.type = VertexType.END
                return 0
        else:
            thisV.type = VertexType.REGULAR
            return 0

