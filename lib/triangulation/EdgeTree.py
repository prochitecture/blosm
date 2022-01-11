from .RBTree import RedBlackTree
from .HalfEdge import SearchEdge

# The class EdgeData is the data class stored in the red-black tree. The task
# ot this class is to provide a comparison operator, so that the edges remain
# sorted from left to right at different y-values of the sweep line. The
# edges of simple polygons with holes never intersect each other. When we insert
# a new edge, the order of the edges already in the tree remains unchanged. The
# comparison is therefore made by computing the x-value of the intersections at
# the current y-value currentY of the swepp line. This is required log(n) times
# for every insertion, deletion and search action, when n is the number of edges
# stred in the tree.
class EdgeData():
    currentY = None
    def __init__(self,edge):
        self.edge = edge

    def __lt__(self,other):
        sX = self.edge.intersectX(EdgeData.currentY)
        oX = other.edge.intersectX(EdgeData.currentY)
        return sX < oX

    def __le__(self,other):
        sX = self.edge.intersectX(EdgeData.currentY)
        oX = other.edge.intersectX(EdgeData.currentY)
        return sX <= oX

    def __gt__(self,other):
        sX = self.edge.intersectX(EdgeData.currentY)
        oX = other.edge.intersectX(EdgeData.currentY)
        return sX > oX

    def __ge__(self,other):
        sX = self.edge.intersectX(EdgeData.currentY)
        oX = other.edge.intersectX(EdgeData.currentY)
        return sX >= oX

# To store the active edges with a complexity of O(log n), a red-black binary search tree
# is used. Its task is to sort the active edges from left to right and, beside insertion
# and deletion, to find the edge immediately left of a vertex. The required comparison of
# of the edges requires a tricky data class EdgeData described above.
class EdgeTree():
    def __init__(self):
        self.sortedEdges = RedBlackTree()
        self.currentY = None
        self.searchEdge = SearchEdge()

    def insertEdge(self,edge):
        sEdge = EdgeData(edge)
        self.sortedEdges.add(sEdge)

    def removeEdge(self,edge):
        sEdge = EdgeData(edge)
        self.sortedEdges.remove(sEdge)

    def findToLeft(self,vertex):
        self.searchEdge.x = vertex[0]
        sEdge = EdgeData(self.searchEdge)
        lEdge = self.sortedEdges.floor(sEdge)
        return lEdge.edge
