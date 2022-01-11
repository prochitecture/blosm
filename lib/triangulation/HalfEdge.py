from math import inf

# Represents an directed edge connecting a vertex to another vertex. Each
# HalfEdge has a 'twin', which is the same edge but in the opposite direction.
# The class is a node of a doubly linked edge list.
class HalfEdge():
    def __init__(self):
        self.orig = None    # The start vertex of this edge
        self.twin = None    # The edge in the opposite direction
        self.succ = None    # The next edge (successor)
        self.pred = None    # The previous edge (predecessor)
        self.triangulated = False
        self.xSlopeCached = None
        self.helper = None

    def iter(self):
        currEdge = self
        while True:
            yield currEdge
            currEdge = currEdge.succ
            if currEdge == self:
                break

    def xSlope(self):
        dy = (self.succ.orig[1] - self.orig[1])
        if not dy:
            return inf  # horizontal edge
        return (self.succ.orig[0] - self.orig[0]) / dy

    def intersectX(self, y):
        """
            Calculates the y-coordinate for a given x-coordinate on a infinite
            line through this edge.
        """
        if not self.xSlopeCached:
            self.xSlopeCached = self.xSlope()
        if self.xSlopeCached == inf:
            return self.orig[0]    # should never happen
        else:
            return self.orig[0] + (y - self.orig[1]) * self.xSlopeCached 

    # create a circular doubly linked list of nodes HalfEdge from a polygon,
    # given as a list of vertices of class Vertex.
    @staticmethod
    def createCircularLinkedList(polygon):
        # create first node
        edge = HalfEdge()
        twin = HalfEdge()
        edge.orig = polygon[0]
        edge.succ = edge
        edge.pred = edge
        edge.twin = twin
        twin.twin = edge

        head = edge
        tail = edge

        # insert edges at end
        for vertex in polygon[1:]:
            edge = HalfEdge()
            twin = HalfEdge()
            edge.orig = vertex
            edge.succ = head
            edge.pred = tail
            edge.twin = twin
            tail.succ = edge
            head.pred = edge
            tail = edge

        # initiate and link twins
        curr = head
        while True:
            curr.orig.edge = curr
            curr.orig.pred = curr.pred.orig
            curr.twin.orig = curr.succ.orig
            curr.twin.pred = curr.succ.twin
            curr.twin.succ = curr.pred.twin
            curr = curr.succ
            if curr == head:
                break

        return head

# Helper class used for the search for the edge left of a vertex in EdgeTree.
# See there.
class SearchEdge(HalfEdge):
    def __init__(self):
        self.x = None

    def intersectX(self, y):
        return self.x