from .HalfEdge import HalfEdge

class DiagonalsList():
    def __init__(self):
        self.diagonals = []

    def insertDiagonal(self,v1, helperEdge):
        v2 = helperEdge.helper
        assert v2 is not None

        v2Succ = self.findNextEdge(helperEdge, v2)
        v2Pred = v2Succ.pred
        v1Succ = v1.edge
        v1Pred = v1Succ.pred

        diag = HalfEdge()
        diagTwin = HalfEdge()

        diag.twin = diagTwin
        diagTwin.twin = diag

        diag.orig = v1
        diag.pred = v1Pred
        diag.succ = v2Succ

        v1Pred.succ = diag
        v2Succ.pred = diag

        diagTwin.orig = v2
        diagTwin.pred = v2Pred
        diagTwin.succ = v1Succ

        v2Pred.succ = diagTwin
        v1Succ.pred = diagTwin

        self.diagonals.append(diag)
        self.diagonals.append(diagTwin)

    def findNextEdge(self,helperEdge, v):
        if helperEdge.orig == v:
            return helperEdge
 
        curr = helperEdge.pred
        while curr.orig != v and curr != helperEdge:
            curr = curr.pred

        if curr == helperEdge:
            raise Exception('Vertex not reachable from edge.')
        else:
            return curr
