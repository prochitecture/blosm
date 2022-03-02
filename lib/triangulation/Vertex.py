from mathutils import Vector

class VertexType():
    START = 0
    SPLIT = 1
    END = 2
    MERGE = 3
    REGULAR = 4

class ChainType():
    START = 0,
    LEFT = 1,
    RIGHT = 2

class Vertex(Vector):
    ID = 0
    def __init__(self,xy):
        self.vType = None
        self.pred = None
        self.edge = None
        self.chainType = None
        self.id = Vertex.ID
        Vertex.ID += 1

    def __lt__(self,other):
        if self[1] == other[1]:
            return self[0] < other[0]
        else:
            return self[1] > other[1]

    def __gt__(self,other):
        if self[1] == other[1]:
            return self[0] > other[0]
        else:
            return self[1] < other[1]

    def __eq__(self,other):
        return self[0] == other[0] and self[1] == other[1]

    def __hash__(self):
        return self.id
