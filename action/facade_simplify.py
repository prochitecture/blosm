from math import atan2
import numpy as np
# import matplotlib.pyplot as plt
import re
from defs.facade_classification import *
# from building import BldgEdge, BldgVector


class FacadeSimplification:

    def __init__(self):
        pass

    def replaceSequences(self, sequences, polygon):
        if sequences:
            for sequence in sequences:
                vector1 = sequence[0]
                vector2 = sequence[1]
                cl = sequence[2]
                head = vector1
                head.edge.pattern = cl
                cur = head.next
                while cur is not vector2.next:
                    cur.edge.pattern = cl
                    cur = cur.next
                # proxyEdge = BldgEdge(vector1.prev.edge.id1, vector1.prev.edge.v1, vector2.edge.id2, vector2.edge.v2)
                # proxyVector = BldgVector(proxyEdge,True,polygon)
                # proxyVector.prev = vector1.prev
                # proxyVector.next = vector2.next
                # vector1.prev.next = proxyVector
                # vector2.next.prev = proxyVector 
                # proxyEdge.addVector(proxyVector)   

    def do(self, manager):
        buildings = manager.buildings
        
        # Create an instance of <util.polygon.Polygon> for each <building>,
        # remove straight angles for them and calculate the total number of vertices
        for building in buildings:
            if not building.polygon:
                building.initPolygon(manager)
        
        for building in buildings:
            building.polygon.processStraightAngles(manager)
        
        for building in buildings:
            building.polygon.processStraightAnglesExtra(manager)

        for building in buildings:

            # prepare numpy matrix with data used for pattern detection
            nEdges = building.polygon.numEdges
            vectorData = np.zeros((nEdges,4))    # length, sine left, sine right, class
            vectors = [vector for vector in building.polygon.getVectors()]
            vectorData[:,:2] = [(vector.edge.length,np.cross(vector.prev.unitVector, vector.unitVector)) for vector in vectors]
            vectorData[:,2] = np.roll(vectorData[:,1],-1)
            vectorData[:,3] = PatternClass.unclassified

            # detect curvy sequences
            curvyEdges = self.detectCurvyEdges(vectorData,vectors)

            # a primitive filter to avoid spiky edge detection for all buildings
            nLongEdges = np.where( vectorData[:,0]>=lengthThresh  )[0].shape[0]
            nShortEdges = np.where( vectorData[:,0]<lengthThresh )[0].shape[0]
            smallPatterns = None
            if (nLongEdges and nShortEdges > nLongEdges) or nShortEdges > 10:
                smallPatterns = self.detectSmallPatterns(vectorData,vectors)

            self.replaceSequences(curvyEdges, building.polygon)
            self.replaceSequences(smallPatterns, building.polygon)
            test =1 

    # Detects curvy sequences.
    # Returns:
    #   (firstEdge,lastEdge,patternClass) : Tuple of first edge of sequence and last edge of sequence and the pattern class.
    #                                       If firstEdge is lastEdge, the whole building polygon is a curvy sequence 
    #   None                              : No curvy sequence found 
    def detectCurvyEdges(self,vectorData,vectors):
        lowAngles = ((abs(vectorData[:,1])>sin_lo) & (abs(vectorData[:,1])<sin_hi)) | \
                    ((abs(vectorData[:,2])>sin_lo) & (abs(vectorData[:,2])<sin_hi)) 
        if not np.any(lowAngles):
            return None

        # estimate a length threshold
        curvyLengthThresh = np.mean(vectorData[lowAngles,0]) * curvyLengthFactor

        # pattern character sequence. edges with angle between 5° and 30° on either end 
        # and a length below <curvyLengthThresh> get a 'C', else a '0'
        sequence =  "".join( np.where( (vectorData[:,0]<curvyLengthThresh) & lowAngles,'C','0') )

        # a sequence of four or more 'C' matches as curvy sequence
        pattern = re.compile(r"(C){4,}")
        matches = [c for c in pattern.finditer(sequence+sequence)]  # adjacent sequence for circularity
        curvyEdges = []
        if matches:
            N = len(sequence)
            for curvySeg in matches:
                s = curvySeg.span()
                if s[0] < N and s[0] >= 0:
                    curvyEdges.append( ( vectors[s[0]], vectors[(s[1]-1)%N], PatternClass.curvy ) )

        return curvyEdges

    # Detects small patterns (spikes and balconies).
    # Returns:
    #   (firstEdge,lastEdge,patternClass) : Tuple of first edge of sequence and last edge of sequence and the pattern class.
    #                                       If firstEdge is lastEdge, the whole building polygon is a curvy sequence 
    #   None                              : No patterns found 
    def detectSmallPatterns(self,vectorData,vectors):
        # characters in sequence:
        # 'S': short edge, turning at left at least at one end
        # 'B': long edge, turning at left at both ends
        # '0': all other edges
        sequence =  "".join( np.where( vectorData[:,0]<lengthThresh , \
                             np.where( (vectorData[:,1]>sin_hi) | (vectorData[:,2]>sin_hi), 'S', '0') , \
                             np.where( (vectorData[:,1]>sin_hi) & (vectorData[:,2]>sin_hi), 'B', '0') ) )

        smallPatterns = []
        pattern = re.compile(r"(S){2,}")
        group_repl = lambda m: ('#' * len(m.group()))
        matches = [r for r in pattern.finditer(sequence+sequence)]
        if matches:
            N = len(sequence)
            for spikyySeg in matches:
                s = spikyySeg.span()
                if s[0] < N and s[0] >= 0:
                    smallPatterns.append( ( vectors[s[0]], vectors[(s[1]-1)%N], PatternClass.spike ) )
                    sequence = re.sub(pattern, group_repl, sequence)

        pattern = re.compile(r"(SBS){1,}")
        matches = [r for r in pattern.finditer(sequence+sequence)]
        if matches:
            N = len(sequence)
            for balkonySeg in matches:
                s = balkonySeg.span()
                if s[0] < N and s[0] >= 0:
                    smallPatterns.append( ( vectors[s[0]], vectors[(s[1]-1)%N], PatternClass.balcony ) )

        if smallPatterns:
            return smallPatterns
        else:
            return None