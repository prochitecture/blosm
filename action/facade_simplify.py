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
            vectorData = np.zeros((nEdges,4))    # shared, length, sine start, sine end, class
            vectors = [vector for vector in building.polygon.getVectors()]
            vectorData[:,:3] = [(vector.edge.hasSharedBldgVectors(), vector.edge.length, np.cross(vector.prev.unitVector, vector.unitVector)) for vector in vectors]
            vectorData[:,3] = np.roll(vectorData[:,2],-1)

            # detect curvy sequences
            curvyEdges = self.detectCurvyEdges(vectorData,vectors)

            # a primitive filter to avoid spiky edge detection for all buildings
            nLongEdges = np.where( vectorData[:,0]>=lengthThresh  )[0].shape[0]
            nShortEdges = np.where( vectorData[:,0]<lengthThresh )[0].shape[0]
            smallPatterns = None
            if (nLongEdges and nShortEdges > 2) or nShortEdges > 5:
                smallPatterns = self.detectSmallPatterns(vectorData,vectors)

            self.replaceSequences(curvyEdges, building.polygon)
            self.replaceSequences(smallPatterns, building.polygon)

    # Detects curvy sequences.
    # Returns:
    #   (firstEdge,lastEdge,patternClass) : Tuple of first edge of sequence and last edge of sequence and the pattern class.
    #                                       If firstEdge is lastEdge, the whole building polygon is a curvy sequence 
    #   None                              : No curvy sequence found 
    def detectCurvyEdges(self,vectorData,vectors):
        shared = vectorData[:,0]
        sineStart = vectorData[:,2]
        sineEnd = vectorData[:,3]
        lowAngles = (shared==0.0) & ( \
                        ((abs(sineStart)>sin_lo) & (abs(sineStart)<sin_me)) | \
                        ((abs(sineEnd)>sin_lo) & (abs(sineEnd)<sin_me)) 
                    )
        if not np.any(lowAngles):
            return None

        # estimate a length threshold
        curvyLengthThresh = np.mean(vectorData[lowAngles,1]) * curvyLengthFactor

        # pattern character sequence. edges with angle between 5° and 30° on either end 
        # and a length below <curvyLengthThresh> get a 'C', else a '0'
        sequence =  "".join( np.where( (vectorData[:,1]<curvyLengthThresh) & lowAngles,'C','0') )

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

    # Detects small patterns (rectangular and triangular).
    # Returns:
    #   (firstEdge,lastEdge,patternClass) : Tuple of first edge of sequence and last edge of sequence and the pattern class.
    #                                       If firstEdge is lastEdge, the whole building polygon is a curvy sequence 
    #   None                              : No patterns found 
    def detectSmallPatterns(self,vectorData,vectors):
        # Shared edges:
        #       'X': edges shared
        # Long edges (>=lengthThresh):
        #       'L': sharp left at both ends
        #       'R': sharp right at both ends
        #       'O': other long edge
        # Short edges:
        #       'l': medium left at both ends
        #       'r': medium right at both ends
        #       '>': alternate medium angle, starting to right
        #       '<': alternate medium angle, starting to left
        #       'o': other short edge
        shared = vectorData[:,0]
        length = vectorData[:,1]
        sineStart = vectorData[:,2]
        sineEnd = vectorData[:,3]
        sequence =  "".join( np.where( (shared!=0.0), 'X', \
                                np.where( (length>=lengthThresh), \
                                ( np.where( ((sineStart>sin_hi) & (sineEnd>sin_hi)), 'L', \
                                  np.where( ((sineStart<-sin_hi) & (sineEnd<-sin_hi)), 'R', 'O') )  \
                                ), \
                                ( np.where( ((sineStart>sin_me) & (sineEnd>sin_me)), 'l', \
                                  np.where( ((sineStart<-sin_me) & (sineEnd<-sin_me)), 'r',  \
                                  np.where( (sineStart<-sin_me)&(sineEnd>sin_me), '>', \
                                  np.where( (sineStart>sin_me)&(sineEnd<-sin_me), '<', 'o') ) ) )\
                                ) ) )
                            )

        N = len(sequence)
        sequence = sequence+sequence # allow cyclic pattern
        smallPatterns = []

        # convex rectangular pattern
        pattern = re.compile(r"(>[L|l]<)")
        matches = [r for r in pattern.finditer(sequence)]
        if matches:
            for spikyySeg in matches:
                s = spikyySeg.span()
                if s[0] < N and s[0] >= 0:
                    smallPatterns.append( ( vectors[s[0]], vectors[(s[1]-1)%N], PatternClass.rectangular ) )
        sequence = re.sub(pattern, lambda m: ('#' * len(m.group())), sequence)

        # triangular pattern
        pattern = re.compile(r">(>|<|l){1,}")
        matches = [r for r in pattern.finditer(sequence)]
        if matches:
            for spikyySeg in matches:
                s = spikyySeg.span()
                if s[0] < N and s[0] >= 0:
                    smallPatterns.append( ( vectors[s[0]], vectors[(s[1]-1)%N], PatternClass.triangular ) )
        sequence = re.sub(pattern, lambda m: ('#' * len(m.group())), sequence)

        # concave rectangular pattern
        pattern = re.compile(r"(<[R,r]>)")
        matches = [r for r in pattern.finditer(sequence)]
        if matches:
            for spikyySeg in matches:
                s = spikyySeg.span()
                if s[0] < N and s[0] >= 0:
                    smallPatterns.append( ( vectors[s[0]], vectors[(s[1]-1)%N], PatternClass.rectangular ) )
        sequence = re.sub(pattern, lambda m: ('#' * len(m.group())), sequence)

        # triangular pattern
        pattern = re.compile(r"<(>|<|r){1,}")
        matches = [r for r in pattern.finditer(sequence)]
        if matches:
            for spikyySeg in matches:
                s = spikyySeg.span()
                if s[0] < N and s[0] >= 0:
                    smallPatterns.append( ( vectors[s[0]], vectors[(s[1]-1)%N], PatternClass.triangular ) )
        sequence = re.sub(pattern, lambda m: ('#' * len(m.group())), sequence)

        if smallPatterns:
            return smallPatterns
        else:
            return None
