from math import atan2
import numpy as np
import matplotlib.pyplot as plt
import re

class PatternClass:
    unclassified = 0
    curvy = 1
    spike = 2


sin05 = abs(np.sin(np.pi/180.*5.))
sin30 = abs(np.sin(np.pi/180.*30))
tolFact = 1.5
lengthThresh = 5.

class FacadePatterns:

    def __init__(self):
        pass
    
    def do(self, manager):
        for building in manager.buildings:
            nEdges = building.polygon.numEdges
            buildData = np.zeros((nEdges,8))    # v1[0], v1[1], v2[0], v2[1], length, cross-, cross+, patternType 
            indx = 0
            for vector in building.polygon.getVectors():
                edge, v1, v2 = vector.edge, vector.v1, vector.v2
                buildData[indx,0:2] = v1
                buildData[indx,2:4] = v2
                buildData[indx,4] = edge.length
                buildData[indx,5] = np.cross(vector.prev.unitVector, vector.unitVector)
                buildData[indx,6] = np.cross(vector.unitVector, vector.next.unitVector)
                indx += 1
            buildData[:,7] = PatternClass.unclassified

            self.detectCurvyEdges(buildData)

            longEdges = np.where( (buildData[:,4]>=lengthThresh)&( buildData[:,7] == PatternClass.unclassified) )[0].shape[0]
            shortEdges = np.where( (buildData[:,4]<lengthThresh)&( buildData[:,7] == PatternClass.unclassified) )[0].shape[0]
            if (longEdges and shortEdges > longEdges) or shortEdges > 10:
                # print(longEdges,shortEdges)
                self.detectSpikes(buildData,lengthThresh)

            indx = 0
            for vector in building.polygon.getVectors():
                vector.edge.pattern = buildData[indx,7]
                indx += 1

            # for indx in range(nEdges):
            #     color = 'red' if buildData[indx,7] == PatternClass.curvy else ('blue' if buildData[indx,7] == PatternClass.spike else 'black' )
            #     linewidth = 3 if buildData[indx,7] == PatternClass.curvy else (2 if buildData[indx,7] == PatternClass.spike else 0.5 )
            #     plt.plot( [buildData[indx,0],buildData[indx,2]], [buildData[indx,1],buildData[indx,3]], linewidth = linewidth, color = color )
            #     plt.plot(buildData[indx,0],buildData[indx,1],'k.', markersize=1)
            # # else:
            #     for vector in building.polygon.getVectors():
            #         edge, v1, v2 = vector.edge, vector.v1, vector.v2
            #         plt.plot( (v1[0], v2[0]), (v1[1], v2[1]), linewidth = 0.5, color = 'black' )
            #         plt.plot(v1[0], v1[1], 'k.', markersize=2.)

        # plt.gca().axis('equal')
        # plt.show()

    def detectCurvyEdges(self,buildData):
        sequence = "".join( np.where( ((abs(buildData[:,5])>sin05) & (abs(buildData[:,5])<sin30)) | \
                                      ((abs(buildData[:,6])>sin05) & (abs(buildData[:,6])<sin30)),'C','0') )
        pattern = re.compile(r"(C){4,}")
        matches = [r for r in pattern.finditer(sequence+sequence)]

        if matches:
            N = len(sequence)
            curvySegmentsIndx = []
            for curvySeg in matches:
                s = curvySeg.span()[0]
                if s < N and s >= 0:
                    curvySegmentsIndx.extend( [ i%N for i in range(*curvySeg.span())] )

            meanCurvyEdgesLength = np.mean(buildData[curvySegmentsIndx,4])
            for indx in curvySegmentsIndx:
                if buildData[indx,4]<tolFact*meanCurvyEdgesLength:
                    buildData[indx,7] = PatternClass.curvy

    # def detectSpikes(self,buildData):
    #     sequence = "".join( np.where( ((buildData[:,5]>sin30) | (buildData[:,6]>sin30)),'L','0') )
    #     pattern = re.compile(r"(L){2,}")
    #     matches = [r for r in pattern.finditer(sequence+sequence)]

    #     if matches:
    #         N = len(sequence)
    #         curvySegmentsIndx = []
    #         for curvySeg in matches:
    #             s = curvySeg.span()[0]
    #             if s < N and s >= 0:
    #                 curvySegmentsIndx.extend( [ i%N for i in range(*curvySeg.span())] )

    #         meanCurvyEdgesLength = np.mean(buildData[curvySegmentsIndx,4])
    #         for indx in curvySegmentsIndx:
    #             if buildData[indx,4]<5 and buildData[indx,7] == PatternClass.unclassified:
    #                 buildData[indx,7] = PatternClass.spike

    def detectSpikes(self,buildData,lengthThresh):
        sequence =  "".join( np.where( buildData[:,4]<lengthThresh , \
                             np.where( (buildData[:,5]>sin30) | (buildData[:,6]>sin30), 'L', '0') , \
                             np.where( (buildData[:,5]>sin30) & (buildData[:,6]>sin30), 'T', 'Q') ) )

        spikySegmentsIndx = []

        pattern = re.compile(r"(L){2,}")
        matches = [r for r in pattern.finditer(sequence+sequence)]
        if matches:
            N = len(sequence)
            for spikyySeg in matches:
                s = spikyySeg.span()[0]
                if s < N and s >= 0:
                    spikySegmentsIndx.extend( [ i%N for i in range(*spikyySeg.span())] )

        pattern = re.compile(r"(LTL){1,}")
        matches = [r for r in pattern.finditer(sequence+sequence)]
        if matches:
            N = len(sequence)
            for balkonySeg in matches:
                s = balkonySeg.span()[0]
                if s < N and s >= 0:
                    spikySegmentsIndx.extend( [ i%N for i in range(*balkonySeg.span())] )

        spikySegmentsIndx = list(dict.fromkeys(spikySegmentsIndx))

        for indx in spikySegmentsIndx:
            if buildData[indx,7] == PatternClass.unclassified:
                buildData[indx,7] = PatternClass.spike
