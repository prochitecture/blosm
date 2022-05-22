from lib.CompGeom.PolyLine import PolyLine

roadTypes = {
    "motorway" : 'driving',
    "motorway_link" : 'driving',
    "trunk" : 'driving',
    "trunk_link" : 'driving',
    "primary" : 'driving',
    "primary_link" : 'driving',
    "secondary" : 'driving',
    "secondary_link" : 'driving',
    "tertiary" : 'driving',
    "tertiary_link" : 'driving',
    "unclassified" : 'default',
    "residential" : 'driving',
    "living_street" : 'driving',
    "service" : 'driving',
    "pedestrian" : 'walking',
    "track" : 'walking',
    "escape" : 'driving',
    "raceway" : 'driving',
    "steps" : 'walking',
    "footway" : 'walking',
    "path" : 'walking',
    "cycleway" : 'cycling',
    "bridleway" : 'walking'
}

def getRoadWidth(tags):
    rType = tags.get('highway',None)
    if rType:
        rType = roadTypes[rType]
        if rType == 'driving':
            lanes = tags.get('lanes',None)
            if lanes:
                nrOfLanes = int(lanes)
            else:
                nrOfLanes = 1
            if not tags.get('oneway','no') == 'yes':
                nrOfLanes *= 2
            return nrOfLanes * 2.5#2.5
        elif rType == 'cycling':
            return 1.8
        elif rType == 'walking':
            return 3.0
        else:
            return 1.0  # ????
    else:
        return 1.0

    pass

class WaySection():
    ID = 0
    def __init__(self,net_section):
        self.originalSection = net_section
        self.leftWidth = getRoadWidth(net_section.tags)/2.
        self.rightWidth = getRoadWidth(net_section.tags)/2.
        self.polyline = PolyLine(net_section.path)
        self.isLoop = net_section.s == net_section.t
        self._sV = None  # vector of first segment
        self._tV = None  # vector of last segment
        self.trimS = 0.    # trim factor for start
        self.trimT = 0.    # trim factor for starget
        self.id = WaySection.ID
        WaySection.ID += 1

    @property
    def sV(self):
        if self._sV is None:
            self._sV = self.polyline.verts[1] - self.polyline.verts[0]
        return self._sV

    @property
    def tV(self):
        if self._tV is None:
            self._tV = self.polyline.verts[-2] - self.polyline.verts[-1]
        return self._tV

    def fwd(self):
        return [v for v in self.polyline]

    def rev(self):
        return [v for v in self.polyline.reversed()]
