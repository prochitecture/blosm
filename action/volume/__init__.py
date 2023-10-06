from .. import Action
from item.footprint import Footprint
import parse

# import roof generators
from .roof_flat import RoofFlat
from .roof_flat_multi import RoofFlatMulti
from .roof_generatrix import RoofGeneratrix
from .roof_profile import RoofProfile, roofDataGabled, roofDataRound, roofDataGambrel, roofDataSaltbox
from .roof_hipped import RoofHipped
from .roof_hipped_multi import RoofHippedMulti


class Volume(Action):
    """
    This action creates a building volume out of its footprint
    """
    
    defaultRoofShape = "flat"
    
    def __init__(self, manager, itemStore, itemRenderers):
        super().__init__(manager.app, manager.data, itemStore)
        self.manager = manager
        if itemRenderers:
            self.setVolumeGenerators(manager.data, itemRenderers)
    
    def preprocess(self, buildingsP):
        # <buildingsP> means "buildings from the parser"
        return
    
    def setRenderer(self, renderer):
        """
        Sets a renderer.
        For the 3D mode it's typically a Facade renderer.
        For the 2D mode it's a dedicated 3D renderer that creates footprints only
        """
        self.renderer = renderer
    
    def prepareFootprint(self, footprint, building, buildingStyle):
        footprint.buildingStyle = buildingStyle
        
        if footprint.styleBlock:
            # the footprint has been generated
            footprint.calculateFootprint()
        else:
            # the footprint is defined in the external data (e.g. OpenStreetMap)
            footprint.calculateStyling()
        # Check if have one or more footprints are defined in the markup definition,
        # it actually means, that those footprints are to be generated
        styleBlocks = footprint.styleBlock.styleBlocks.get("Footprint")
        if styleBlocks:
            for styleBlock in styleBlocks:
                _footprint = Footprint(None, building, styleBlock)
                _footprint.parent = footprint
                _footprint.buildingStyle = buildingStyle
                self.itemStore.add(_footprint)
    
    def do(self, building, buildingStyle, globalRenderer):
        itemStore = self.itemStore
        while itemStore.hasItems(Footprint):
            footprint = itemStore.getItem(Footprint)
            self.prepareFootprint(footprint, building, buildingStyle)
            
            element = footprint.element
            if element.t is parse.multipolygon:
                # check if the multipolygon has holes
                if element.hasInner():
                    if footprint.getStyleBlockAttr("roofShape") in ("hipped", "gabled"):
                        self.volumeGeneratorMultiHipped.do(footprint)
                    else:
                        self.volumeGeneratorMultiFlat.do(footprint)
                else:
                    # That's a quite rare case
                    # We treat each polygon of the multipolygon as a single polygon
                    # <footprint> won't be used in this case, a new footprint will be created
                    # for each polygon of the multipolygon
                    
                    # A note.
                    # If <building.footprint> is set, for example, in a <building.area()> call,
                    # we don't need to care about it in any way. Since <building.footprint>
                    # won't be used. A new footprint will be created for each polygon of
                    # the multipolygon.
                    
                    ls = element.ls
                    for _l in ls:
                        footprint = Footprint(None, building)
                        footprint.element = element
                        footprint.createPolygon(_l, self.manager)
                        self.prepareFootprint(footprint, building, buildingStyle)
                        self.generateVolume(footprint)
                    element.ls = ls
            else:
                self.generateVolume(footprint)
    
    def generateVolume(self, footprint):
        volumeGenerator = self.volumeGenerators.get(
            footprint.getStyleBlockAttr("roofShape"),
            self.volumeGenerators[Volume.defaultRoofShape]
        )
        volumeGenerator.do(footprint)
    
    def setVolumeGenerators(self, data, itemRenderers):
        self.volumeGenerators = {
            'flat': RoofFlat("RoofFlat", data, self, itemRenderers),
            'gabled': RoofProfile(roofDataGabled, data, self, itemRenderers),
            'pyramidal': RoofGeneratrix("RoofPyramidal", data, self, itemRenderers),
            #'skillion': RoofSkillion(),
            'hipped': RoofHipped(data, self, itemRenderers),
            'dome': RoofGeneratrix("RoofDome", data, self, itemRenderers),
            'half-dome': RoofGeneratrix("RoofHalfDome", data, self, itemRenderers),
            'onion': RoofGeneratrix("RoofOnion", data, self, itemRenderers),
            'round': RoofProfile(roofDataRound, data, self, itemRenderers),
            #'half-hipped': RoofHalfHipped(),
            'gambrel': RoofProfile(roofDataGambrel, data, self, itemRenderers),
            'saltbox': RoofProfile(roofDataSaltbox, data, self, itemRenderers)
            #'mansard': RoofMansard()
        }
        self.volumeGeneratorMultiFlat = RoofFlatMulti(data, self, itemRenderers)
        self.volumeGeneratorMultiHipped = RoofHippedMulti(data, self, itemRenderers)