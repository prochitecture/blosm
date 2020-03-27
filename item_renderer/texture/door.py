

_materialTemplateFilename = "building_material_templates.blend"


class Door:
    """
    The Door renderer is the special case of the <item_renderer.level.Level> when
    a door in the only element in the level markup
    
    A mixin class for Door texture based item renderers
    """
    
    def __init__(self):
        self.facadeMaterialTemplateFilename = _materialTemplateFilename
        # do we need to initialize <self.facadePatternInfo>
        self.initFacadePatternInfo = False
        self.facadePatternInfo = dict(Door=1)