from ..building.manager import BuildingManager


class RealisticBuildingManager(BuildingManager):
    
    def process(self):
        for _ in self.process_steps():
            pass

    def process_steps(self, session=None):
        yield from super().process_steps(session)
        for action in self.renderer.actions:
            action.preprocess(self.buildings)
            yield None
    
    def render(self):
        renderer = self.renderer
        
        for building in self.buildings:
            renderer.render(building, self.data)
        
        if not self.app.renderAfterExtrude:
            for building in self.buildings:
                renderer.renderExtrudedVolumes(building, self.data)