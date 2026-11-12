from ..building.manager import BuildingManager
from . import Building


class RealisticBuildingManager(BuildingManager):
    
    def process(self):
        for _ in self.process_steps():
            pass

    def process_steps(self, session=None):
        yield from super().process_steps(session)
        for action in Building.actions:
            action.preprocess(self.buildings)
            yield None