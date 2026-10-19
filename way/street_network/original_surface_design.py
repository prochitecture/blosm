"""Replaceable physical design policy for original-span surface authoring.

Topology preparation and GN construction are independent of these choices.
The default policy retains the normal import's category radii, common datum
and profile rise. An authored-input adapter can supply existing per-port and
per-Section choices without replacing the shared surface solver.
"""
from .physical_junctions import PreparedPhysicalJunctions
from ..way_properties import wayCategoryProps


class OriginalSurfaceDesign:
    def __init__(self, design, *, retain_unresolved_profiles=False):
        self.design = design
        self.retain_unresolved_profiles = retain_unresolved_profiles

    def ports(self, inventory, plan):
        return PreparedPhysicalJunctions(None, None, wayCategoryProps,
            rise=self.design.profile_rise, cell_size=self.design.cell_size
        ).prepare_port_designs(inventory, plan, retain_unresolved=self.retain_unresolved_profiles)

    def span(self, inventory, span):
        """Return physical (profile rise, datum) for an independent run."""
        return self.design.profile_rise, 0.

    def interpolation(self, section):
        return self.design.interpolation
