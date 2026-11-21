"""Propagate a prepared correction only within one original dimensional run."""
import math
from .item import Intersection, Section
from .source_placement import placement_signature


def source_width(section):
    value = getattr(section, 'sourceWidth', None)
    return section.width if value is None else value


def compatible(reference, section):
    """A real lane/width change at a minor junction is a propagation boundary.

    A SideLane's optional pixel correction can change its prepared width. Use
    the original resolved width, recorded before that correction, to distinguish
    an equal-width continuation from another source lane count or width tag.
    """
    return (getattr(reference, 'totalLanes', None) == getattr(section, 'totalLanes', None)
            and placement_signature(getattr(reference, 'tags', {}))
                == placement_signature(getattr(section, 'tags', {}))
            and (not placement_signature(getattr(reference, 'tags', {}))
                 or reference.tags.get('width:lanes') == section.tags.get('width:lanes'))
            and math.isclose(source_width(reference), source_width(section), rel_tol=0., abs_tol=1e-7))


class TransitionDimensionPropagation:
    def apply(self, transition, *, offset=False):
        reference = transition.succ
        if reference is None:
            return
        width, correction = reference.width, getattr(reference, 'offset', 0.)
        item, seen = reference, {id(reference)}
        while item.succ is not None:
            item = item.succ
            if id(item) in seen:
                raise ValueError('Dimension propagation encountered a cyclic original item chain')
            seen.add(id(item))
            if isinstance(item, Section):
                if not compatible(reference, item):
                    break
                item.width = width
                if offset:
                    item.offset = correction
                item.widthOrigin = ('Propagated after SideLane by Blosm preparation' if offset else
                                    'Propagated after SymLane by Blosm preparation')
            elif not (isinstance(item, Intersection) and item.isMinor):
                break
