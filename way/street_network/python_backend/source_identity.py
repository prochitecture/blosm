"""Canonical numerical resources across equal immutable input reconstructions."""
from .model import Key, Section


def matches_current_section(context, section):
    """Keep Path provenance strict while allowing an equal reread Section.

    UpdateCoordinator deliberately reuses native Paths when their source digest
    is unchanged. A new scene read can allocate an equal Section record; that
    record's Python address is not its persistent source identity. The consumer
    must still use the context's canonical native resource, and every authored
    Section value must match the currently supplied input.
    """
    if type(section) is not Section:
        return False
    current = context.network.section_index.get(section.id)
    path = context.products.get(Key('path', section.id))
    return (type(current) is Section and path is not None
            and path.section is section and section == current)
