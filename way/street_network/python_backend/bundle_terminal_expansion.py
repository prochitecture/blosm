"""Close complete consumed Bundle members through original terminal owners."""
from .model import End


def missing_terminals(context, bundle, consumed, members, available):
    """Return only proven missing opposite original physical terminals.

    Consuming an entire short member cannot leave an unsupported partial paired
    body.  Acquire the opposite original terminal owners before constructing a
    source-cycle exterior.  Existing member/pass budgets and all final physical
    coverage/separation guards still apply to that expanded group.
    """
    if (len(bundle.members) != 2 or not all(m.complete for m in bundle.members)
            or set(bundle.sections) <= consumed):
        return ()
    physical = {j.key: j for j in context.network.intersections if len(j.ends) > 2}
    incident = {end: j.key for j in physical.values() for end in j.ends}
    terminal_pairs = []
    for member in bundle.members:
        first, last = member.path[0], member.path[-1]
        terminal_pairs.append((End(first.section_id, 0 if first.forward else 1),
                               End(last.section_id, 1 if last.forward else 0)))
    complete = [role for role, member in enumerate(bundle.members)
                if {v.section_id for v in member.path} <= consumed
                and all(incident.get(e) in members for e in terminal_pairs[role])]
    if len(complete) != 1:
        return ()
    opposite = 1-complete[0]
    terminals = terminal_pairs[opposite]
    owners = tuple(incident.get(e) for e in terminals)
    if any(owner is None or owner not in available for owner in owners):
        raise ValueError('Whole-member Bundle expansion needs available original physical terminal owners')
    sections = {v.section_id for v in bundle.members[opposite].path}
    for sid in sections:
        source = context.network.section_index[sid]
        if source.scope_ends or source.selection_ends:
            raise ValueError('Whole-member Bundle expansion cannot cross a source-scope or selection boundary')
    # A long opposite chain with another unowned physical attachment needs its
    # own bounded ownership preparation. Do not jump over that original site.
    for end, owner in incident.items():
        if end.section_id in sections and owner not in members and end not in terminals:
            raise ValueError('Whole-member Bundle expansion cannot skip an unowned original attachment')
    return tuple(sorted(set(owners)-set(members)))
