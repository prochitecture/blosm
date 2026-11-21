"""Acquire a paired partial terminal through its original neighboring route."""
from .model import End
from .routes import RouteTopology


def _local_source_closure(context, topology, sites, owner, members, available, width):
    """At most four real route hops inside a finite terminal contact budget."""
    from .model import Key
    radius = max(sites[key].radius for key in (owner, *members) if key in sites)
    maximum = max(4.*width, 2.*radius+width)
    pending = [(owner, (owner,), 0.)]
    candidates = []
    shortest, visited = None, 0
    lengths = {}
    while pending:
        current, chain, distance = pending.pop(0)
        visited += 1
        if visited > 64:raise ValueError('Paired partial terminal exceeds its original-route search budget')
        if len(chain) > 4 or (shortest is not None and len(chain) > shortest):continue
        for entry in sorted(sites[current].ends):
            route = topology.trace(entry); neighbor = route.terminal_site
            if (neighbor is None or neighbor in chain or neighbor not in sites
                    or neighbor not in available or route.terminal_scope or route.terminal_selection):continue
            total = distance
            for visit in route.traversals:
                sid = visit.section_id
                if sid not in lengths:
                    path = getattr(context, 'products', {}).get(Key('path', sid))
                    if path is None:
                        from .paths import Path
                        path = Path(context.network.section_index[sid])
                    lengths[sid] = path.length
                total += lengths[sid]
            if total > maximum+1e-8:continue
            if neighbor in members:
                shortest = len(chain) if shortest is None else shortest
                candidates.append(chain)
            elif len(chain) < 4:
                pending.append((neighbor, (*chain, neighbor), total))
    closures = {tuple(sorted(chain)) for chain in candidates if len(chain) == shortest}
    if len(closures) > 1:
        raise ValueError('Paired partial terminal has ambiguous original neighboring source routes')
    return next(iter(closures), ())


def missing_terminal(context,bundle,consumed,members,available):
    """Preserve corresponding HEAD/TAIL ownership and real half-profile cuts.

    Only an existing opposite terminal, connected directly to this compound
    through original source incidence, may expand it. Scope cuts, incomplete
    members, disconnected intervals and implicit coordinate joins do not pass.
    """
    if len(bundle.members)!=2 or not all(member.complete for member in bundle.members):return ()
    sites={j.key:j for j in context.network.intersections if len(j.ends)>2}
    incident={end:j.key for j in sites.values() for end in j.ends}
    topology=RouteTopology(context.network);result=set();terminals=[]
    for member in bundle.members:
        for previous,following in zip(member.path,member.path[1:]):
            exit=End(previous.section_id,int(previous.forward))
            entry=End(following.section_id,int(not following.forward))
            owner=topology.incidents.get(exit)
            if owner is None or owner!=topology.incidents.get(entry):
                raise ValueError('Partial Bundle terminal ownership needs original consecutive source incidences')
        first,last=member.path[0],member.path[-1]
        terminals.append((End(first.section_id,0 if first.forward else 1),
                          End(last.section_id,1 if last.forward else 0)))
    for role,member in enumerate(bundle.members):
        indices=[i for i,visit in enumerate(member.path) if visit.section_id in consumed]
        if not indices or len(indices)==len(member.path) or indices!=list(range(indices[0],indices[-1]+1)):continue
        side=0 if indices[0]==0 else 1 if indices[-1]==len(member.path)-1 else None
        if side is None or incident.get(terminals[role][side]) not in members:continue
        end=terminals[1-role][side];owner=incident.get(end)
        if owner in members:continue
        if owner is None or owner not in available:
            raise ValueError('Paired partial terminal needs its available original physical owner')
        source=context.network.section_index[end.section_id]
        if end.end in source.scope_ends or end.end in source.selection_ends:
            raise ValueError('Paired partial terminal cannot cross a source-scope or selection boundary')
        linked=False
        for entry in sites[owner].ends:
            route=topology.trace(entry)
            if (route.terminal_site in members and not route.terminal_scope
                    and not route.terminal_selection):linked=True;break
        if linked:
            result.add(owner)
        else:
            closure = _local_source_closure(context, topology, sites, owner, members, available,
                max(context.network.section_index[v.section_id].width for m in bundle.members for v in m.path))
            if not closure:
                raise ValueError('Paired partial terminal lacks an original neighboring source route')
            result.update(closure)
    return tuple(sorted(result))
