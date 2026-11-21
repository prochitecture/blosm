"""Recover a paired source chain from a fragmented legacy Bundle inventory.

Legacy grouping can list every overlap fragment at both ends of a Bundle. Two
opposed one-way OSM paths that meet at their original end junctions then appear
to change lateral order. Original physical incidences and source ownership can
identify the complete two chains without joining controls by their positions.
This only maps topology; offset and smoothed-path validity remain surface checks.
"""
import math


def shared_terminal_order(bundle,streets,origins,tolerance):
    """Order two complete members when one original common end is a sort tie."""
    head,tail=tuple(bundle['head']),tuple(bundle['tail'])
    if origins is None or len(head)!=2 or set(head)!=set(tail):return None
    ends={}
    for sid in head:
        source=streets.get(sid)
        items=[i for i in source['items'] if i['kind']=='Section'] if source is not None else []
        if not items:return None
        ends[sid]=((items[0]['id'],tuple(items[0]['points'][0][:2])),
                   (items[-1]['id'],tuple(items[-1]['points'][-1][:2])))
    def site(sid,point):
        found=[]
        for original in origins.get('physical',()):
            if math.dist(point[:2],original['location'][:2])>tolerance:continue
            for connector in original['ring']:
                if connector.get('item_kind')!='Street' or len(connector['endpoints'])!=1:continue
                endpoint=connector['endpoints'][0]
                if endpoint['street_id']!=sid or endpoint['end'] not in(0,1):continue
                terminal,location=ends[sid][endpoint['end']]
                if (endpoint['terminal_kind']=='Section' and endpoint['terminal_id']==terminal
                        and math.dist(point[:2],location)<=tolerance):found.append(original['junction_id'])
        return found[0] if len(found)==1 else None
    points={name:{sid:tuple(point[:2]) for sid,point in zip(ids,bundle[name+'_locations'])}
            for name,ids in (('head',head),('tail',tail))}
    common=[]
    for name in('head','tail'):
        ids=[site(sid,points[name][sid]) for sid in head]
        common.append(ids[0] is not None and ids[0]==ids[1])
    if common.count(True)!=1:return None
    # Both extents must retain full original endpoint incidences; an unrelated
    # coincident control or a partial overlap cannot create this terminal tie.
    if any(not any(math.dist(points[name][sid],p)<=tolerance for _,p in ends[sid])
           for name in('head','tail') for sid in head):return None
    first,last=(tuple(sum(points[name][sid][k] for sid in head)*.5 for k in(0,1))
                for name in('head','tail'))
    direction=(last[0]-first[0],last[1]-first[1])
    separated=points['tail' if common[0] else 'head']
    delta=tuple(separated[head[0]][k]-separated[head[1]][k] for k in(0,1))
    side=direction[0]*delta[1]-direction[1]*delta[0]
    if abs(side)<=tolerance*max(1.,math.hypot(*direction)):return None
    order=head if side>0. else head[::-1]
    return dict(bundle,head=order,tail=order,
        head_locations=tuple(points['head'][sid] for sid in order),
        tail_locations=tuple(points['tail'][sid] for sid in order))


def paired_source_chains(bundle, streets, origins, tolerance):
    """Return normalized terminals only for two unambiguous original paths."""
    owned = tuple(bundle.get('owned_streets', ()))
    if (origins is None or len(owned) < 3 or len(set(owned)) != len(owned)
            or set(bundle['head']) != set(owned) or set(bundle['tail']) != set(owned)):
        return None
    groups, sections, names, categories = {}, {}, set(), set()
    for sid in owned:
        street = streets.get(sid)
        if street is None or street.get('bundle_id') != bundle['id']:
            return None
        items = street['items']
        if not items or any(item['kind'] != 'Section' for item in items):
            return None
        signatures = set()
        for item in items:
            source = tuple(item.get('source_way_ids', ()))
            tags = item.get('tags', {})
            segments = item.get('source_segments', ())
            if (len(source) != 1 or str(tags.get('oneway', '')).lower() not in ('yes', '1', 'true')
                    or not tags.get('name') or not segments):
                return None
            if (any(len(segment) < 5 or segment[0] != source[0]
                    or segment[3] != 0. or segment[4] != 1. for segment in segments)
                    or any(a[2] != b[1] for a, b in zip(segments, segments[1:]))):
                return None
            signatures.add(source[0]); names.add(tags['name']); categories.add(item.get('category'))
        if len(signatures) != 1:
            return None
        groups.setdefault(signatures.pop(), []).append(sid)
        sections[sid] = items
    if len(groups) != 2 or len(names) != 1 or len(categories) != 1:
        return None

    incidence, locations = {}, {}
    for site in origins.get('physical', ()):
        jid = site['junction_id']
        if jid in locations:
            return None
        locations[jid] = tuple(site['location'][:2])
        for connector in site['ring']:
            relevant = [e for e in connector['endpoints'] if e['street_id'] in sections]
            if relevant and (connector.get('item_kind') != 'Street' or len(connector['endpoints']) != 1):
                return None
            for endpoint in relevant:
                sid, end = endpoint['street_id'], endpoint['end']
                if end not in (0, 1) or endpoint['terminal_kind'] != 'Section' or (sid, end) in incidence:
                    return None
                section = sections[sid][-1 if end else 0]
                point = section['points'][-1 if end else 0]
                if (endpoint['terminal_id'] != section['id']
                        or math.dist(point[:2], locations[jid]) > tolerance):
                    return None
                incidence[sid, end] = jid
    if set(incidence) != {(sid, end) for sid in owned for end in (0, 1)}:
        return None

    chains = []
    for source, ids in sorted(groups.items()):
        incoming, outgoing = {}, {}
        for sid in ids:
            start, stop = incidence[sid, 0], incidence[sid, 1]
            if start == stop or start in outgoing or stop in incoming:
                return None
            outgoing[start], incoming[stop] = sid, sid
        starts, stops = set(outgoing)-set(incoming), set(incoming)-set(outgoing)
        if len(starts) != 1 or len(stops) != 1:
            return None
        start, stop = starts.pop(), stops.pop()
        chain, sites, seen, current = [], [start], set(), start
        while current in outgoing:
            sid = outgoing[current]
            if sid in seen:
                return None
            seen.add(sid); chain.append(sid)
            current = incidence[sid, 1]; sites.append(current)
        if current != stop or seen != set(ids):
            return None
        # Source node identities must agree with each directed physical join.
        for before, after in zip(chain, chain[1:]):
            a = sections[before][-1]['source_segments'][-1]
            b = sections[after][0]['source_segments'][0]
            if a[0] != source or b[0] != source or a[2] != b[1]:
                return None
        chains.append((chain, sites))
    a, b = chains
    if (a[1][0] != b[1][-1] or a[1][-1] != b[1][0]
            or set(a[1][1:-1]) & set(b[1]) or set(b[1][1:-1]) & set(a[1])):
        return None
    head, tail = sorted((a[1][0], a[1][-1]))
    oriented = []
    for chain, sites in chains:
        forward = sites[0] == head
        ordered = chain if forward else list(reversed(chain))
        points = []
        for sid in ordered:
            local = [tuple(p[:2]) for item in sections[sid] for p in item['points']]
            if not forward:
                local.reverse()
            points.extend(local[1:] if points else local)
        oriented.append((ordered, forward, points))
    # A simple control polygon establishes lateral order. Its endpoints alone
    # may touch, and only because their original physical owners agree above.
    first, second = oriented[0][2], oriented[1][2]
    rim = first + list(reversed(second))[1:-1]
    edges = list(zip(rim, rim[1:]+rim[:1]))
    if len(edges)*(len(edges)-1)//2>200000:
        return None
    def cross(u, v):
        return u[0]*v[1]-u[1]*v[0]
    def cuts(a, b, c, d):
        u, v, w = (b[0]-a[0], b[1]-a[1]), (d[0]-c[0], d[1]-c[1]), (c[0]-a[0], c[1]-a[1])
        lu, lv = math.hypot(*u), math.hypot(*v)
        if min(lu, lv) <= tolerance:
            return True
        determinant = cross(u, v)
        if abs(determinant) <= tolerance*max(lu, lv):
            if abs(cross(w, u)) > tolerance*lu:
                return False
            t = [((p[0]-a[0])*u[0]+(p[1]-a[1])*u[1])/(lu*lu) for p in (c, d)]
            return max(min(t), 0.) <= min(max(t), 1.)+tolerance/lu
        t, s = cross(w, v)/determinant, cross(w, u)/determinant
        return -tolerance/lu <= t <= 1.+tolerance/lu and -tolerance/lv <= s <= 1.+tolerance/lv
    for i, (a, b) in enumerate(edges):
        for j in range(i+1, len(edges)):
            if j == i+1 or (i == 0 and j == len(edges)-1):
                continue
            if cuts(a, b, *edges[j]):
                return None
    origin = rim[0]
    area = sum(cross((a[0]-origin[0], a[1]-origin[1]), (b[0]-origin[0], b[1]-origin[1])) for a, b in edges)
    if abs(area) <= tolerance*tolerance:
        return None
    if area > 0.:
        oriented.reverse()
    return dict(bundle, head=tuple(chain[0] for chain, _, _ in oriented),
        tail=tuple(chain[-1] for chain, _, _ in oriented),
        head_locations=(locations[head], locations[head]), tail_locations=(locations[tail], locations[tail]))
