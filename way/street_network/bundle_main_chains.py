"""Explicit original carriageway continuations inside broad Bundle candidates.

Parallel preparation also owns frontage roads and attaching links. Those are
not automatically additional lateral members. A unique highest-class pair at
the declared ends can retain its source-proven longitudinal continuations;
all other original Streets still need independent interior/exterior review.
"""
import math


_CLASSES=('motorway','trunk','primary','secondary','tertiary','unclassified',
          'residential','living_street','service')


def normalize_main_chains(bundle,streets,origins,tolerance):
    if origins is None:return None
    declared=set(bundle.get('owned_streets',()))|set(bundle['head'])|set(bundle['tail'])
    if not declared or any(sid not in streets for sid in declared):return None
    items={sid:[i for i in streets[sid]['items'] if i['kind']=='Section'] for sid in declared}
    if any(not parts for parts in items.values()):return None
    def category(sid):
        lengths={}
        for item in items[sid]:
            kind=item.get('category');length=sum(math.dist(a[:2],b[:2]) for a,b in zip(item['points'],item['points'][1:]))
            lengths[kind]=lengths.get(kind,0.)+length
        return max(lengths,key=lengths.get)
    kinds={sid:category(sid) for sid in declared}
    def rank(kind):
        base=kind[:-5] if kind and kind.endswith('_link') else kind
        return (_CLASSES.index(base),bool(kind.endswith('_link'))) if base in _CLASSES else (100,False)
    best=min((kinds[sid] for sid in bundle['head']),key=rank)
    if rank(best)[0]==100:return None
    heads=[(sid,point) for sid,point in zip(bundle['head'],bundle['head_locations']) if kinds[sid]==best]
    tails=[(sid,point) for sid,point in zip(bundle['tail'],bundle['tail_locations']) if kinds[sid]==best]
    derive_tail=not tails and len(heads)==2
    if len(heads)!=2 or (len(tails)!=2 and not derive_tail) or heads[0][0]==heads[1][0]:return None
    if tails and tails[0][0]==tails[1][0]:return None
    if derive_tail:
        names=[{item.get('tags',{}).get('name') for item in items[sid]}-{None,''} for sid,_ in heads]
        if not names[0] or names[0]!=names[1]:return None
    # A frontage or link can only bridge complete main-class source Streets.
    # It cannot supply an invented lateral terminal of a different category.
    core={sid for sid,kind in kinds.items() if kind==best}
    links={sid for sid,kind in kinds.items() if kind==best+'_link'} if not best.endswith('_link') else set()
    eligible=core|links
    if any(streets[sid].get('bundle_id')!=bundle['id'] for sid in eligible):return None
    positions={sid:(tuple(parts[0]['points'][0][:2]),tuple(parts[-1]['points'][-1][:2])) for sid,parts in items.items()}
    terminals={sid:(parts[0]['id'],parts[-1]['id']) for sid,parts in items.items()}
    def endpoint(sid,point):
        ends=[end for end in(0,1) if math.dist(point[:2],positions[sid][end])<=tolerance]
        return (sid,ends[0]) if len(ends)==1 else None
    starts=[endpoint(*record) for record in heads];stops=[endpoint(*record) for record in tails]
    if None in starts or None in stops or len(set(starts+stops))!=len(starts+stops):return None
    sites={};incidence={}
    for site in origins.get('physical',()):
        ring=[]
        for connector in site['ring']:
            for e in connector['endpoints']:
                sid,end=e['street_id'],e['end']
                if sid not in eligible:continue
                if (connector.get('item_kind')!='Street' or len(connector['endpoints'])!=1 or end not in(0,1)
                        or e['terminal_kind']!='Section' or e['terminal_id']!=terminals[sid][end]
                        or math.dist(positions[sid][end],site['location'][:2])>tolerance
                        or (sid,end) in incidence):return None
                incidence[sid,end]=site['junction_id'];ring.append((sid,end))
        sites[site['junction_id']]=tuple(ring)
    def tangent(port,arriving):
        sid,end=port;points=items[sid][0 if end==0 else -1]['points']
        a,b=(points[0],points[1]) if end==0 else (points[-1],points[-2])
        v=(b[0]-a[0],b[1]-a[1]);length=math.hypot(*v)
        if length<=tolerance:return None
        return tuple(((-1. if arriving else 1.)/length)*value for value in v)
    def source(port):
        sid,end=port
        return set(items[sid][0 if end==0 else -1].get('source_way_ids',()))
    def candidates(exit,options):
        options=[p for p in options if p[1]!=exit[1]]
        same=[p for p in options if source(exit)&source(p)]
        if len(same)==1:return same
        incoming=tangent(exit,True)
        if incoming is None:return []
        aligned=[]
        for p in options:
            outgoing=tangent(p,False)
            if outgoing is not None and sum(a*b for a,b in zip(incoming,outgoing))>=math.cos(math.radians(20.)):
                aligned.append(p)
        main=[p for p in aligned if p[0] in core]
        return main if main else aligned
    def continuation(exit,options,forbidden,remaining):
        options=candidates(exit,options)
        if len(options)==1:return options[0]
        def reaches(port):
            pending=[port];seen=set()
            for current in pending:
                if len(seen)>512:return False
                sid,end=current
                if sid in forbidden or current in seen:continue
                seen.add(current);stop=(sid,1-end)
                if stop in remaining:return True
                ring=sites.get(incidence.get(stop),())
                pending.extend(p for p in candidates(stop,[p for p in ring if p!=stop
                    and p[0] not in forbidden and p not in starts and p not in stops]) if p not in seen)
            return False
        possible=[p for p in options if reaches(p)]
        return possible[0] if len(possible)==1 else None
    used=set();chains=[];ordered_stops=[]
    for start in starts:
        current=start;chain=[]
        for _ in range(128):
            sid,end=current
            if sid in used:return None
            used.add(sid);chain.append(sid)
            exit=(sid,1-end)
            if exit in stops:
                ordered_stops.append(exit);break
            ring=sites.get(incidence.get(exit),())
            options=[p for p in ring if p!=exit and p[0] not in used]
            if derive_tail and not candidates(exit,options):
                if sid not in core:return None
                ordered_stops.append(exit);break
            following=continuation(exit,options,used,set(stops)-set(ordered_stops))
            if following is None or following in starts or following in stops:return None
            current=following
        else:return None
        chains.append(tuple(chain))
    if len(set(ordered_stops))!=2:return None
    if derive_tail:
        if not core<=used:return None
        end_points=[positions[sid][end] for sid,end in ordered_stops]
        width=max(item.get('width',0.) for sid in core for item in items[sid])
        if not width>0. or math.dist(*end_points)>4.*width:return None
        tail_points={port:positions[port[0]][port[1]] for port in ordered_stops}
    else:tail_points={endpoint(sid,point):tuple(point[:2]) for sid,point in tails}
    normalized=dict(bundle,head=tuple(p[0] for p in starts),tail=tuple(p[0] for p in ordered_stops),
        head_locations=tuple(tuple(point[:2]) for _,point in heads),
        tail_locations=tuple(tail_points[p] for p in ordered_stops),owned_streets=tuple(sorted(used)),
        source_extent_derived=derive_tail)
    return normalized,tuple(sorted(declared-used))
