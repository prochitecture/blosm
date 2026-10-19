"""Restore complete near-boundary OSM junctions in the prepared input layer.

Original source incidence establishes connectivity. Existing Street/Section IDs,
controls and dimensions survive; additions get ordinary network IDs and explicit
OSM provenance. No fillets, mesh samples or surface geometry are calculated here.
The shared native constructor remains responsible for geometric acceptance.
"""
from copy import deepcopy
from dataclasses import dataclass
import math

from .physical_ports import GroundSurfacePolicy


class BoundaryContextError(ValueError):
    pass


@dataclass(frozen=True)
class BoundarySupportPolicy:
    """Conservative straight-ray support, followed by native construction guards.

    Adjacent offset curb rays give a tangent-distance estimate. A dimensional
    margin leaves room beyond it for the full approach profile. Reflex and
    exactly opposing pairs use the dimensional baseline; unequal opposing
    streets still require the native transition constructor. This is a bounded
    source-support heuristic, not proof of a feasible corner or terrain fit.
    """
    clearance: float = 1.4

    def lengths(self, arms):
        if not math.isfinite(self.clearance) or self.clearance <= 0:
            raise BoundaryContextError('Boundary support clearance must be positive and finite')
        result = [a['width']+2*a['radius']+self.clearance for a in arms]
        order = sorted(range(len(arms)), key=lambda i:arms[i]['angle'])
        for ia, ib in zip(order, order[1:]+order[:1]):
            a,b = arms[ia],arms[ib]
            gap = (b['angle']-a['angle']) % math.tau
            if min(gap, math.tau-gap) < 1e-6:
                raise BoundaryContextError('Outside junction has coincident source directions')
            if gap >= math.pi-1e-6:
                continue
            ha = a['width']*.5+abs(a['offset'])
            hb = b['width']*.5+abs(b['offset'])
            radius = max(a['radius'],b['radius'])
            sa = (hb+radius+math.cos(gap)*(ha+radius))/math.sin(gap)
            sb = (ha+radius+math.cos(gap)*(hb+radius))/math.sin(gap)
            result[ia] = max(result[ia],sa+self.clearance)
            result[ib] = max(result[ib],sb+self.clearance)
        if any(not math.isfinite(v) or v <= 0 for v in result):
            raise BoundaryContextError('Outside junction needs unsupported source support')
        return result


class BoundaryJunctionAdapter:
    """Augment a snapshot atomically per authentic outside OSM node.

    The first implementation needs each missing leg's dimensions to be proven
    by a retained Section of the same OSM way, category and tags. It never
    guesses a new way's width, extends a Bundle independently, or treats a
    grade-separated crossing as a ground junction. Insufficient adjacent source
    support is diagnosed rather than extrapolated or truncated to fit.
    """
    policy = 'COMPLETE_ORIGINAL_BOUNDARY_JUNCTIONS_V1'

    def __init__(self, design, category_properties, *, support=None, ground_passages=False):
        self.design = design
        self.categories = set(design.categories)
        self.properties = category_properties
        self.support = support or BoundarySupportPolicy()
        self.ground_passages = ground_passages

    @staticmethod
    def _all_streets(snapshot):
        return {**snapshot.streets, **snapshot.retained_physical_streets}

    @staticmethod
    def _sections(street):
        return [s for s in street['items'] if s['kind']=='Section']

    @classmethod
    def _next_ids(cls, snapshot):
        streets = cls._all_streets(snapshot)
        sections = [s['id'] for street in streets.values() for s in cls._sections(street)]
        origins = snapshot.junction_origins
        sites = set(snapshot.intersections) | set(snapshot.junctions)
        sites.update(i['id'] for s in streets.values() for i in s['items'] if i['kind']=='Intersection')
        sites.update(r['junction_id'] for r in origins.get('physical', ()))
        sites.update(r['junction_id'] for r in origins.get('seams', ()))
        sites.update(r['junction']['junction_id'] for r in origins.get('logical', ()))
        connectors = [r['connector_id'] for site in origins.get('physical', ()) for r in site['ring']]
        connectors.extend(c.source_id for j in snapshot.junctions.values() for c in j.connectors)
        connectors.extend(a['connector_id'] for s in streets.values() for i in s['items']
                          for a in i.get('attachments', ()))
        connectors.extend(c['connector_id'] for r in origins.get('logical', ()) for c in r['junction']['ring'])
        return dict(street=max(streets,default=-1)+1, section=max(sections,default=-1)+1,
                    junction=max(sites,default=-1)+1, connector=max(connectors,default=-1)+1)

    def _dimensions(self, snapshot, leg):
        matches = []
        for street in self._all_streets(snapshot).values():
            for section in self._sections(street):
                if (section['category'] != leg['category'] or section.get('tags',{}) != leg['tags']
                        or leg['source_way_id'] not in map(str,section.get('source_way_ids',()))):
                    continue
                signs = {1 if r[4]>r[3] else -1 for r in section.get('source_segments',())
                         if str(r[0]) == leg['source_way_id']}
                if len(signs) != 1:
                    continue
                sign, = signs
                matches.append((float(section['width']),float(section.get('offset') or 0.)*sign))
        if not matches or len(set(matches)) != 1:
            raise BoundaryContextError('Outside leg lacks unambiguous retained dimensions for its original OSM way')
        width,offset = matches[0]
        if not math.isfinite(width) or width<=0 or not math.isfinite(offset):
            raise BoundaryContextError('Outside leg has invalid physical dimensions')
        # A newly added source runs outward from the junction. Its offset must
        # follow that direction, regardless of the original OSM way direction.
        return width,offset*(1 if leg['source_end']==0 else -1)

    def _plan(self, snapshot, node, ids):
        if node['topology_kind'] != 'JUNCTION' or not node.get('open_ends'):
            return None
        ends = node['open_ends']
        if any(e['owner_scope'] != 'ACTIVE' for e in ends):
            raise BoundaryContextError('Outside junction includes retained physical ownership needing a shared recovery plan')
        if len({e['street_id'] for e in ends}) != len(ends):
            raise BoundaryContextError('Outside recovery would close a Street and needs a cyclic path provider')
        keys = [(r['source_way_id'],r['other_node_id']) for r in node['legs']]
        if len(keys)!=len(set(keys)) or len(keys)<3:
            raise BoundaryContextError('Outside junction needs at least three unique original incidences')
        position = tuple(node['position'])+(0.,)
        cut_by_leg = {}
        for end in ends:
            first,last = end['source_nodes']
            if str(node['osm_node_id']) not in (first,last):
                raise BoundaryContextError('Outside node does not belong to its source incidence')
            key = end['source_way_id'],last if first==str(node['osm_node_id']) else first
            if key in cut_by_leg:
                raise BoundaryContextError('An outside incidence has multiple prepared endpoint owners')
            cut_by_leg[key] = end
        if set(cut_by_leg)-set(keys):
            raise BoundaryContextError('Prepared cut has no original outside leg')
        arms = []
        for key,leg in zip(keys,node['legs']):
            if (leg['category'] not in self.categories or GroundSurfacePolicy.reason(
                    leg['tags'],building_passages=self.ground_passages)
                    or str(leg['tags'].get('level','0')) not in ('0','0.0')):
                raise BoundaryContextError('Outside junction includes a category or elevation needing another feature provider')
            if leg.get('source_end') not in (0,1) or len(leg.get('source_nodes',()))!=2:
                raise BoundaryContextError('Outside source support lacks original direction provenance')
            if (str(leg['source_nodes'][leg['source_end']]) != str(node['osm_node_id'])
                    or str(leg['source_nodes'][1-leg['source_end']]) != str(leg['other_node_id'])):
                raise BoundaryContextError('Outside source direction differs from its original node incidence')
            delta = tuple(leg['other_position'][i]-position[i] for i in range(2))
            length = math.hypot(*delta)
            if length <= 1e-6:
                raise BoundaryContextError('Outside source leg has zero length')
            end = cut_by_leg.get(key)
            if end is not None:
                street = snapshot.streets[end['street_id']]
                if street['status'] != 'SUPPORTED' or street.get('bundle_id') is not None:
                    raise BoundaryContextError('Outside junction needs supported unbundled source Streets')
                section = self._sections(street)[0 if end['end']==0 else -1]
                if section['id']!=end['section_id'] or street['endpoints'][end['end']] is not None:
                    raise BoundaryContextError('Outside endpoint no longer has its recorded terminal ownership')
                if section.get('clip_start' if end['end']==0 else 'clip_end',0.):
                    raise BoundaryContextError('Outside recovery conflicts with an explicit selection trim')
                record = section['source_segments'][0 if end['end']==0 else -1]
                point = section['points'][0 if end['end']==0 else -1]
                cuts = [c for c in node['incidences'] if
                        (c['source_way_id'],*c['source_nodes'],c['parameter']) ==
                        (end['source_way_id'],*end['source_nodes'],end['parameter'])]
                if (tuple(map(str,record[:3])) != (end['source_way_id'],*end['source_nodes'])
                        or tuple(end['source_nodes']) != tuple(leg['source_nodes'])
                        or record[3 if end['end']==0 else 4] != end['parameter']
                        or len(cuts)!=1 or math.dist(point[:2],cuts[0]['clipped_point'])>.001):
                    raise BoundaryContextError('Outside incidence differs from the retained source interval')
                width,offset = float(section['width']),float(section.get('offset') or 0.)
            else:
                width,offset = self._dimensions(snapshot,leg)
                if any(n['osm_node_id']==leg['other_node_id'] and n['topology_kind']=='JUNCTION'
                       for n in snapshot.clipped_source_context['nodes']):
                    raise BoundaryContextError('Neighboring outside junctions need coordinated source support')
            radius = float(self.properties[leg['category']]['radius'])
            if not all(math.isfinite(v) for v in (width,offset,radius)) or min(width,radius)<=0:
                raise BoundaryContextError('Outside junction has invalid width or corner radius')
            arms.append(dict(leg=leg,end=end,width=width,offset=offset,radius=radius,
                             angle=math.atan2(delta[1],delta[0])%math.tau,length=length))
        support = self.support.lengths(arms)
        # Only recover a node within its dimensional influence region. Distant
        # outside intersections do not cause an unbounded import expansion.
        influence = max(s+a['width']*.5+abs(a['offset']) for a,s in zip(arms,support))
        if node['boundary_distance_m'] > influence:
            return dict(status='OUTSIDE_INFLUENCE',osm_node_id=node['osm_node_id'],influence_m=influence)
        used = dict(ids)
        jid = used['junction']; used['junction'] += 1
        changed = {}; added = {}; ring = []; source_extensions = []
        for arm,distance in zip(arms,support):
            leg,end = arm['leg'],arm['end']
            if end is not None:
                sid = end['street_id']
                street = changed.setdefault(sid,deepcopy(snapshot.streets[sid]))
                section = self._sections(street)[0 if end['end']==0 else -1]
                points = list(section['points'])
                # Keep the old clipped control as an ordinary authored control.
                points.insert(0,position) if end['end']==0 else points.append(position)
                section['points'] = tuple(points)
                target = float(leg['source_end'])
                interval = (target,end['parameter']) if end['end']==0 else (end['parameter'],target)
                extension = (end['source_way_id'],*end['source_nodes'],*interval)
                records = list(section['source_segments'])
                records.insert(0,extension) if end['end']==0 else records.append(extension)
                section['source_segments'] = tuple(records)
                endpoints = list(street['endpoints']); endpoints[end['end']] = jid; street['endpoints'] = tuple(endpoints)
                section_id = section['id']; terminal_end = end['end']
            else:
                # Do not invent support beyond the adjacent retained segment.
                if distance >= arm['length']-1e-4:
                    raise BoundaryContextError('Outside leg needs more source context than its retained adjacent segment')
                sid = used['street']; used['street'] += 1
                section_id = used['section']; used['section'] += 1
                fraction = distance/arm['length']
                far = tuple(leg['other_position'])+(0.,)
                start = float(leg['source_end']); stop = 1.-start
                output_stop = start+fraction*(1 if start==0 else -1)
                trim = arm['length']-distance
                section = dict(kind='Section',id=section_id,points=(position,far),category=leg['category'],
                    width=arm['width'],offset=arm['offset'],tags=deepcopy(leg['tags']),
                    source_way_ids=[leg['source_way_id']],
                    source_segments=((leg['source_way_id'],*leg['source_nodes'],start,stop),),
                    clip_end=trim,
                    width_origin='Retained same-way dimensions for bounded outside support')
                added[sid] = dict(id=sid,items=[section],status='SUPPORTED',reason=None,
                    endpoints=(jid,None),bundle_id=None,source_role='BOUNDARY_SUPPORT')
                terminal_end = 0
                source_extensions.append(dict(street_id=sid,section_id=section_id,
                    source_way_id=leg['source_way_id'],source_nodes=leg['source_nodes'],
                    source_interval=(start,stop),output_source_interval=(start,output_stop),
                    support_length_m=distance,trim_end_m=trim))
            cid = used['connector']; used['connector'] += 1
            ring.append((arm['angle'],dict(connector_id=cid,item_id=sid,item_kind='Street',
                leaving=terminal_end==0,endpoints=[dict(street_id=sid,end=terminal_end,
                terminal_id=section_id,terminal_kind='Section')])))
        if max(used.values()) >= 2**31:
            raise BoundaryContextError('Outside recovery exhausts network identity capacity')
        ring.sort(key=lambda row:row[0])
        return dict(status='RESTORED',osm_node_id=node['osm_node_id'],junction_id=jid,
                    position=position,changed=changed,added=added,ids=used,
                    ring=[r for _,r in ring],ends=ends,extensions=source_extensions,
                    influence_m=influence)

    def _junction(self, snapshot, plan):
        from ..gn.kernel import inputs
        sources,connectors,profiles = [],[],[]
        all_streets = {**snapshot.streets,**plan['changed'],**plan['added']}
        for rank,row in enumerate(plan['ring']):
            endpoint, = row['endpoints']
            street = all_streets[endpoint['street_id']]
            section = next(s for s in self._sections(street) if s['id']==endpoint['terminal_id'])
            end = endpoint['end']
            points = section['points'][:2] if end==0 else section['points'][-2:]
            # Full source paths remain in the snapshot. Only the local frame
            # uses these two exact neighboring controls, as in ordinary intake.
            sources.append(inputs.StreetInput(street['id'],(inputs.SectionInput(
                section['id'],points,section['width'],section.get('offset') or 0.,section['category']),)))
            profiles.append(inputs.ProfileInput(self.design.profile_rise,
                inputs.SMOOTH_CROWN if self.design.profile_rise else inputs.FLAT))
            connectors.append(inputs.ConnectorInput(row['connector_id'],
                (inputs.MemberInput(rank,end==0),),rank,
                radius=self.properties[section['category']]['radius']))
        result = inputs.PreparedJunction(plan['junction_id'],plan['position'],tuple(sources),
                                        tuple(connectors),tuple(profiles),0.,self.design.cell_size)
        inputs.validate(result)
        return result

    def apply(self, snapshot):
        context = snapshot.clipped_source_context
        if not context:
            return
        if 'boundary_junction_recovery' in snapshot.preparation_audit:
            raise BoundaryContextError('Boundary junction recovery already ran on this snapshot')
        report = dict(policy=self.policy,output_extent='COMPLETE_JUNCTION_WITH_BOUNDED_OUTSIDE_SUPPORT',
                      nodes=[],geometry_qualified=False)
        snapshot.preparation_audit['boundary_junction_recovery'] = report
        if context.get('version')!=1 or context.get('issues'):
            report['status'] = 'UNRESOLVED_CONTEXT'
            snapshot.issue('BoundaryJunction',-1,'BOUNDARY_CONTEXT_UNRESOLVED',
                           'Outside source context needs unambiguous original incidence before recovery')
            return
        ids = self._next_ids(snapshot)
        for node in sorted(context['nodes'],key=lambda n:str(n['osm_node_id'])):
            try:
                plan = self._plan(snapshot,node,ids)
                if plan is None:
                    continue
                if plan['status'] != 'RESTORED':
                    report['nodes'].append(plan); continue
                # All validation precedes this commit. A failed node cannot
                # leave half an intersection or consume persistent identities.
                jid = plan['junction_id']
                junction = self._junction(snapshot,plan)
                # Physical incidence is authoritative even when the older
                # local kernel's stricter angle domain cannot represent it.
                # This matches ordinary Blosm preparation: the original-span
                # planner can use the physical ring and prove its own guards.
                from ..gn.kernel import inputs
                local_reason = None
                try:inputs.validate_initial_scope(junction)
                except ValueError as exc:local_reason=str(exc)
                origins = deepcopy(snapshot.junction_origins)
                origins['physical'] = list(origins.get('physical',()))+[
                    dict(junction_id=jid,location=plan['position'],ring=plan['ring'],
                         source_osm_node_id=node['osm_node_id'],source_role='RECOVERED_ORIGINAL_OSM')]
                replaced = {(e['street_id'],e['section_id'],e['end']) for e in plan['ends']}
                origins['free_endpoints'] = [r for r in origins.get('free_endpoints',()) if
                    (r['endpoint']['street_id'],r['endpoint']['terminal_id'],r['endpoint']['end']) not in replaced]
                origins['selection_boundaries'] = list(origins.get('selection_boundaries',()))
                for street in plan['added'].values():
                    section, = street['items']
                    # The actual OSM segment continues. Keep its full adjacent
                    # control and declare an explicit metric context cut; this
                    # endpoint is not an invented dead end or free-end alias.
                    origins['selection_boundaries'].append(dict(endpoint=dict(
                        street_id=street['id'],section_id=section['id'],end=1),trim_m=section['clip_end']))
                snapshot.streets.update(plan['changed']); snapshot.streets.update(plan['added'])
                if local_reason is None:snapshot.junctions[jid] = junction
                snapshot.intersections[jid] = dict(id=jid,location=plan['position'],
                    order=len(plan['ring']),status='SUPPORTED' if local_reason is None else 'UNSUPPORTED',
                    streets=tuple(r['item_id'] for r in plan['ring']),reason=local_reason,
                    source_osm_node_id=node['osm_node_id'],source_role='BOUNDARY_RECOVERY')
                snapshot.junction_origins = origins
                ids = plan['ids']
                record = {k:plan[k] for k in ('status','osm_node_id','junction_id','influence_m','extensions')}
                record['local_kernel_scope_reason'] = local_reason
                report['nodes'].append(record)
            except (BoundaryContextError,KeyError,TypeError,ValueError) as exc:
                report['nodes'].append(dict(status='UNSUPPORTED',osm_node_id=node['osm_node_id'],reason=str(exc)))
                snapshot.issue('BoundaryJunction',-1,'BOUNDARY_JUNCTION_UNSUPPORTED',
                               str(node['osm_node_id'])+': '+str(exc))
        report['status'] = 'PREPARED_WITH_UNSUPPORTED_CONTEXT' if any(
            n['status']=='UNSUPPORTED' for n in report['nodes']) else 'PREPARED'
        # Preserve previous attachment/remapping evidence while keeping the
        # authoritative Section inventory counts current after additions.
        if hasattr(snapshot,'ownership_audit'):
            sections = [s['id'] for street in snapshot.streets.values() for s in self._sections(street)]
            snapshot.ownership_audit.update(section_count=len(set(sections)),
                                            section_occurrence_count=len(sections))
