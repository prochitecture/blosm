"""Identify source cross-links inside aggregated Bundle-end intersections.

This is a semantic review, not a surface-coverage certificate. Both original
connector incidences and every source Section remain in the prepared snapshot.
A later provider must prove that their complete footprint is covered by the
junction surface before replacing their separate carriageway generation.
"""
from dataclasses import asdict, dataclass
import math


@dataclass(frozen=True)
class JunctionInternalLink:
    junction_id: int
    street_id: int
    connector_ids: tuple[int, int]  # Authored start/end, not circular ring order.
    section_ids: tuple[int, ...]
    bundle_id: int
    bundle_side: int
    endpoint_members: tuple[int, int]  # Member order at that Bundle end.

    def report(self):
        return dict(status='PENDING_FOOTPRINT', **asdict(self))


@dataclass(frozen=True)
class PreparedCompoundJunction:
    """Complete original sources and a separate external approach contract.

    Planning this input does not mark the source Street generated. A component
    must assign it to this patch exactly once, bind live source controls and
    retain the native full-footprint guard before claiming generated coverage.
    """
    junction: object
    internal_streets: tuple
    links: tuple[JunctionInternalLink, ...]
    original_connector_ids: tuple[int, ...]
    retained: tuple


class PreparedJunctionInteriorAdapter:
    """Reusable source adapter, independent of bpy and particular OSM inputs.

    Inject the contract/connector adapter during candidate qualification; the
    normal provider can use its packaged equivalent once that asset is enabled.
    section_reader is the same provenance-preserving reader used for external
    Streets. No radii, mouth stations, curves or surface vertices are calculated.
    """
    def __init__(self, contract, connector_adapter, section_reader):
        self.contract=contract
        self.connector_adapter=connector_adapter
        self.section_reader=section_reader

    def read(self, snapshot, intersection, *, profiles_by_connector, **options):
        raw=self.connector_adapter._connectors(intersection)
        links,retained=PreparedJunctionInteriorReview().read(snapshot,intersection,raw)
        if not links:return None
        internal={cid for link in links for cid in link.connector_ids}
        if len(internal)!=2*len(links):
            raise ValueError('An internal source incidence has more than one owner')
        external=tuple(c for c in raw if c.id not in internal)
        sources=[]
        for link in links:
            street=next(c.item for c in raw if c.id==link.connector_ids[0])
            sections=tuple(self.section_reader(street))
            if tuple(s.source_id for s in sections)!=link.section_ids:
                raise ValueError('Internal Section order or identity differs from the original source')
            original=snapshot.streets[link.street_id]['items']
            for section,record in zip(sections,original):
                points=tuple(tuple(p) for p in record['points'])
                if (tuple(tuple(p) for p in section.centerline)!=points or
                    section.width!=record['width'] or section.offset!=(record['offset'] or 0.)):
                    raise ValueError('Internal source geometry and dimensions must remain unchanged')
            sources.append(self.contract.StreetInput(link.street_id,sections))
        junction=self.connector_adapter.from_blosm_intersection(intersection,
            profiles_by_connector=profiles_by_connector,section_reader=self.section_reader,
            source_connectors=external,**options)
        return PreparedCompoundJunction(junction,tuple(sources),links,tuple(c.id for c in raw),retained)


class PreparedJunctionInteriorReview:
    """Recognize finite cross-links; never infer one from a duplicate ID alone."""
    tolerance = 1e-5

    def _straight_chain(self, sections):
        """Return endpoints only for a connected, finite, monotone source chain.

        Collinearity alone also accepts disconnected pieces and paths which
        double back. Neither may be silently reduced to one internal link.
        Duplicate samples at an authored Section seam are expected; duplicate
        edges inside a Section are not.
        """
        if not sections or len({s['id'] for s in sections}) != len(sections):
            return None
        points = []
        previous = None
        for section in sections:
            chain = section['points']
            if len(chain) < 2 or any(
                len(p) < 2 or not all(math.isfinite(v) for v in p) for p in chain
            ):
                return None
            if previous is not None and math.dist(previous[:2], chain[0][:2]) > self.tolerance:
                return None
            if any(math.dist(a[:2], b[:2]) <= self.tolerance for a, b in zip(chain, chain[1:])):
                return None
            points.extend(chain)
            previous = chain[-1]
        a, b = points[0][:2], points[-1][:2]
        delta = tuple(y-x for x, y in zip(a, b))
        length = math.hypot(*delta)
        if length <= self.tolerance:
            return None
        direction = tuple(v/length for v in delta)
        previous_station = -math.inf
        for point in points:
            station = sum((point[i]-a[i])*direction[i] for i in (0, 1))
            nearest = tuple(a[i]+max(0., min(length, station))*direction[i] for i in (0, 1))
            if (math.dist(point[:2], nearest) > self.tolerance
                or station < previous_station-self.tolerance):
                return None
            previous_station = station
        return a, b

    def read(self, snapshot, intersection, connectors):
        grouped={};bundle_ends=[]
        for connector in connectors:
            item=connector.item
            if hasattr(item,'streetsHead') and hasattr(item,'streetsTail'):
                points=item.headLocs if connector.leaving else item.tailLocs
                if len(points)==2:
                    bundle_ends.append((item.id,0 if connector.leaving else 1,
                                        tuple(tuple(float(v) for v in p[:2]) for p in points)))
            else:
                grouped.setdefault(item.id,[]).append(connector)
        links=[];retained=[]
        for sid,incidences in sorted(grouped.items()):
            if len(incidences)<2:continue
            reason=None;source=snapshot.streets.get(sid)
            if len(incidences)!=2 or {c.leaving for c in incidences}!={False,True}:
                reason='Repeated incidences are not the two distinct authored Street ends'
            elif source is None or source.get('endpoints')!=(intersection.id,intersection.id):
                reason='Both original Street ends must belong to this logical junction'
            elif source.get('status')!='SUPPORTED' or any(i['kind']!='Section' for i in source['items']):
                reason='The cross-link needs a separate feature or item treatment'
            if reason:
                retained.append(dict(street_id=sid,status='RETAINED_APPROACH',reason=reason));continue
            sections=source['items'];chain=self._straight_chain(sections)
            matches=[]
            if chain is not None:
                a,b=chain
                for bid,side,ends in bundle_ends:
                    for first,last in ((0,1),(1,0)):
                        if math.dist(a,ends[first])<=self.tolerance and math.dist(b,ends[last])<=self.tolerance:
                            matches.append((bid,side,(first,last)))
            if len(matches)!=1:
                retained.append(dict(street_id=sid,status='RETAINED_APPROACH',
                    reason='No unique finite cross-link between the two members of one Bundle end'))
                continue
            bid,side,members=matches[0]
            ends={int(not c.leaving):c.id for c in incidences}
            links.append(JunctionInternalLink(intersection.id,sid,(ends[0],ends[1]),
                tuple(s['id'] for s in sections),bid,side,members))
        return tuple(links),tuple(retained)
