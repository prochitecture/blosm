"""Original bank provenance and bounded physical junction reconciliation."""
from dataclasses import dataclass,replace
import math
import numpy as np
from .model import End,Key
from .paths import Path
from .routes import Route
from .owned_prefixes import OwnedPrefixRoute
from .intersections import Approach
from .compounds import rim_label
from .coincident_source_curbs import _intervals,source_station,curve_parameter
from .source_curb_coordinates import CurbCoordinates
from .source_curb_rounding import FiniteSourceBank,fit_banks
from .generated_junction_curbs import GeneratedJunctionArea

SOURCE_BANK_POLICY='7-finite-local-crossing-of-unequal-outer-banks'


class SourceBankReconciliationError(ValueError):
    """A typed source-bank repair did not satisfy its original owner gates."""


class SourceBankExposureError(SourceBankReconciliationError):
    """Preserve actual finite source witnesses when an arc is not exposed."""
    def __init__(self,owner,proposal,fitted,spans,reason):
        self.owner,self.proposal,self.fitted=owner,proposal,fitted
        self.exposed_intervals=tuple((float(span.first),float(span.last)) for span in spans)
        first,last=proposal.first.token,proposal.last.token
        super().__init__(f'Source-bank fillet {owner.text()} {first} → {last}, '
            f'kind {proposal.kind}, radius {fitted.corner.radius:g} m: {reason}; '
            f'exposed intervals {self.exposed_intervals}')


def owned_ranges(context,group):
    """Actual disjoint native ownership; no hull across a returning middle."""
    values={claim.section_id:[(0.,context[Key('path',claim.section_id)].length)] for claim in group.claims}
    for arm,station in zip(group.physical_layout.approaches,group.physical_layout.outward_stations):
        if type(arm) is not Approach:continue
        if type(arm.path) is OwnedPrefixRoute:
            from .owned_prefix_corner_authority import verified_owned_prefix_ranges
            proof=verified_owned_prefix_ranges(context,group,arm,station)
            for sid,ranges in proof.ranges:
                values.setdefault(sid,[]).extend(ranges)
            continue
        if type(arm.path) is not Route:continue
        for (sid,end),at in arm.path.trims(station).items():
            length=context[Key('path',sid)].length
            first,last=(0.,at) if end==0 else (at,length)
            if last>first:values.setdefault(sid,[]).append((first,last))
    result={};count=0
    for sid,ranges in values.items():
        path=context[Key('path',sid)];merged=[]
        for first,last in sorted(ranges):
            if not (np.isfinite((first,last)).all() and 0.<=first<last<=path.length):
                raise SourceBankReconciliationError('Source bank escaped its original finite ownership')
            if merged and first<=merged[-1][1]:merged[-1]=(merged[-1][0],max(last,merged[-1][1]))
            else:merged.append((float(first),float(last)))
        count+=len(merged)
        if count>8192:raise SourceBankReconciliationError('Source banks exceed their finite ownership budget')
        result[sid]=tuple(merged)
    return result


@dataclass(frozen=True)
class BankOrigin:
    site: object
    endpoint: End
    outward_side: int


@dataclass(frozen=True)
class NativeBankWitness:
    path: object
    side: int
    station: float
    intervals: tuple
    origins: tuple
    owner: object
    index: int
    proof: str

    @property
    def token(self):return self.path.section.id,self.side


@dataclass(frozen=True)
class CircularCornerWitness:
    site: object
    first: NativeBankWitness
    last: NativeBankWitness
    corner: object
    owner: object
    index: int
    tangent: tuple
    entry_parameter: float
    exit_parameter: float


@dataclass(frozen=True)
class BankSwitch:
    loop: int
    vertex: int
    kind: str
    first: NativeBankWitness
    last: NativeBankWitness
    site: object
    radius: float
    mode: str
    sectors: tuple
    point: tuple


@dataclass(frozen=True)
class PreparedBankFillet(GeneratedJunctionArea):
    proposal: BankSwitch
    fitted: object


class SourceBankDiscovery:
    """One bounded scan of actual neighboring physical boundary owners."""
    def __init__(self,context,group):
        self.context,self.group=context,group
        self.sources={source.region.owner:source for source in group.footprints}
        self.sites={site.key:site for site in context.network.intersections}
        self.by_end={end:site.key for site in context.network.intersections for end in site.ends}
        self.owned=owned_ranges(context,group)
        self.internal=set()
        for claim in group.claims:
            sid=claim.section_id;path=context[Key('path',sid)]
            if (any(a==0. and b==path.length for a,b in self.owned[sid])
                    and all(self.by_end.get(End(sid,end)) in group.members for end in (0,1))):
                self.internal.update((End(sid,0),End(sid,1)))

    def origin(self,arm,side):
        endpoint=End(arm.source.id,arm.end);site=self.by_end.get(endpoint)
        if site not in self.group.members:return ()
        return (BankOrigin(site,endpoint,int(side)),)

    def witness(self,path,side,station,origins,owner,index,proof):
        if path is not self.context[Key('path',path.section.id)]:
            raise ValueError('Discovery changed an original native Path identity')
        if path.section.id not in self.owned:return None
        intervals=self.owned[path.section.id]
        if not any(a<=station<=b for a,b in intervals):return None
        if not origins:return None
        return NativeBankWitness(path,int(side),float(station),intervals,tuple(origins),owner,index,proof)

    def corner_bank(self,source,corner,first):
        arm=corner.first if first else corner.second
        station=corner.stations[0 if first else 1];outward_side=-1 if first else 1
        if type(arm) is not Approach or not 0.<=station<=arm.length:return None
        route=arm.path
        if type(route) is Route:
            i=int(np.clip(np.searchsorted(route.stations,station,side='right')-1,0,len(route.pieces)-1))
            piece=route.pieces[i];path=piece.path;local=station-route.stations[i]
            direction=1 if piece.forward else -1
            native=local if piece.forward else path.length-local
        elif type(route) is Path:
            path=route;direction=1 if arm.end==0 else -1
            native=station if arm.end==0 else path.length-station
        else:return None
        if (path.section.width!=arm.width or path.section.offset!=direction*arm.offset):return None
        origins=self.origin(arm,outward_side)
        site=getattr(source.layout,'original_junction',None)
        if not origins or origins[0].site!=site:return None
        return self.witness(path,direction*outward_side,native,origins,
            source.region.owner,-1,'EXPLICIT_ORIGINAL_CORNER_TANGENT_BANK')

    def direct(self,source,curve,label,tail,curve_intervals=None):
        coordinates=getattr(curve,'source_coordinates',None)
        pieces=[]
        if type(coordinates) is CurbCoordinates:
            intervals = None if curve_intervals is None else curve_intervals.get(id(curve))
            if intervals is None:
                coordinates.validate(curve.points,self.context)
                knots = curve.stations / curve.length
                # Match the actual floating-point expression for the end of
                # an edge. At fraction 1 it can differ from tb by one ulp.
                ends = knots[:-1] + np.diff(knots)
                intervals = knots, ends
                if curve_intervals is not None:
                    curve_intervals[id(curve)] = intervals
            knots, ends = intervals
            for edge in np.flatnonzero((knots[:-1] < label.last) & (ends > label.first)):
                ta,tb=knots[edge:edge+2]
                for piece in coordinates.edge_sources(edge):
                    low=ta+piece.first_fraction*(tb-ta);high=ta+piece.last_fraction*(tb-ta)
                    a,b=max(low,label.first),min(high,label.last)
                    if b<=a or piece.authored_side is None:continue
                    t=b if tail else a
                    station=piece.first_station+(t-low)/(high-low)*(piece.last_station-piece.first_station)
                    direction=piece.last_station-piece.first_station
                    if direction*piece.authored_side<=0.:continue
                    if coordinates.approach is not None:
                        origins=self.origin(coordinates.approach,int(coordinates.sides[edge]))
                    else:
                        origins=tuple(BankOrigin(self.by_end[End(piece.path.section.id,end)],End(piece.path.section.id,end),
                            piece.authored_side*(1 if end==0 else -1)) for end in (0,1)
                            if self.by_end.get(End(piece.path.section.id,end)) in self.group.members)
                    witness=self.witness(piece.path,piece.authored_side,station,origins,
                        source.region.owner,label.index,'EXPLICIT_PREPARED_NATIVE_VERTEX_PARAMETERS')
                    if witness is not None:pieces.append((a,b,witness))
        if pieces:
            piece=max(pieces,key=lambda row:row[1]) if tail else min(pieces,key=lambda row:row[0])
            target=label.last if tail else label.first
            if (piece[1] if tail else piece[0])==target:return piece[2]
        # An analytical composite can retain a finite OffsetCurb bank directly.
        target=label.last if tail else label.first
        for row in _intervals(curve,label.first,label.last,.5*(label.first+label.last)):
            a,b=sorted((curve_parameter(row,row['low']),curve_parameter(row,row['high'])))
            if not a-1e-12<=target<=b+1e-12:continue
            section=row['path'].section
            sides=[side for side in (-1,1) if row['offset']==section.offset+.5*section.width*side]
            if len(sides)!=1:continue
            corner=next(iter(getattr(source.layout,'corners',())),None)
            if corner is None:continue
            origins=[]
            for first in (True,False):
                bank=self.corner_bank(source,corner,first)
                if bank is not None and bank.path is row['path'] and bank.side==sides[0]:origins.extend(bank.origins)
            return self.witness(row['path'],sides[0],source_station(row,target),origins,
                source.region.owner,label.index,'EXPLICIT_FINITE_ANALYTICAL_OFFSET_MAP')
        return None

    @staticmethod
    def circular_tangent(curve,corner,parameter):
        """Use a declared analytical arc phase, never a polygon chord turn."""
        if corner.center is None:return None
        if hasattr(curve,'points'):
            points=np.asarray(curve.points);count=len(points)
            angles=np.linspace(corner.start_angle,corner.start_angle+corner.sweep,count)
            expected=corner.center+corner.radius*np.column_stack((np.cos(angles),np.sin(angles)))
            if np.max(np.linalg.norm(points-expected,axis=1))>1e-8:return None
            angle=float(np.interp(parameter,curve.stations/curve.length,angles))
            entry,exit=0.,1.
        elif getattr(curve,'corner',None) is corner and len(getattr(curve,'lengths',()))==4:
            entry,exit=map(float,curve.lengths[1:3]/curve.length)
            if not entry-1e-12<=parameter<=exit+1e-12 or exit<=entry:return None
            angle=corner.start_angle+(parameter-entry)/(exit-entry)*corner.sweep
        else:return None
        sign=1. if corner.sweep>0. else -1.
        return (sign*-math.sin(angle),sign*math.cos(angle)),entry,exit

    def describe(self,span,tail,curve_intervals=None):
        source=self.sources[span.owner];label=rim_label(span,self.sources)
        if label.kind!='curb':return None
        curve=source.layout.boundaries[label.index]
        bank=self.direct(source,curve,label,tail,curve_intervals)
        if bank is not None:return bank
        site=getattr(source.layout,'original_junction',None)
        corners=getattr(source.layout,'corners',())
        if site not in self.group.members or len(corners)!=1:return None
        generated=getattr(source.layout,'source_bank_fillet',None)
        if type(generated) is PreparedBankFillet:
            p=generated.proposal;fit=generated.fitted
            banks=[]
            for bank,arm,station in zip((p.first,p.last),(fit.corner.first,fit.corner.second),fit.corner.stations):
                banks.append(replace(bank,station=station if arm.end==0 else arm.length-station))
            tangent=self.circular_tangent(curve,fit.corner,label.last if tail else label.first)
            if tangent is None:return None
            return CircularCornerWitness(p.site,*banks,fit.corner,span.owner,label.index,*tangent)
        corner,=corners
        if corner.center is None:return None  # Polynomial envelopes have a separate ancestry contract.
        first,last=self.corner_bank(source,corner,True),self.corner_bank(source,corner,False)
        if first is None or last is None:return None
        tangent=self.circular_tangent(curve,corner,label.last if tail else label.first)
        if tangent is None:return None
        return CircularCornerWitness(site,first,last,corner,span.owner,label.index,*tangent)

    @staticmethod
    def same_bank(a,b):return a.path is b.path and a.side==b.side

    def tangent(self,bank):
        low,high=next((a,b) for a,b in bank.intervals if a<=bank.station<=b)
        step=min(1e-4,.001*(high-low));lo=max(low,bank.station-step);hi=min(high,bank.station+step)
        offset=bank.path.section.offset+.5*bank.path.section.width*bank.side
        tangent=bank.side*(bank.path.point(hi,offset)-bank.path.point(lo,offset))/(hi-lo)
        return tangent/np.linalg.norm(tangent)

    @staticmethod
    def tangent_continuous(a,b):
        return (float(a@b)>=1.-1e-10 and abs(float(a[0]*b[1]-a[1]*b[0]))<=1e-5)

    def continuous(self,a,b):
        if self.same_bank(a,b):return True
        tangents=[];heights=[];slopes=[]
        for bank in (a,b):
            tangents.append(self.tangent(bank))
            profile=self.context[Key('profile',bank.path.section.id)]
            arm=Approach(bank.path,0,profile)
            heights.append(profile.datum+float(profile.value(bank.side)))
            slopes.append(float(np.linalg.norm(arm.physical_gradient(bank.station,bank.side))))
        return (abs(heights[0]-heights[1])<=1e-7 and abs(slopes[0]-slopes[1])<=1e-7
            and self.tangent_continuous(*tangents))

    def _incident_bank_crossing(self, a, b, first, last):
        """Confirm the actual crossing of the two original endpoint offsets."""
        centers, starts, tangents, points, reaches = [], [], [], [], []
        tolerance = min(.001, .05*self.context.network.settings.curve_tolerance)
        local_reach = .5*(a.path.section.width+b.path.section.width)
        for bank, origin in ((a, first), (b, last)):
            path, end = bank.path, origin.endpoint.end
            section = path.section
            cuts = set(section.scope_ends) | set(section.selection_ends)
            if end in cuts:
                return False
            root = path.length if end else 0.
            far = path.length-root
            distance = abs(bank.station-root)
            if distance > local_reach:
                return False
            # A distant import cut can coexist with a complete junction
            # prefix. The prefix must stop before that cut.
            if not any(lo <= root <= hi and lo < bank.station < hi and
                       (1-end not in cuts or far < lo or far > hi)
                       for lo,hi in bank.intervals):
                return False
            offset = section.offset + .5*section.width*bank.side
            center,tangent,_ = path.frame(root)
            centers.append(center)
            starts.append(path.point(root,offset))
            tangents.append(tangent*(1-2*end))
            points.append(path.point(bank.station,offset))
            reaches.append(distance)
        matrix = np.column_stack((tangents[0],-tangents[1]))
        if np.linalg.norm(centers[0]-centers[1]) > 1e-7 or abs(np.linalg.det(matrix)) < 1e-5:
            return False
        stations = np.linalg.solve(matrix, starts[1]-starts[0])
        if np.any(stations < 0.) or np.max(np.abs(stations-reaches)) > tolerance:
            return False
        crossing = starts[0]+stations[0]*tangents[0]
        return all(np.linalg.norm(point-crossing) <= tolerance for point in points)

    def _exposed_profile_join(self, a, b, first, last, walk, site):
        """Two exposed outer banks can meet behind one live side entrance."""
        if len(walk) != 2 or walk[1] in self.internal:
            return False
        if any(site.corner_radius(end)[1] != 'STYLE' for end in walk):
            return False
        layout = self.group.physical_layout
        if not any(End(arm.source.id, arm.end) == walk[1] for arm in layout.approaches):
            return False
        profiles = []
        crossfalls = []
        needs_crossing_proof = False
        for bank, origin in ((a, first), (b, last)):
            section = bank.path.section
            if bank.proof != 'EXPLICIT_PREPARED_NATIVE_VERTEX_PARAMETERS':
                return False
            distance = bank.station if origin.endpoint.end == 0 else bank.path.length-bank.station
            if section.scope_ends or section.selection_ends or distance > .5*section.width:
                needs_crossing_proof = True
            profile = self.context[Key('profile', section.id)]
            if profile.value(bank.side) != 0.:
                return False
            profiles.append(profile)
            crossfalls.append(abs(float(profile.derivative(bank.side))) * 2. / section.width)
        if (bool(profiles[0].half) == bool(profiles[1].half)
                or abs(profiles[0].datum-profiles[1].datum) > 1e-10
                or abs(crossfalls[0]-crossfalls[1]) <= 1e-7):
            return False
        incoming, outgoing = self.tangent(a), self.tangent(b)
        if incoming[0]*outgoing[1]-incoming[1]*outgoing[0] >= -1e-5:
            return False
        if not needs_crossing_proof and incoming @ outgoing > 0.:
            return True
        return self._incident_bank_crossing(a,b,first,last)

    def policy(self,a,b,kind,*,exterior=False):
        choices=[]
        for first in a.origins:
            for last in b.origins:
                if first.site!=last.site or first.endpoint==last.endpoint:continue
                if (first.site not in self.group.members or first.site not in self.sites
                        or self.by_end.get(first.endpoint)!=first.site
                        or self.by_end.get(last.endpoint)!=last.site):continue
                site=self.sites[first.site]
                # This is the same angular order as the source-cycle producer,
                # over real original incidences, not source-coordinate matching.
                def angle(end):
                    path=self.context[Key('path',end.section_id)]
                    tangent=path.frame(0. if end.end==0 else path.length)[1]*(1 if end.end==0 else -1)
                    return math.atan2(tangent[1],tangent[0])%math.tau
                ring=tuple(sorted(site.ends,key=angle));i,j=ring.index(first.endpoint),ring.index(last.endpoint)
                sectors=[]
                if ring[(i+1)%len(ring)]==last.endpoint:
                    sectors=[first.endpoint]
                elif ring[(j+1)%len(ring)]==first.endpoint:
                    sectors=[last.endpoint]
                elif first.outward_side==-1 and last.outward_side==1:
                    walk=[];at=i
                    while at!=j:
                        walk.append(ring[at]);at=(at+1)%len(ring)
                    if (all(end in self.internal for end in walk[1:])
                            or exterior and kind=='DIRECT_ORIGINAL_BANK_SWITCH'
                            and self._exposed_profile_join(a,b,first,last,walk,site)):
                        sectors=walk
                if not sectors:continue
                records=tuple((end,*site.corner_radius(end)) for end in sectors)
                strict={radius for end,radius,mode in records if mode=='STRICT'}
                style=[radius for end,radius,mode in records if mode=='STYLE']
                if len(strict)>1:continue
                if strict:
                    radius=next(iter(strict));mode='STRICT'
                    if style and radius>min(style):continue
                else:radius=min(style);mode='STYLE'
                choices.append((first.site,radius,mode,records))
        unique={(site,radius,mode,sectors) for site,radius,mode,sectors in choices}
        return next(iter(unique)) if len(unique)==1 else None

    def natural_separator_tips(self, loop, switches):
        """Recognize a closed plateau between two original Bundle banks."""
        from .corner_radii import CornerContactError, CornerRadiusPolicy
        from .profiles import memberships
        if (not loop.hole or len(switches) != 2 or
                not all(p.mode == 'STYLE' and p.kind == 'DIRECT_ORIGINAL_BANK_SWITCH' for p in switches)):
            return ()
        start, end = (np.asarray(p.point) for p in switches)
        axis = end - start
        length = np.linalg.norm(axis)
        if length == 0.:
            return ()
        # The projection bounds the width of the entire loop, including
        # between samples because the actual rim consists of straight chords.
        width = np.ptp((loop.xy - start) @ (np.array((-axis[1], axis[0])) / length))
        if width >= 2. * min(CornerRadiusPolicy.minimum_gore_radius, *(p.radius for p in switches)):
            return ()
        proofs = getattr(self.group.physical_layout, 'source_cycle_proofs', {})
        if not proofs:
            return ()
        owned = memberships(self.context.network)
        bundle_id = None
        roles = []
        plateau = None
        for span in loop.spans:
            first, last = (self.describe(span, tail) for tail in (False, True))
            if not (isinstance(first, NativeBankWitness) and isinstance(last, NativeBankWitness)
                    and self.same_bank(first, last)):
                return ()
            sid = first.path.section.id
            if sid not in owned:
                return ()
            bundle, role, forward = owned[sid]
            side = (1 if role == 0 else -1) * (1 if forward else -1)
            profile = self.context[Key('profile', sid)]
            if (bundle.id not in proofs or len(bundle.members) != 2 or not all(m.complete for m in bundle.members)
                    or first.side != side or profile.half != side or profile.derivative(side) != 0.):
                return ()
            level = profile.datum + profile.value(side)
            if plateau is not None and abs(level - plateau) > 1e-10:
                return ()
            plateau = level
            if bundle_id is not None and bundle.id != bundle_id:
                return ()
            bundle_id = bundle.id
            roles.append(role)
        changes = {i for i, role in enumerate(roles) if role != roles[i-1]}
        if changes != {p.vertex for p in switches}:
            return ()
        tips = []
        for p in switches:
            try:
                fit_exposed_switch(self.context, self.group, p)
            except CornerContactError:
                tips.append(p.vertex)
        return tips

    def sampled_street_continuation(self, proposal):
        """A fitted style arc can be indistinguishable from one original Street seam."""
        from .corner_radii import CornerContactError
        from .regions import RegionUnion
        if proposal.kind != 'DIRECT_ORIGINAL_BANK_SWITCH' or proposal.mode != 'STYLE':
            return False
        a, b = proposal.first, proposal.last
        sections = a.path.section, b.path.section
        if (sections[0].street_id != sections[1].street_id or sections[0].id == sections[1].id
                or sections[0].width != sections[1].width
                or a.side*sections[0].offset != b.side*sections[1].offset):
            return False
        endpoints = []
        profiles = []
        # This bounds the diameter of the entire junction edit in XY,
        # including longitudinal displacement along a nearly straight bank.
        budget = min(.001, self.context.network.settings.curve_tolerance)
        for bank in (a,b):
            origins = [o.endpoint.end for o in bank.origins if o.site == proposal.site]
            if len(origins) != 1 or bank.proof != 'EXPLICIT_PREPARED_NATIVE_VERTEX_PARAMETERS':
                return False
            end, = origins
            section = bank.path.section
            if (end in section.scope_ends or end in section.selection_ends
                    or abs(bank.station - end*bank.path.length) > budget):
                return False
            endpoints.append(section.points[-1 if end else 0])
            profiles.append(self.context[Key('profile', section.id)].oriented(int(bank.side < 0)))
        if endpoints[0] != endpoints[1] or profiles[0] != profiles[1]:
            return False
        # The two tangents point along the same traversed rim. This estimate
        # only avoids unnecessary fits; the actual fitted circle decides.
        ta,tb = self.tangent(a),self.tangent(b)
        if ta @ tb <= 0. or abs(ta[0]*tb[1]-ta[1]*tb[0])*proposal.radius > budget:
            return False
        try:
            fit = fit_exposed_switch(self.context,self.group,proposal)
        except CornerContactError:
            return False
        corner = fit.corner
        sagitta = 2.*corner.radius*math.sin(.25*corner.sweep)**2
        if (corner.radius != proposal.radius or abs(corner.sweep) >= math.pi/2.
                or sagitta > RegionUnion().epsilon):
            return False
        # Bound the entire circular arc and both trimmed bank runs relative
        # to the emitted seam. The disk contains the endpoint chord by
        # convexity. Expanding it by the sagitta also contains the exact arc.
        pivot = np.asarray(proposal.point)
        normal = np.array((-ta[1],ta[0]))
        transverse_budget = min(budget, .05*self.context.network.settings.curve_tolerance)
        for arm,side,cut,contact in zip((corner.first,corner.second),(-1,1),fit.intersection,corner.stations):
            if type(arm.path) is not Path or arm.path.exact:
                return False
            offset = arm.offset + side*.5*arm.width
            crossing = arm.point(cut,offset)
            point = arm.point(contact,offset)
            curvature = arm.path.curvature_bound
            if abs(offset)*curvature >= 1.:
                return False
            length = abs(contact-cut)*(1.+abs(offset)*curvature)
            bank_error = curvature/(1.-abs(offset)*curvature)*length**2/8.
            if (np.linalg.norm(crossing-pivot)+length > .5*budget
                    or np.linalg.norm(point-pivot)+sagitta > .5*budget
                    or max(abs((crossing-pivot) @ normal),abs((point-pivot) @ normal))
                    + max(bank_error,sagitta) > .5*transverse_budget):
                return False
        return True

    def discover(self):
        proposals=[]
        for li,loop in enumerate(self.group.loops):
            start = len(proposals)
            for index,last_span in enumerate(loop.spans):
                first_span=loop.spans[index-1]
                if first_span.owner==last_span.owner:continue
                first,last=self.describe(first_span,True),self.describe(last_span,False)
                kind=None
                if isinstance(first,NativeBankWitness) and isinstance(last,NativeBankWitness):
                    if self.continuous(first,last):continue
                    a,b=first,last;kind='DIRECT_ORIGINAL_BANK_SWITCH'
                elif isinstance(first,CircularCornerWitness) and isinstance(last,NativeBankWitness):
                    if self.same_bank(first.last,last):continue
                    if self.tangent_continuous(np.array(first.tangent),self.tangent(last)):continue
                    if not self.visible_anchor(loop,index-1,first,entry=True):continue
                    a,b=first.first,last;kind='CLIPPED_ORIGINAL_CORNER_ENTRY_TO_BANK'
                elif isinstance(first,NativeBankWitness) and isinstance(last,CircularCornerWitness):
                    if self.same_bank(first,last.first):continue
                    if self.tangent_continuous(self.tangent(first),np.array(last.tangent)):continue
                    if not self.visible_anchor(loop,index,last,entry=False):continue
                    a,b=first,last.last;kind='BANK_TO_CLIPPED_ORIGINAL_CORNER_EXIT'
                else:continue
                if self.same_bank(a,b):continue
                policy=self.policy(a,b,kind,exterior=not loop.hole)
                if policy is None:continue
                site,radius,mode,sectors=policy
                proposal=BankSwitch(li,index,kind,a,b,site,radius,mode,sectors,tuple(loop.xy[index]))
                if self.sampled_street_continuation(proposal):
                    continue
                proposals.append(proposal)
                if len(proposals)>64:raise ValueError('Source-bank discovery exceeds its finite transition budget')
            switches = proposals[start:]
            tips = self.natural_separator_tips(loop, switches)
            if tips:
                proposals[start:] = [p for p in switches if p.vertex not in tips]
        return tuple(proposals)

    def visible_anchor(self,loop,index,corner,*,entry):
        spans=loop.spans;step=-1 if entry else 1;at=index%len(spans)
        for _ in range(len(spans)):
            target=(at+step)%len(spans);span=spans[target];label=rim_label(span,self.sources)
            if (span.owner,label.kind,label.index)!=(corner.owner,'curb',corner.index):
                boundary=rim_label(spans[at],self.sources)
                retained=boundary.first if entry else boundary.last
                expected=corner.entry_parameter if entry else corner.exit_parameter
                if abs(retained-expected)>1e-9:return False
                witness=self.describe(span,entry)
                return (isinstance(witness,NativeBankWitness)
                    and self.same_bank(witness,corner.first if entry else corner.last))
            at=target
        return False


def visible_bank_stations(discovery, proposal, *, incoming):
    """Follow one directed Section bank up to its next change of source."""
    spans = discovery.group.loops[proposal.loop].spans
    bank = proposal.first if incoming else proposal.last
    direction = -1 if incoming else 1
    at = (proposal.vertex-1 if incoming else proposal.vertex) % len(spans)
    # One synchronous walk over a fixed physical footprint. Revalidate on
    # the next walk, including after any refinement of its original curves.
    curve_intervals = {}
    for _ in range(len(spans)):
        adjacent = (at+direction) % len(spans)
        witness = discovery.describe(spans[adjacent],incoming,curve_intervals)
        if not isinstance(witness,NativeBankWitness) or not discovery.same_bank(witness,bank):
            break
        at = adjacent
    else:
        raise SourceBankReconciliationError('One original bank occupies the entire junction boundary')
    witness = discovery.describe(spans[at],not incoming,curve_intervals)
    if not isinstance(witness,NativeBankWitness) or not discovery.same_bank(witness,bank):
        return None
    start, stop = (witness.station,bank.station) if incoming else (bank.station,witness.station)
    low,high = sorted((start,stop))
    if bank.side*(stop-start) <= 0. or not any(a <= low < high <= b for a,b in bank.intervals):
        return None
    return start,stop


def fit_switch(context,proposal):
    return fit_banks(context,FiniteSourceBank(proposal.first.path,proposal.first.side,proposal.first.intervals),
        FiniteSourceBank(proposal.last.path,proposal.last.side,proposal.last.intervals),
        proposal.radius,proposal.mode,np.array(proposal.point))


def fit_exposed_switch(context,group,proposal):
    from .corner_radii import CornerContactError
    try:
        return fit_switch(context,proposal)
    except CornerContactError:
        if proposal.kind != 'DIRECT_ORIGINAL_BANK_SWITCH':
            raise
        from .local_bank_contacts import fit
        return fit(context,group,proposal)


def is_generated_bank_area(area):
    return type(getattr(area.layout,'source_bank_fillet',None)) is PreparedBankFillet


def signature(group):
    record=getattr(group.physical_layout,'source_bank_reconciliation',None)
    handoffs=getattr(group.physical_layout,'source_bank_handoff_canonicalizations',())
    fillets = getattr(group.physical_layout,'native_original_circle_fillets',())
    if fillets:
        return (SOURCE_BANK_POLICY, tuple(record['records']), tuple(f.report() for f in fillets),
                tuple(handoff['construction'] for handoff in handoffs))
    # Unrelated scene owners can change exclusion counters without changing
    # this geometry. The local original bank and sector records carry its
    # actual construction inputs; the final PLAN guard must still pass.
    if handoffs:
        return (SOURCE_BANK_POLICY,tuple(record['records']) if record else (),
                tuple(handoff['construction'] for handoff in handoffs))
    return (SOURCE_BANK_POLICY,tuple(record['records'])) if record else None


def reconcile(planner,group,partitions):
    from .updates import drain_steps
    return drain_steps(reconcile_steps(planner,group,partitions))


def reconcile_steps(planner,group,partitions):
    from .source_bank_preparation import SourceBankReconciler
    reconciler = SourceBankReconciler(planner,group,partitions)
    try:
        return (yield from reconciler.build_steps())
    except SourceBankExposureError as failure:
        if failure.proposal.kind != 'BANK_TO_CLIPPED_ORIGINAL_CORNER_EXIT':
            raise
        from .bank_circle_preparation import recover_steps
        return (yield from recover_steps(reconciler, failure))


def validate_prepared_neighbors(context,groups,layouts):
    from .source_bank_neighbors import validate
    return validate(context,groups,layouts)
