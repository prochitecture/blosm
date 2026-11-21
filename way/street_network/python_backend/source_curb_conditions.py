"""Boundary jets from exposed original source curbs.

A closed road block is not necessarily a median plateau. Original strip and
port curbs supply physical height/crossfall anchors. Clipped corner fragments
join the surviving anchors; a swallowed internal contact adds no phantom seam.
The optional caller retains the original profiles, XY, radii and source spans.
"""
from dataclasses import replace
from types import SimpleNamespace
import math
import numpy as np
from .compounds import rim_label
from .conditions import CurbEdge
from .model import Key
from .paths import right


class SourceCurbConditionUnavailable(ValueError):
    """A physical curb lacks the explicit source provenance for this policy."""


class _CurbJoinFailure(ValueError):
    def at_join(self, spans, points, index):
        self.join=(spans[index-1].owner,spans[index].owner,tuple(points[index]))
        return self


class SourceAnchorConflict(_CurbJoinFailure):
    """Two retained original boundary anchors impose incompatible values."""


class SourceJetDirectionError(_CurbJoinFailure):
    """A retained source jet points outside its generated corner chain."""


def _needs_exact_station_ledger(path):
    paths = tuple(piece.path for piece in path.pieces) if hasattr(path, 'pieces') else (path,)
    return any(value.exact and np.max(abs(np.linalg.norm(value.normals, axis=1)-1.)) > 1e-10
               for value in paths)


def _corner_slopes(corner):
    return tuple(float(np.linalg.norm(approach.physical_gradient(max(0., station), q)))
                 for approach, station, q in zip((corner.first, corner.second), corner.stations, (-1., 1.)))


class _Anchor:
    def __init__(self, height, gradient):
        self.height, self.gradient = float(height), np.asarray(gradient)

    def sample(self, t):
        return self.height, self.gradient


class _SourceRule:
    def __init__(self, curve, first, last, z0, z1, k0, k1, phase=(0., 1.), *, corner_fragment=False):
        self.first, self.last = float(first), float(last)
        self.z0, self.z1, self.k0, self.k1 = map(float, (z0, z1, k0, k1))
        self.phase = tuple(map(float, phase))
        self.corner_fragment = bool(corner_fragment)
        self.source_anchor = (not corner_fragment or last <= phase[0]+1e-10 or first >= phase[1]-1e-10
            or abs(z0-z1) <= 1e-10 and abs(k0-k1) <= 1e-10)
        if (not all(map(math.isfinite, (*phase, first, last, z0, z1, k0, k1)))
                or not 0. <= first < last <= 1.+1e-9
                or not 0. <= phase[0] < phase[1] <= 1.+1e-9 or min(k0, k1) < 0.):
            raise ValueError('Source curb condition has invalid finite height/jet parameters')
        self.length = float(curve.length)
        points = np.asarray(curve(np.asarray((first, last))))
        tangent = points[1]-points[0]
        if not self.length > 0. or not np.isfinite(points).all() or np.linalg.norm(tangent) <= 1e-10:
            raise ValueError('Source curb condition lost its finite physical chord')
        self.tangent = tangent/np.linalg.norm(tangent)

    def polynomial_breaks(self):
        # Clamping the quintic at the contacts creates a constant branch
        # on each lead, even when a source chord straddles the contact.
        return tuple((contact-self.first)/(self.last-self.first) for contact in self.phase
                     if self.first < contact < self.last)

    def sample(self, t):
        at = self.first+np.asarray(t)*(self.last-self.first)
        first, last = self.phase
        v = np.clip((at-first)/(last-first), 0., 1.)
        weight = v**3*(10.+v*(-15.+6.*v))
        height = self.z0+weight*(self.z1-self.z0)
        inward = self.k0+weight*(self.k1-self.k0)
        along = (self.z1-self.z0)*30.*v*v*(1.-v)**2/((last-first)*self.length)
        return height, along[..., None]*self.tangent-inward[..., None]*right(self.tangent)


class _FieldRule:
    """An explicit prepared merged entrance already owns its physical field."""
    def __init__(self, curve, first, last, field, *, original=None):
        self.curve, self.first, self.last, self.field = curve, first, last, field
        self.original = original

    def native_values(self, parameters):
        if self.original is None:
            raise SourceCurbConditionUnavailable('Merged source curb needs its original bank ancestry outside the common field')
        context, approach, bank = self.original
        path = approach.paths[bank.role]
        original = approach.original_approaches[bank.role]
        if original.source.id != bank.section or path.section.id != bank.section:
            raise ValueError('Terminal source-bank field changed its original Section identity')
        at = bank.first+parameters*(bank.last-bank.first)
        outward = at if original.end == 0 else path.length-at
        q = bank.side*(1 if original.end == 0 else -1)
        profile = context[Key('profile', bank.section)]
        return (np.full(len(at),profile.datum+profile.value(bank.side)),
                np.asarray([original.physical_gradient(float(value), q) for value in outward]))

    def sample(self, t):
        t = np.asarray(t)
        parameters = self.first+t*(self.last-self.first)
        xy = self.curve(parameters).reshape(-1, 2)
        station = self.field.frame.local(xy)[:, 0]
        low = max(float(curve.stations[0]) for curve in self.field.curves)
        high = min(float(curve.stations[-1]) for curve in self.field.curves)
        inside = (station >= low) & (station <= high)
        height, gradient = np.empty(len(xy)), np.empty((len(xy), 2))
        if inside.any():
            height[inside] = self.field.values(xy[inside])
            gradient[inside] = self.field.gradient(xy[inside])
        if (~inside).any():
            height[~inside],gradient[~inside] = self.native_values(parameters.reshape(-1)[~inside])
        return height.reshape(t.shape), gradient.reshape(t.shape+(2,))


class _LinearTerminalRule(_FieldRule):
    """One-sided gradients along a straight clipped edge of the original paired banks."""
    def __init__(self, curve, first, last, field, *, original):
        super().__init__(curve, first, last, field, original=original)
        a,b=field.frame.local(curve(np.array((first,last))))[:,0]
        self.common=(max(ray.stations[0] for ray in field.curves),
                     min(ray.stations[-1] for ray in field.curves))
        self.stations={0.:a,1.:b}
        low,high=sorted((a,b))
        for ray in field.curves:
            i,j=np.searchsorted(ray.stations,(low,high))
            for station in ray.stations[i:j]:
                if low<station<high and self.common[0]<=station<=self.common[1]:
                    self.stations[(station-a)/(b-a)]=station
        self.bounds=tuple(sorted(self.stations))

    def polynomial_breaks(self):
        return self.bounds[1:-1]

    def sample_piece(self,t,first,last):
        i=np.searchsorted(self.bounds,first)
        if i>=len(self.bounds)-1 or self.bounds[i:i+2]!=(first,last):
            raise ValueError('Linear terminal prediction must use one original chord interval')
        a,b=self.stations[first],self.stations[last]
        t=np.asarray(t)
        u=first+t*(last-first)
        if a==b:
            return super().sample(u)
        p=self.first+u*(self.last-self.first)
        if min(a,b)>=self.common[0] and max(a,b)<=self.common[1]:
            xy=self.curve(p).reshape(-1,2)
            h,g=self.field.values(xy),self.field.gradient(xy,chord_interval=(a,b))
        else:
            h,g=self.native_values(p.reshape(-1))
        return h.reshape(t.shape),g.reshape(t.shape+(2,))

    def sample(self,t):
        t=np.asarray(t)
        if len(self.bounds)==2:
            return self.sample_piece(t,0.,1.)
        parameters=t.reshape(-1)
        indices=np.clip(np.searchsorted(self.bounds,parameters,side='right')-1,0,len(self.bounds)-2)
        height,gradient=np.empty(len(parameters)),np.empty((len(parameters),2))
        for i in np.unique(indices):
            selected=indices==i
            a,b=self.bounds[i:i+2]
            height[selected],gradient[selected]=self.sample_piece((parameters[selected]-a)/(b-a),a,b)
        return height.reshape(t.shape),gradient.reshape(t.shape+(2,))


class _OffsetRule:
    """The declared finite original approach/side supplies its own exact jet."""
    def __init__(self, curve, first, last):
        self.curve, self.first, self.last = curve, first, last
        self.q = 2.*(curve.offset-curve.approach.offset)/curve.approach.width
        if abs(abs(self.q)-1.) > 1e-8:
            raise ValueError('Retained source offset condition is not an actual carriageway curb')

    def sample(self, t):
        t = self.first+np.asarray(t)*(self.last-self.first)
        station = np.interp(t*self.curve.length, self.curve.cumulative, self.curve.parameters)
        approach = self.curve.approach
        height = np.full(np.asarray(t).shape, approach.profile.datum+approach.profile.value(self.q))
        gradient = np.asarray([approach.physical_gradient(float(value), self.q)
                               for value in np.asarray(station).reshape(-1)])
        return height, gradient.reshape(np.asarray(t).shape+(2,))


def _native_gradient(path,profile,station,side,*,interval=None):
    """Exact native profile covector, with an explicit one-sided Exact domain.

    A Smooth offset has tangent parallel to its centreline and a unit normal.
    Exact normals interpolate the authored miter vectors; their derivative is
    constant on each original segment. A stencil across a control knot would
    mix different source Jacobians rather than evaluate either source piece.
    """
    at=np.asarray(station)
    normal=path.frame(at)[2]
    if not path.exact:
        return profile.derivative(side)*normal/(.5*path.section.width)
    if interval is None:
        index=np.clip(np.searchsorted(path.knots,at,side='right')-1,0,len(path.chords)-1)
    else:
        low,high=sorted(map(float,interval))
        index=int(np.clip(np.searchsorted(path.knots,.5*(low+high),side='right')-1,
                          0,len(path.chords)-1))
        tolerance=8.*np.spacing(max(1.,path.length))
        if low<path.knots[index]-tolerance or high>path.knots[index+1]+tolerance:
            raise SourceCurbConditionUnavailable('Exact native jet crosses an unretained original segment')
    offset=path.section.offset+.5*path.section.width*side
    velocity=(path.directions[index]+offset*(path.normals[index+1]-path.normals[index])/
              np.asarray(path.chords[index])[...,None])
    covector=right(velocity)
    denominator=np.sum(covector*(.5*path.section.width*normal),axis=-1)
    if np.any(abs(denominator)<=1e-10):
        raise ValueError('Original native bank has a singular physical profile frame')
    return profile.derivative(side)*covector/denominator[...,None]


class _NativeCoordinateRule:
    """Actual bank sides/stations from an explicit physical vertex ledger.

    A closed strip is split into convenient curve intervals. Either interval
    can include parts of both banks, so its curve index is not a bank identity.
    """
    corner_fragment=False
    source_anchor=True

    def __init__(self,context,curve,first,last,*,source_knots=None):
        from .source_curb_coordinates import CurbCoordinates
        coordinates=getattr(curve,'source_coordinates',None)
        if type(coordinates) is not CurbCoordinates:
            raise SourceCurbConditionUnavailable('Native curb needs its explicit original coordinate ledger')
        coordinates.validate(curve.points,context)
        self.first,self.last=float(first),float(last)
        if (not np.isfinite((first,last)).all() or not 0.<=first<last<=1.+1e-9):
            raise ValueError('Native curb needs one finite increasing boundary interval')
        self.curve=curve
        low,high=np.asarray((first,last))*curve.length
        begin=max(0,int(np.searchsorted(curve.stations,low,side='right')-1))
        end=min(len(curve.stations)-1,int(np.searchsorted(curve.stations,high,side='left')))
        records=[]
        for edge in range(begin,end):
            a,b=curve.stations[edge:edge+2]
            for piece in coordinates.edge_sources(edge):
                p0=a+piece.first_fraction*(b-a);p1=a+piece.last_fraction*(b-a)
                if min(high,p1)<=max(low,p0)+1e-12:continue
                side=piece.authored_side
                expected=piece.path.section.offset+.5*piece.path.section.width*(side or 0)
                if (side not in (-1,1) or abs(piece.first_offset-expected)>1e-8
                        or abs(piece.last_offset-expected)>1e-8):
                    raise SourceCurbConditionUnavailable('Physical curb interval crosses a source cap or non-bank offset')
                sid=piece.path.section.id
                if piece.path is not context[Key('path',sid)]:
                    raise ValueError('Native curb changed original Path identity')
                records.append((float(p0),float(p1),piece.path,context[Key('profile',sid)],
                                piece.first_station,piece.last_station,side))
        if not records or records[0][0]>low+1e-9 or records[-1][1]<high-1e-9:
            raise SourceCurbConditionUnavailable('Native curb interval has incomplete finite parameter coverage')
        if any(abs(a[1]-b[0])>1e-9 for a,b in zip(records,records[1:])):
            raise SourceCurbConditionUnavailable('Native curb interval skips a real source parameter range')
        levels=[profile.datum+float(profile.value(side)) for _,_,_,profile,_,_,side in records]
        if max(levels)-min(levels)>1e-10:
            raise SourceCurbConditionUnavailable('Adjacent native curb pieces have incompatible source heights')
        self.records=tuple(records)
        self.ends=np.array([record[1] for record in records])
        self.z0=self.z1=levels[0]
        cache={} if source_knots is None else source_knots
        breaks=[]
        # Multiplication by curve.length can recover a stored chord endpoint
        # a few ulps away. Only this representation roundoff is canonicalized;
        # no physical edge or original source control is removed or resampled.
        roundoff=8.*np.spacing(max(1.,curve.length))
        for p0,p1,path,_,s0,s1,_ in records:
            if s0==s1:
                raise SourceCurbConditionUnavailable('Native bank interval has no finite source progress')
            if path not in cache:
                knots=np.asarray(path.station(path.knots))
                knots.flags.writeable=False
                cache[path]=knots
            knots=cache[path]
            a,b=sorted((s0,s1))
            inside=knots[(knots>a)&(knots<b)]
            metric=p0+(inside-s0)*(p1-p0)/(s1-s0)
            for value in (*metric,p0,p1):
                if low+roundoff<value<high-roundoff:
                    breaks.append(float((value-low)/(high-low)))
        self._polynomial_breaks=tuple(sorted(set(breaks)))

    def polynomial_breaks(self):
        """Exact source-piece boundaries in this physical chord's local t."""
        return self._polynomial_breaks

    def sample_piece(self,t,first,last):
        """Sample one declared branch, including its one-sided endpoint jets."""
        boundaries=(0.,*self._polynomial_breaks,1.)
        if (first not in boundaries or last not in boundaries or not first<last
                or any(first<value<last for value in self._polynomial_breaks)):
            raise SourceCurbConditionUnavailable('Native polynomial branch must retain every declared source piece')
        t=np.asarray(t)
        if not np.isfinite(t).all() or np.any(t<0.) or np.any(t>1.):
            raise ValueError('Native polynomial samples need finite local branch parameters')
        metric=(self.first+(first+t*(last-first))*(self.last-self.first))*self.curve.length
        middle=(self.first+.5*(first+last)*(self.last-self.first))*self.curve.length
        index=int(np.clip(np.searchsorted(self.ends,middle,side='right'),0,len(self.records)-1))
        p0,p1,path,profile,s0,s1,side=self.records[index]
        tolerance=8.*np.spacing(max(1.,self.curve.length))
        if np.any(metric<p0-tolerance) or np.any(metric>p1+tolerance):
            raise SourceCurbConditionUnavailable('Native polynomial branch crosses an original Route seam')
        station=s0+np.clip((metric-p0)/(p1-p0),0.,1.)*(s1-s0)
        ends=(self.first+np.asarray((first,last))*(self.last-self.first))*self.curve.length
        interval=s0+np.clip((ends-p0)/(p1-p0),0.,1.)*(s1-s0)
        height=np.full(t.shape,profile.datum+profile.value(side))
        return height,_native_gradient(path,profile,station,side,interval=interval)

    def sample(self,t):
        if not self._polynomial_breaks:
            return self.sample_piece(t,0.,1.)
        shape=np.asarray(t).shape
        metric=(self.first+np.asarray(t)*(self.last-self.first)).reshape(-1)*self.curve.length
        indices=np.clip(np.searchsorted(self.ends,metric,side='right'),0,len(self.records)-1)
        heights=np.empty(len(metric));gradients=np.empty((len(metric),2))
        for i in np.unique(indices):
            selected=indices==i
            first,last,path,profile,s0,s1,side=self.records[int(i)]
            station=s0+np.clip((metric[selected]-first)/(last-first),0.,1.)*(s1-s0)
            heights[selected]=profile.datum+profile.value(side)
            gradients[selected]=_native_gradient(path,profile,station,side)
        return heights.reshape(shape),gradients.reshape(shape+(2,))


def _owned_connector_rule(resolver,source,curve,first,last,record):
    """Actual native contact jets, distinct from logical connector incidence."""
    from .source_bank_reconciliation import owned_ranges
    from .junction_bank_transitions import original_corner
    if (getattr(source.layout,'generated_junction_area',None) is not record
            or len(record.contacts)!=2 or len(record.source_contact_walks)!=2
            or record.proposal.site not in resolver.group.members
            or record.connector_walk.through_site not in resolver.group.members):
        raise SourceCurbConditionUnavailable('Owned connector lost its typed original source ancestry')
    from .successor_bank_contacts import SuccessorBankFillet
    original_site = record.proposal.site
    if isinstance(record,SuccessorBankFillet):
        original_site = record.connector_walk.through_site
        sites = {site.key:site for site in resolver.context.network.intersections}
        if any(key not in sites or end not in sites[key].ends or mode != 'STYLE'
               or sites[key].corner_radius(end) != (radius,mode)
               for key,end,radius,mode in record.connector_walk.sector_policies):
            raise SourceCurbConditionUnavailable('Successor bank lost its original corner radius policy')
    node=next((node for node in resolver.context.network.intersections if node.key==original_site),None)
    if node is None or not original_corner(resolver,node,record.original_corner.corner):
        raise SourceCurbConditionUnavailable('Owned connector changed its logical original radius/profile witness')
    claimed={claim.section_id for claim in resolver.group.claims}
    if (not record.connector_walk.pieces
            or any(piece.path is not resolver.context[Key('path',piece.section.id)]
                   or piece.section.id not in claimed for piece in record.connector_walk.pieces)):
        raise SourceCurbConditionUnavailable('Owned connector lacks its actual completely consumed original road')
    ranges=owned_ranges(resolver.context,resolver.group)
    corner=record.corner;heights=[];slopes=[];contacts=[]
    for contact,bank,arm,at,walk,facing in zip(record.contacts,
            (record.proposal.first,record.proposal.last),(corner.first,corner.second),
            corner.stations,record.source_contact_walks,(-1,1)):
        path=contact.path;sid=path.section.id
        profile=resolver.context[Key('profile',sid)]
        native=float(at if arm.end==0 else path.length-at)
        q=contact.side*(1 if arm.end==0 else -1)
        if (path is not resolver.context[Key('path',sid)] or bank.path is not path
                or bank.side!=contact.side or contact.side not in (-1,1) or q!=facing
                or arm.path is not path or arm.source is not path.section
                or arm.width!=path.section.width
                or arm.offset!=path.section.offset*(1 if arm.end==0 else -1)
                or arm.profile!=profile.oriented(arm.end) or walk.approach is not arm
                or not np.isfinite(contact.station) or not 0.<=contact.station<=path.length
                or abs(native-contact.station)>1e-9
                or not contact.owned_intervals or contact.owned_intervals!=ranges.get(sid)
                or not any(a<=contact.station<=b for a,b in contact.owned_intervals)):
            raise SourceCurbConditionUnavailable('Owned connector contact lacks its sealed finite native profile domain')
        heights.append(profile.datum+float(profile.value(contact.side)))
        slopes.append(float(np.linalg.norm(_native_gradient(path,profile,contact.station,contact.side))))
        contacts.append(path.point(contact.station,path.section.offset+.5*path.section.width*contact.side))
    if not hasattr(curve,'points'):
        raise SourceCurbConditionUnavailable('Owned connector arc needs its explicit original contact knots')
    points=np.asarray(curve.points)
    stations=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
    ids=[int(np.argmin(np.linalg.norm(points-point,axis=1))) for point in contacts]
    if ids[0]==ids[1] or max(np.linalg.norm(points[i]-point) for i,point in zip(ids,contacts))>1e-7:
        raise SourceCurbConditionUnavailable('Owned connector boundary lost its actual native contacts')
    if ids[0]>ids[1]:ids.reverse();heights.reverse();slopes.reverse()
    result=_SourceRule(curve,first,last,*heights,*slopes,
        phase=(float(stations[ids[0]]/stations[-1]),float(stations[ids[1]]/stations[-1])),
        corner_fragment=True)
    result.source_anchor=False
    return result


def _legacy_bank_side(context,source,index,first,last):
    """Only unambiguous original bank runs may omit a station ledger.

    Endpoint checks compare the two sides of this already declared source;
    they never search other roads, create a connection, or infer a source End.
    """
    from .intersections import Approach
    from .paths import Path
    from .routes import Route
    curve=source.layout.boundaries[index]
    if not np.isfinite((first,last)).all() or not 0.<=first<last<=1.+1e-9:
        raise ValueError('Legacy native bank needs a finite increasing interval')
    runs=[]
    for label in sorted((s for s in source.segments if s.index==index),key=lambda s:s.first):
        if label.kind!='curb':continue
        if runs and abs(runs[-1][1]-label.first)<=1e-12:
            runs[-1]=(runs[-1][0],label.last)
        else:runs.append((label.first,label.last))
    retained=[run for run in runs if run[0]-1e-10<=first<last<=run[1]+1e-10]
    if len(retained)!=1:
        raise SourceCurbConditionUnavailable('Legacy curb interval crosses a cap or an unknown native bank')
    if source.region.owner.kind=='cycle_strip':
        sid=source.region.owner.id;path=context[Key('path',sid)];profile=context[Key('profile',sid)]
        if type(path) is not Path:raise SourceCurbConditionUnavailable('Legacy strip needs one original native Path')
        low,high=0.,path.length;offset,width=path.section.offset,path.section.width
        point=path.point
    else:
        approaches=getattr(source.layout,'approaches',())
        if len(approaches)!=1 or type(approaches[0]) is not Approach:
            raise SourceCurbConditionUnavailable('Legacy port needs one original Approach')
        a=approaches[0];path,profile=a.path,a.profile
        if (type(path) not in (Path,Route) or hasattr(path,'mouth_piece')
                or profile!=context[Key('profile',a.source.id)].oriented(a.end)):
            raise SourceCurbConditionUnavailable('Legacy port needs an explicit original profile/route contract')
        low,high=0.,float(source.layout.outward_stations[0]);offset,width=a.offset,a.width;point=a.point
        if type(path) is Route and high>path.pieces[0].path.length:
            raise SourceCurbConditionUnavailable('Legacy port cannot cross an unrecorded native Route seam')
    endpoints=curve(np.asarray(retained[0]));matches=[]
    for side in (-1,1):
        expected=point(np.array((low,high)),offset+.5*width*side)
        if (np.max(np.linalg.norm(endpoints-expected,axis=1))<=1e-8
                or np.max(np.linalg.norm(endpoints-expected[::-1],axis=1))<=1e-8):matches.append(side)
    if len(matches)==1:return matches[0]
    # A full symmetric crown has identical outer-bank values and crossfall
    # magnitudes on both sides. No side identity is claimed in this fallback.
    if (not profile.half and profile.value(-1.)==profile.value(1.)
            and abs(profile.derivative(-1.))==abs(profile.derivative(1.))):return 1
    raise SourceCurbConditionUnavailable('Legacy half-profile bank needs explicit native side coordinates')


class _MappedRule:
    """Use a normalized carrier's explicit original span/vertex ledger.

    Only the already certified polyline parameter map changes. Jets use the
    current physical chord's one-sided normal, as ordinary CurbEdge does.
    There is no nearest-source or geometric ownership inference.
    """
    def __init__(self, resolver, curve, first, last, owner, index):
        self.curve, self.first, self.last = curve, float(first), float(last)
        cached = resolver.mapped.get((owner, index))
        if cached is None:
            original = np.asarray(curve.original_polyline)
            stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(original, axis=0), axis=1))]
            positions, lookup = [], {}
            for i, point in enumerate(original):
                key = tuple(point)
                lookup[key] = i if key not in lookup else None
            for point in curve.points:
                index_at = lookup.get(tuple(point))
                if index_at is None:
                    raise ValueError('Normalized source condition lost an explicit retained source vertex')
                positions.append(index_at)
            if (positions != sorted(positions) or positions[0] != 0
                    or positions[-1] != len(original)-1
                    or len(curve.source_spans) != len(original)-1):
                raise ValueError('Normalized source condition changed original point/span order')
            layout = resolver.sources[owner].layout
            before = after = None
            if layout.approaches and not getattr(curve, 'on_separator_rim', False):
                if index >= len(layout.approaches):
                    raise SourceCurbConditionUnavailable('Mapped exterior curb lost its two bounding mouths')
                a, b = layout.approaches[index], layout.approaches[(index+1) % len(layout.approaches)]
                before = _Anchor(a.profile.datum+a.profile.value(-1.),
                    a.physical_gradient(layout.outward_stations[index], -1.))
                after = _Anchor(b.profile.datum+b.profile.value(1.),
                    b.physical_gradient(layout.outward_stations[(index+1) % len(layout.approaches)], 1.))
            # A hole has no physical mouth, even when the same layout also
            # carries exterior approaches. Its chain may reuse its
            # explicit source-span ledger when that ledger supplies all real
            # anchors. A missing/mixed end remains an explicit failure.
            rules = resolver.chain_rules(curve.source_spans, original, False, before=before, after=after)
            cached = stations, stations[positions], rules
            resolver.mapped[owner, index] = cached
        self.old_stations, self.old_at, self.rules = cached
        middle = .5*(first+last)*curve.length
        edge = int(np.clip(np.searchsorted(curve.stations, middle, side='right')-1, 0, len(curve.points)-2))
        delta = curve.points[edge+1]-curve.points[edge]
        self.tangent = delta/np.linalg.norm(delta)
        self.arc_scale = ((self.old_at[edge+1]-self.old_at[edge])/
                          (curve.stations[edge+1]-curve.stations[edge]))

    def sample(self, t):
        shape = np.asarray(t).shape
        at = (self.first+np.asarray(t)*(self.last-self.first))*self.curve.length
        old = np.interp(at, self.curve.stations, self.old_at).reshape(-1)
        edge = np.clip(np.searchsorted(self.old_stations, old, side='right')-1, 0, len(self.old_stations)-2)
        # A declared original knot has two source jets. At a physical edge's
        # ends use the branch inside this edge, retaining the old interior
        # arithmetic. Inverting its published normalized parameter by another
        # multiply must not accidentally select the opposite original span.
        self.polynomial_breaks()
        values=np.asarray(t).reshape(-1);bounds=self._prediction_bounds
        for endpoint,inside in ((0.,bounds[1]),(1.,bounds[-2])):
            selected=values==endpoint
            if not np.any(selected):continue
            anchor=self._prediction_anchors[endpoint][1]
            interior=self._prediction_anchors[inside][1]
            # Adjacent floats can have a midpoint that rounds to the next
            # span. The lower bound of this declared piece owns its right jet.
            branch=np.clip(np.searchsorted(self.old_stations,min(anchor,interior),side='right')-1,
                           0,len(self.old_stations)-2)
            changed=selected & (edge!=branch)
            old[changed]=anchor;edge[changed]=branch
        heights, gradients = np.empty(len(old)), np.empty((len(old), 2))
        for i in np.unique(edge):
            selected = edge == i
            fraction = ((old[selected]-self.old_stations[i])/
                        (self.old_stations[i+1]-self.old_stations[i]))
            heights[selected], gradient = self.rules[int(i)].sample(fraction)
            tangent = ((self.curve.original_polyline[i+1]-self.curve.original_polyline[i])/
                       (self.old_stations[i+1]-self.old_stations[i]))
            along, inward = np.asarray(gradient)@tangent, np.asarray(gradient)@-right(tangent)
            gradients[selected] = (along[:, None]*self.arc_scale*self.tangent
                                   -inward[:, None]*right(self.tangent))
        return heights.reshape(shape), gradients.reshape(shape+(2,))

    def polynomial_breaks(self):
        """Declare the carrier's actual span map and native rule-piece seams.

        These are construction-table knots, not an adaptive integration table.
        Keep the old sampler unchanged when one existing polynomial suffices.
        """
        cached=getattr(self,'_prediction_bounds',None)
        if cached is not None:return cached[1:-1]
        low,high=self.first*self.curve.length,self.last*self.curve.length
        if not low<high:raise ValueError('Mapped source jet needs one positive finite carrier interval')
        old_low,old_high=map(float,np.interp((low,high),self.curve.stations,self.old_at))
        # Exact published table parameters identify their original construction
        # knots. No distance or near-equality test creates an extra alias.
        table=list(zip(self.curve.stations,self.old_at))
        table.extend((float(np.interp(original,self.old_at,self.curve.stations)),float(original))
                     for original in self.old_stations)
        for current,original in table:
            parameter=float(current/self.curve.length)
            if self.first==parameter:low,old_low=float(current),float(original)
            if self.last==parameter:high,old_high=float(current),float(original)
        anchors={0.:(low,old_low,None),1.:(high,old_high,None)}

        def add(current,original,child=None):
            if not low<current<high:return
            parameter=float((current-low)/(high-low))
            if not 0.<parameter<1.:return
            previous=anchors.get(parameter)
            if previous is not None and previous[:2]!=(current,original):
                raise ValueError('Mapped source rule has unrepresentable distinct construction knots')
            anchors[parameter]=(float(current),float(original),child)

        for current,original in zip(self.curve.stations,self.old_at):
            add(float(current),float(original))
        for original in self.old_stations:
            if old_low<original<old_high:
                add(float(np.interp(original,self.old_at,self.curve.stations)),float(original))
        begin=max(0,int(np.searchsorted(self.old_stations,old_low,side='right')-1))
        finish=min(len(self.rules),int(np.searchsorted(self.old_stations,old_high,side='left'))+1)
        for index in range(begin,finish):
            rule=self.rules[index];declared=getattr(rule,'polynomial_breaks',None)
            breaks=() if declared is None else declared()
            if (not isinstance(breaks,tuple)
                    or any(not np.isfinite(value) or not 0.<value<1. for value in breaks)
                    or any(a>=b for a,b in zip(breaks,breaks[1:]))):
                raise ValueError('Mapped native source rule has invalid declared construction pieces')
            a,b=self.old_stations[index:index+2]
            for parameter in breaks:
                original=float(a+parameter*(b-a))
                if old_low<original<old_high:
                    add(float(np.interp(original,self.old_at,self.curve.stations)),original,(index,parameter))
        bounds=tuple(sorted(anchors))
        if len(bounds)>20001:raise ValueError('Mapped source jets exceed the existing prediction branch budget')
        self._prediction_bounds,self._prediction_anchors=bounds,anchors
        return bounds[1:-1]

    def sample_piece(self,t,first,last):
        """Evaluate one declared map piece with its true one-sided native jets."""
        self.polynomial_breaks()
        bounds=self._prediction_bounds
        if not any(a==first and b==last for a,b in zip(bounds,bounds[1:])):
            raise ValueError('Mapped source jet leaves its declared consecutive construction piece')
        start,end=self._prediction_anchors[first],self._prediction_anchors[last]
        current,original=start[:2]
        carrier=int(np.clip(np.searchsorted(self.curve.stations,current,side='right')-1,
                            0,len(self.curve.points)-2))
        index=int(np.clip(np.searchsorted(self.old_stations,original,side='right')-1,
                          0,len(self.old_stations)-2))
        a,b=self.old_stations[index:index+2]
        parameters=[]
        for anchor in (start,end):
            native=anchor[2]
            parameters.append(native[1] if native is not None and native[0]==index
                              else float((anchor[1]-a)/(b-a)))
        begin,finish=parameters
        if not 0.<=begin<finish<=1.:
            raise ValueError('Mapped source jet changed its finite original rule domain')
        values=np.asarray(t);shape=values.shape;rule=self.rules[index]
        piece=getattr(rule,'sample_piece',None)
        local=begin+values*(finish-begin)
        if piece is not None:
            declared=getattr(rule,'polynomial_breaks',lambda:())()
            native=(0.,*declared,1.)
            branch=int(np.clip(np.searchsorted(native,begin,side='right')-1,0,len(native)-2))
            low,high=native[branch:branch+2]
            if begin<low or finish>high:
                raise ValueError('Mapped source jet crosses an undeclared native rule piece')
            height,gradient=piece((local-low)/(high-low),low,high)
        else:height,gradient=rule.sample(local)
        gradient=np.asarray(gradient).reshape(-1,2)
        old_tangent=(self.curve.original_polyline[index+1]-self.curve.original_polyline[index])/(b-a)
        current_length=self.curve.stations[carrier+1]-self.curve.stations[carrier]
        tangent=(self.curve.points[carrier+1]-self.curve.points[carrier])/current_length
        scale=(self.old_at[carrier+1]-self.old_at[carrier])/current_length
        along,inward=gradient@old_tangent,gradient@-right(old_tangent)
        mapped=along[:,None]*scale*tangent-inward[:,None]*right(tangent)
        return np.asarray(height).reshape(shape),mapped.reshape(shape+(2,))


class SourceCurbConditions:
    def __init__(self, context, group, *, maximum_slope=None):
        if maximum_slope is not None and (isinstance(maximum_slope, (bool, np.bool_))
                or not isinstance(maximum_slope, (int, float, np.integer, np.floating))
                or not math.isfinite(maximum_slope) or maximum_slope <= 0.):
            raise ValueError('Source curb transition needs its positive finite physical slope limit')
        self.context, self.group = context, group
        self.sources = {source.region.owner: source for source in group.footprints}
        self.cache, self.mapped, self.native_knots = {}, {}, {}
        self.one_sided_profiles = {}
        self.native_ends = {}
        self.maximum_slope, self.transitions = maximum_slope, []

    def _endpoint(self, rule, end):
        # Each native piece owns its one-sided gradient even where two curbs
        # meet at the same XY. The cache lives only for this boundary assembly.
        if type(rule) is not _NativeCoordinateRule:
            return rule.sample(end)
        key = rule, end
        if key not in self.native_ends:
            self.native_ends[key] = rule.sample(end)
        return self.native_ends[key]

    def derived_rule(self, curve, first, last, depth=0):
        """Follow only declared finite piece maps; never fit another source."""
        from .compound_curbs import CompositeCurb, OffsetCurb, CircularCurb, SamplingBudgetCurb
        from .bundle_merged_rim import RetainedCurbPiece
        from .compound_dimensional_corners import _CertifiedContinuation
        from .free_pads import EnvelopeBoundary
        if depth >= 16:
            raise ValueError('Source curb condition exceeds its finite wrapper depth')
        if type(curve) is EnvelopeBoundary:
            from .free_pad_curb_conditions import validate
            domain=validate(self.context,self.group,self.sources,curve)
            rule=_SourceRule(curve,first,last,*domain.heights,*domain.inward_slopes,corner_fragment=True)
            rule.source_anchor=False
            return rule
        from .retained_arc_fillets import OriginalAnalyticalCircleCurve
        if isinstance(curve,OriginalAnalyticalCircleCurve):
            corner=curve.corner
            a,b=corner.first,corner.second
            z0,z1=a.profile.datum+float(a.profile.value(-1.)),b.profile.datum+float(b.profile.value(1.))
            k0,k1=_corner_slopes(corner)
            if curve.circle_interval is not None:
                from .source_offset_rules import SourceCurvePieces,SourceOffsetRule
                declared=curve.circle_interval;declared.validate(curve.original_prepared_curve,corner)
                bounds=sorted({first,last,*(value for value in declared.prepared_parameters if first<value<last),
                    *(value for lead in declared.auxiliary_leads for value in lead.prepared_parameters if first<value<last)})
                pieces=[];low,high=declared.prepared_bounds
                for u,v in zip(bounds,bounds[1:]):
                    if low<=.5*(u+v)<=high:
                        start,end=map(float,declared.circle_parameter(np.clip((u,v),low,high)))
                        rule=_SourceRule(CircularCurb(corner),start,end,z0,z1,k0,k1,corner_fragment=True)
                    else:
                        leads=[lead for lead in declared.auxiliary_leads
                               if lead.prepared_bounds[0]-1e-12<=u<v<=lead.prepared_bounds[1]+1e-12]
                        if len(leads)!=1:
                            raise SourceCurbConditionUnavailable('Original circular source lost its explicit auxiliary lead domain')
                        lead=leads[0];begin,end=lead.prepared_bounds
                        actual=lead.original_approach
                        mapped=SimpleNamespace(approach=actual,parameters=np.asarray(lead.outward_stations),
                            cumulative=(np.asarray(lead.prepared_parameters)-begin)*curve.length,
                            offset=actual.offset+.5*actual.width*lead.bank_side,length=(end-begin)*curve.length)
                        rule=SourceOffsetRule(self.context,mapped,max(0.,(u-begin)/(end-begin)),
                                              min(1.,(v-begin)/(end-begin)))
                        # Producer-declared negative auxiliary construction is
                        # generated junction geometry, never a native body bank.
                        rule.corner_fragment=True;rule.source_anchor=False
                    pieces.append((u,v,rule))
                return SourceCurvePieces(first,last,pieces)
            return _SourceRule(curve,first,last,z0,z1,k0,k1,corner_fragment=True)
        if isinstance(curve, SamplingBudgetCurb):
            return self.derived_rule(curve.curve, first, last, depth+1)
        if isinstance(curve, RetainedCurbPiece):
            return self.derived_rule(curve.original, curve.start+first*(curve.stop-curve.start),
                                     curve.start+last*(curve.stop-curve.start), depth+1)
        if isinstance(curve, CompositeCurb):
            low, high = first*curve.length, last*curve.length
            index = int(np.clip(np.searchsorted(curve.lengths, .5*(low+high), side='right')-1,
                                0, len(curve.pieces)-1))
            begin, end = curve.lengths[index:index+2]
            tolerance = 1e-9*max(1., curve.length)
            if low < begin-tolerance or high > end+tolerance:
                from .source_offset_rules import SourceCurvePieces
                pieces=[]
                for i,piece in enumerate(curve.pieces):
                    a,b=max(low,curve.lengths[i]),min(high,curve.lengths[i+1])
                    if not a<b:continue
                    rule=self.derived_rule(piece,(a-curve.lengths[i])/piece.length,
                        (b-curve.lengths[i])/piece.length,depth+1)
                    if rule is None:
                        raise SourceCurbConditionUnavailable('Composite source curb lost an original piece rule')
                    pieces.append((max(first,a/curve.length),min(last,b/curve.length),rule))
                return SourceCurvePieces(first,last,pieces)
            piece = curve.pieces[index]
            return self.derived_rule(piece, max(0., (low-begin)/piece.length),
                                     min(1., (high-begin)/piece.length), depth+1)
        from .local_width_continuation import LocalBank, LocalTransition, TransitionHeight
        if isinstance(curve, LocalBank):
            from .source_offset_rules import SourceOffsetRule
            return SourceOffsetRule(self.context, curve, first, last)
        if isinstance(curve, LocalTransition):
            return TransitionHeight(curve.curve, first, last)
        if isinstance(curve, OffsetCurb):
            return _OffsetRule(curve,first,last)
        if isinstance(curve, (_CertifiedContinuation, CircularCurb)):
            actual = curve.curve if isinstance(curve, _CertifiedContinuation) else curve
            corner = actual.corner
            a, b = corner.first, corner.second
            z0, z1 = a.profile.datum+float(a.profile.value(-1.)), b.profile.datum+float(b.profile.value(1.))
            if isinstance(curve, _CertifiedContinuation) and abs(z0-z1) > 1e-10:
                raise SourceCurbConditionUnavailable('Mixed-height polynomial core needs its explicit physical speed/jet map')
            k0, k1 = _corner_slopes(corner)
            return _SourceRule(curve, first, last, z0, z1, k0, k1, corner_fragment=True)
        from .intersections import BoundaryCurve
        if isinstance(curve, BoundaryCurve):
            from .source_offset_rules import whole_original_curve
            complete=whole_original_curve(self.context,self.group,self.sources,curve)
            corner = curve.corner
            a, b = corner.first, corner.second
            z0, z1 = a.profile.datum+float(a.profile.value(-1.)), b.profile.datum+float(b.profile.value(1.))
            k0, k1 = _corner_slopes(corner)
            phase = ((float(curve.lengths[1]/curve.length), float(curve.lengths[2]/curve.length))
                     if corner.center is not None else (0., 1.))
            if corner.center is not None:
                low, high = first*curve.length, last*curve.length
                for index, begin, end in ((0, 0., curve.lengths[1]),
                                           (1, curve.lengths[2], curve.lengths[3])):
                    if end > begin and low >= begin-1e-10 and high <= end+1e-10:
                        approach, parameters, cumulative, offset = curve.leads[index]
                        lead = SimpleNamespace(approach=approach, parameters=parameters,
                            cumulative=cumulative, offset=offset, length=float(end-begin))
                        from .source_offset_rules import SourceOffsetRule
                        return SourceOffsetRule(self.context,lead, max(0., (low-begin)/lead.length),
                                                min(1., (high-begin)/lead.length))
                if (any(first+1e-10 < contact < last-1e-10 for contact in phase)
                        and (complete or any(hasattr(approach,'field') or _needs_exact_station_ledger(approach.path)
                                             for approach in (a,b)))):
                    from .source_offset_rules import SourceCurvePieces,SourceOffsetRule
                    bounds=(first,*(value for value in phase if first<value<last),last)
                    pieces=[]
                    for u,v in zip(bounds,bounds[1:]):
                        middle=.5*(u+v)
                        if phase[0]<middle<phase[1]:
                            rule=_SourceRule(curve,u,v,z0,z1,k0,k1,phase,corner_fragment=not complete)
                        else:
                            index=0 if middle<phase[0] else 1
                            begin,end=(0.,curve.lengths[1]) if index==0 else tuple(curve.lengths[2:4])
                            approach,parameters,cumulative,offset=curve.leads[index]
                            lead=SimpleNamespace(approach=approach,parameters=parameters,
                                cumulative=cumulative,offset=offset,length=float(end-begin))
                            rule=SourceOffsetRule(self.context,lead,
                                max(0.,(u*curve.length-begin)/lead.length),
                                min(1.,(v*curve.length-begin)/lead.length))
                        pieces.append((u,v,rule))
                    return SourceCurvePieces(first,last,pieces)
            return _SourceRule(curve, first, last, z0, z1, k0, k1, phase, corner_fragment=not complete)
        return None

    def source_rule(self, owner, index, first, last):
        key = owner, index, float(first), float(last)
        if key in self.cache:
            return self.cache[key]
        source = self.sources[owner]
        curve = source.layout.boundaries[index]
        from .retained_arc_fillets import RetainedArcFillet,OriginalAnalyticalCircleCurve
        from .source_bank_reconciliation import PreparedBankFillet
        from .source_owned_connector_fillets import PreparedOwnedConnectorFillet
        from .source_parent_chain_fillets import PreparedParentChainFillet
        bank_fillet=getattr(source.layout,'source_bank_fillet',None)
        owned_connector=getattr(source.layout,'source_owned_connector_fillet',None)
        from .bank_circle_contacts import BankCircleFillet
        from .arc_lead_joins import ArcLeadJoin
        from .exposed_bank_joins import BankJoin
        from .bundle_crest_roots import CrestEnvelope
        crest=getattr(source.layout,'generated_junction_area',None)
        if hasattr(curve, 'free_cap_banks'):
            from .finite_free_caps import OwnedFiniteFreeCap
            claims = [c for c in source.layout.finite_free_pads if isinstance(c, OwnedFiniteFreeCap)]
            if index != 0 or len(claims) != 1 or len(curve.free_cap_banks) != len(curve.stations)-1:
                raise SourceCurbConditionUnavailable('Free-cap curb lost its original two-bank inventory')
            claim, = claims
            at = int(np.searchsorted(curve.stations, .5*(first+last)*curve.length, side='right')-1)
            start, end = curve.stations[at:at+2]/curve.length
            if first < start-1e-9 or last > end+1e-9:
                raise SourceCurbConditionUnavailable('Free-cap source interval crosses an unlabelled bank change')
            bank = curve.free_cap_banks[at]
            if bank is None:
                # The generated zero curb acquires its inward crossfall from
                # the neighboring surviving original banks or mouth profiles.
                datum = claim.pad.profile.datum
                result = _SourceRule(curve, first, last, datum, datum, 0., 0., corner_fragment=True)
                result.source_anchor = False
            else:
                role, side, sa, sb = bank
                arm = claim.approaches[role]
                if arm is not source.layout.approaches[role]:
                    raise SourceCurbConditionUnavailable('Free cap changed its original bank approach')
                height = arm.profile.datum + float(arm.profile.value(side))
                s0, s1 = np.interp((first, last), (start, end), (sa, sb))
                k0, k1 = (float(np.linalg.norm(arm.physical_gradient(float(s), side))) for s in (s0, s1))
                result = _SourceRule(curve, first, last, height, height, k0, k1, phase=(first,last))
        elif isinstance(crest,CrestEnvelope):
            if curve is not crest.curve or source.region.owner!=crest.source.region.owner:
                raise SourceCurbConditionUnavailable('Bundle crest condition lost its original curb')
            crest.validate(self.context,self.group)
            result=_SourceRule(curve,first,last,crest.reference.height,crest.reference.height,0.,0.,corner_fragment=True)
        elif isinstance(curve,(ArcLeadJoin,BankJoin)):
            if index!=0 or source.region.owner!=curve.owner or crest is not curve:
                raise SourceCurbConditionUnavailable('Arc/lead join lost its original boundary identity')
            result=curve.height_rule(self,first,last)
        elif isinstance(curve,BankCircleFillet):
            from .bank_circle_contacts import height_rule
            result=height_rule(self,source,index,curve,first,last)
        elif isinstance(curve,RetainedArcFillet):
            if index!=0 or getattr(source.layout,'generated_junction_area',None) is not curve:
                raise SourceCurbConditionUnavailable('Retained arc lost its typed boundary/source ancestry')
            heights=[];slopes=[]
            for parent in curve.parent_arcs:
                parent_source=self.sources[parent.owner]
                bound=parent_source.layout.boundaries[parent.curve_index]
                if (not isinstance(bound,OriginalAnalyticalCircleCurve)
                        or bound.corner is not parent.original_corner
                        or bound.original_prepared_curve is not parent.original_curve):
                    raise SourceCurbConditionUnavailable('Retained arc changed its exact original circular parents')
                from .displaced_corner_authority import original_parent_corner
                node=next((item for item in self.context.network.intersections
                           if item.key==getattr(parent_source.layout,'original_junction',None)),None)
                if node is None or node.key not in self.group.members or not original_parent_corner(self,node,parent.original_corner):
                    raise SourceCurbConditionUnavailable('Retained arc lost its actual parent source/profile authority')
                low,high=sorted(parent.exposed_circle_parameters)
                if not low<=parent.analytical_contact_parameter<=high:
                    raise SourceCurbConditionUnavailable('Retained arc height contact leaves its finite parent')
                c=parent.original_corner
                heights.extend((c.first.profile.datum+float(c.first.profile.value(-1.)),
                                c.second.profile.datum+float(c.second.profile.value(1.))))
                k0,k1=_corner_slopes(c);t=parent.analytical_contact_parameter
                weight=t**3*(10.+t*(-15.+6.*t));slopes.append(k0+weight*(k1-k0))
            # Constant true parent levels cover zero road blocks and genuine
            # inner peaks. Mixed generated crossing values are never hard.
            z0,z1=(heights[0],heights[0]) if max(heights)-min(heights)<=1e-10 else (heights[0],heights[-1])
            result=_SourceRule(curve,first,last,z0,z1,*slopes,corner_fragment=True)
            result.source_anchor=False
        elif isinstance(curve,PreparedParentChainFillet):
            from .source_parent_chain_conditions import build_rule
            result=build_rule(self,source,index,curve,first,last)
        elif isinstance(owned_connector,PreparedOwnedConnectorFillet):
            result=_owned_connector_rule(self,source,curve,first,last,owned_connector)
        elif isinstance(bank_fillet,PreparedBankFillet):
            if getattr(source.layout,'generated_junction_area',None) is not bank_fillet:
                raise SourceCurbConditionUnavailable('Bank fillet is not a typed geometry-only source area')
            c=bank_fillet.fitted.corner;heights=[];slopes=[];contacts=[]
            for bank,arm,at,cut in zip((bank_fillet.proposal.first,bank_fillet.proposal.last),
                    (c.first,c.second),c.stations,bank_fillet.fitted.intersection):
                sid=bank.path.section.id
                native=float(at if arm.end==0 else bank.path.length-at)
                native_cut=float(cut if arm.end==0 else bank.path.length-cut)
                if (bank.path is not self.context[Key('path',sid)] or arm.path is not bank.path
                        or not any(a<=min(native,native_cut)<=max(native,native_cut)<=b for a,b in bank.intervals)):
                    raise SourceCurbConditionUnavailable('Generated bank arc lost its finite original native support')
                profile=self.context[Key('profile',sid)]
                if (bank.side not in (-1,1) or arm.source is not bank.path.section
                        or arm.width!=bank.path.section.width
                        or arm.offset!=bank.path.section.offset*(1 if arm.end==0 else -1)
                        or arm.profile!=profile.oriented(arm.end)):
                    raise SourceCurbConditionUnavailable('Bank fillet changed its true native profile/frame')
                heights.append(profile.datum+float(profile.value(bank.side)))
                q=bank.side*(1 if arm.end==0 else -1)
                slopes.append(float(np.linalg.norm(arm.physical_gradient(float(at),q))))
                contacts.append(arm.point(float(at),arm.offset+.5*arm.width*q))
            if not hasattr(curve,'points'):
                raise SourceCurbConditionUnavailable('Generated bank arc needs explicit prepared contact knots')
            points=np.asarray(curve.points);stations=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
            ids=[int(np.argmin(np.linalg.norm(points-contact,axis=1))) for contact in contacts]
            if max(np.linalg.norm(points[i]-point) for i,point in zip(ids,contacts))>1e-7 or ids[0]==ids[1]:
                raise SourceCurbConditionUnavailable('Generated bank arc lost its original native contact points')
            if ids[0]>ids[1]:ids.reverse();heights.reverse();slopes.reverse()
            result=_SourceRule(curve,first,last,*heights,*slopes,
                phase=(float(stations[ids[0]]/stations[-1]),float(stations[ids[1]]/stations[-1])),corner_fragment=True)
            result.source_anchor=False
        elif getattr(source.layout,'one_sided_source_preparation',False) and hasattr(curve,'one_sided_edges'):
            from .one_sided_curb_heights import exterior_profile, separator_level
            profile_key = owner,index
            if profile_key not in self.one_sided_profiles:
                pairs = getattr(source.layout, 'one_sided_separator_pairs', ())
                pair = next((indices for indices in pairs if index in indices), None)
                if pair is None:
                    self.one_sided_profiles[profile_key] = exterior_profile(self.context,source.layout,index)
                else:
                    level = separator_level(self.context, source.layout, pair)
                    for i in pair:
                        self.one_sided_profiles[owner,i] = level, 0., 0., (0.,1.)
            z,k0,k1,phase = self.one_sided_profiles[profile_key]
            result = _SourceRule(curve,first,last,z,z,k0,k1,phase)
        elif hasattr(curve, 'source_spans') and hasattr(curve, 'original_polyline'):
            result = _MappedRule(self, curve, first, last, owner, index)
        elif owner.kind in ('cycle_strip','cycle_port') and getattr(curve,'source_coordinates',None) is not None:
            from .owned_prefix_source_curbs import validate_curve
            validate_curve(self,source,index,curve,first,last)
            result = _NativeCoordinateRule(self.context,curve,first,last,source_knots=self.native_knots)
        elif owner.kind == 'cycle_strip':
            if index not in (0, 1):
                raise ValueError('Original source strip lost its two explicit curb sides')
            profile = self.context[Key('profile', owner.id)]
            path = self.context[Key('path', owner.id)]
            if path.section.id != owner.id:
                raise ValueError('Original source strip changed its Path identity')
            q = _legacy_bank_side(self.context,source,index,first,last)
            height = profile.datum+float(profile.value(q))
            slope = abs(float(profile.derivative(q)))*2./path.section.width
            if slope and _needs_exact_station_ledger(path):
                raise SourceCurbConditionUnavailable('Prepared Exact strip jets need their original canonical station ledger')
            result = _SourceRule(curve, first, last, height, height, slope, slope)
        elif owner.kind == 'cycle_port':
            if index not in (0, 1) or len(source.layout.approaches) != 1:
                raise ValueError('Original source port lost its two explicit curb sides')
            approach = source.layout.approaches[0]
            q = _legacy_bank_side(self.context,source,index,first,last)
            height = approach.profile.datum+float(approach.profile.value(q))
            slope = abs(float(approach.profile.derivative(q)))*2./approach.width
            if slope and _needs_exact_station_ledger(approach.path):
                raise SourceCurbConditionUnavailable('Prepared Exact port jets need their original canonical station ledger')
            result = _SourceRule(curve, first, last, height, height, slope, slope)
        elif getattr(source.layout, 'terminal_merged_port', None):
            if len(source.layout.approaches) != 1:
                raise ValueError('Merged source curb lost its unique prepared physical entrance')
            approach = source.layout.approaches[0]
            if not hasattr(approach, 'original_approaches') or not hasattr(approach, 'source_ranges'):
                raise ValueError('Merged source curb lost its explicit original source aliases/ranges')
            from .terminal_curb_ancestry import validate
            ancestry = getattr(source.layout, 'terminal_curb_ancestry', {}).get(index)
            if ancestry is None:
                raise SourceCurbConditionUnavailable('Merged source curb is missing its original bank ancestry')
            validate(approach, ancestry, curve.points)
            from .finite_terminal_entrances import LinearCrown
            rule = _LinearTerminalRule if isinstance(approach.field,LinearCrown) and len(curve.points)==2 else _FieldRule
            result = rule(curve, first, last, approach.field,
                                original=(self.context, approach, ancestry))
        elif hasattr(curve, 'corner') and hasattr(curve.corner, 'first'):
            result = self.derived_rule(curve, first, last)
            if result is None:
                raise SourceCurbConditionUnavailable('Original curb needs an explicit finite lead/core domain')
        elif len(source.layout.corners) == 1:
            corner = source.layout.corners[0]
            a, b = corner.first, corner.second
            z0, z1 = a.profile.datum+float(a.profile.value(-1.)), b.profile.datum+float(b.profile.value(1.))
            k0, k1 = _corner_slopes(corner)
            if not hasattr(curve, 'points'):
                raise SourceCurbConditionUnavailable('Source corner needs its prepared tangent-contact knots')
            points = np.asarray(curve.points)
            stations = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
            if corner.center is None:
                aa = a.point(corner.stations[0], a.offset-.5*a.width)
                bb = b.point(corner.stations[1], b.offset+.5*b.width)
            else:
                aa = corner.center+corner.radius*np.array((np.cos(corner.start_angle), np.sin(corner.start_angle)))
                bb = corner.center+corner.radius*np.array((np.cos(corner.start_angle+corner.sweep), np.sin(corner.start_angle+corner.sweep)))
            ia, ib = (int(np.argmin(np.linalg.norm(points-point, axis=1))) for point in (aa, bb))
            if max(np.linalg.norm(points[ia]-aa), np.linalg.norm(points[ib]-bb)) > 1e-7 or ia == ib:
                raise SourceCurbConditionUnavailable('Source corner lost its distinct original tangent contacts')
            if ia > ib:
                ia, ib, z0, z1, k0, k1 = ib, ia, z1, z0, k1, k0
            phase = float(stations[ia]/stations[-1]), float(stations[ib]/stations[-1])
            result = _SourceRule(curve, first, last, z0, z1, k0, k1, phase, corner_fragment=True)
        else:
            raise SourceCurbConditionUnavailable('No explicit source curb condition for '+owner.text())
        self.cache[key] = result
        return result

    def chain_rules(self, spans, points, closed, *, before=None, after=None, original=None):
        if self.maximum_slope is None:
            return self._chain_rules(spans, points, closed, before=before, after=after, original=original)
        baseline, failure = None, None
        try:
            baseline = self._chain_rules(spans, points, closed, before=before, after=after, original=original)
        except (SourceAnchorConflict, SourceJetDirectionError) as error:
            failure = error
        from .junction_bank_transitions import reconcile
        rules, reports = reconcile(self, spans, points, closed, before=before, after=after,
                                    original=original, baseline=baseline, original_error=failure)
        self.transitions.extend(reports)
        return rules

    def _chain_rules(self, spans, points, closed, *, before=None, after=None, original=None, soften=()):
        count, rules, source_rules = len(spans), [], []
        if not 1 <= count <= 20000:
            raise ValueError('Source curb conditions exceed their finite boundary budget')
        for i, span in enumerate(spans):
            label = rim_label(span, self.sources)
            source = self.sources[span.owner]
            if label.kind == 'mouth':
                if original is None:
                    raise ValueError('Source chain needs its original incoming profile condition')
                rule = original[i]
            elif label.kind == 'curb':
                rule = self.source_rule(span.owner, label.index, label.first, label.last)
                if getattr(rule, 'corner_fragment', False):
                    if not rule.source_anchor:
                        rule = None
                if i in soften:
                    rule = None
            else:
                raise ValueError('Source boundary conditions cannot hide an internal preparation seam')
            source_rules.append(rule if rule is not None else self.source_rule(
                span.owner, label.index, label.first, label.last))
            rules.append(rule)
        for i, following in enumerate(rules):
            if following is None or (i == 0 and not closed):
                continue
            previous = rules[(i-1) % count]
            if previous is None:
                continue
            za, ga = self._endpoint(previous, 1.)
            zb, gb = self._endpoint(following, 0.)
            if (abs(float(za)-float(zb)) > 1e-7
                    or abs(float(np.linalg.norm(ga))-float(np.linalg.norm(gb))) > 1e-7):
                if (self.maximum_slope is not None and abs(float(za)-float(zb))<=1e-10
                        and max(np.linalg.norm(ga),np.linalg.norm(gb))<=self.maximum_slope):
                    from .junction_bank_transitions import terminal_contact_arc
                    arc=terminal_contact_arc(self,spans,rules,i,closed)
                    if arc:
                        for index in arc:
                            rules[index]=None
                        continue
                raise SourceAnchorConflict('Adjacent original source anchors disagree in physical height or crossfall').at_join(spans,points,i)
        if closed and not any(rule is not None for rule in rules):
            # A closed chain containing only constant original corner rules
            # has an unambiguous source level, even when all tangent contacts
            # are clipped. This covers ordinary zero-curb blocks and genuine
            # inner plateaus without inventing a hidden mixed-height contact.
            height = getattr(source_rules[0], 'z0', None)
            if height is not None and all(isinstance(rule, _SourceRule)
                    and abs(rule.z0-height) <= 1e-10 and abs(rule.z1-height) <= 1e-10
                    for rule in source_rules):
                return source_rules
            raise SourceCurbConditionUnavailable('Closed source curb needs an exposed original source anchor')
        starts = [i for i, rule in enumerate(rules) if rule is None
                  and (i == 0 and not closed or rules[(i-1) % count] is not None)]
        for start in starts:
            ids, at = [], start
            while rules[at] is None:
                ids.append(at)
                at = (at+1) % count
                if not closed and at == 0:
                    break
                if len(ids) > count:
                    raise ValueError('Unbounded source corner chain')
            previous = rules[(start-1) % count] if start > 0 or closed else before
            following = rules[at] if at > 0 or closed else after
            if previous is None or following is None:
                raise SourceCurbConditionUnavailable('Source corner fragments lost their surviving boundary anchors')
            za, ga = self._endpoint(previous, 1.)
            zb, gb = self._endpoint(following, 0.)
            vectors = np.asarray([points[(i+1) % len(points)]-points[i] for i in ids])
            lengths = np.linalg.norm(vectors, axis=1)
            units = vectors/lengths[:, None]
            alignments=float(np.asarray(ga)@-right(units[0])),float(np.asarray(gb)@-right(units[-1]))
            if min(alignments)<-1e-6:
                error=SourceJetDirectionError('Exposed original source jet faces outside its retained corner chain')
                at=ids[0] if alignments[0]<=alignments[1] else (ids[-1]+1)%count
                if closed or at!=0:
                    error.at_join(spans,points,at)
                raise error
            # Preserve the original crossfall magnitude. A sampled chord's
            # one-sided normal is used consistently on every retained edge.
            ka, kb = float(np.linalg.norm(ga)), float(np.linalg.norm(gb))
            stations = np.r_[0., np.cumsum(lengths)]
            for j, i in enumerate(ids):
                rules[i] = CurbEdge(float(za), float(zb), ka, kb, float(stations[-1]),
                    float(stations[j]), float(stations[j+1]), tuple(units[j]))
        return rules

    def apply(self, loop, data):
        rules = self.chain_rules(loop.spans, loop.xy, True, original=data.rules)
        values = [self._endpoint(rule, 0.) for rule in rules]
        heights, gradients = np.array([row[0] for row in values]), np.array([row[1] for row in values])
        for port in data.ports:
            heights[list(port.vertices)] = data.heights[list(port.vertices)]
            gradients[list(port.vertices)] = data.gradients[list(port.vertices)]
        jump = max(abs(float(self._endpoint(rules[i-1], 1.)[0])-float(value[0]))
                   for i, value in enumerate(values))
        if jump > 1e-7:
            raise ValueError('Source curb heights do not join continuously: '+str(jump))
        if (not np.isfinite(gradients).all() or heights.min() < data.datum-1e-8
                or heights.max() > data.ceiling+1e-8):
            raise ValueError('Source curb conditions exceed the original physical envelope')
        heights.flags.writeable = gradients.flags.writeable = False
        result = replace(data, heights=heights, gradients=gradients, rules=tuple(rules))
        if self.transitions:
            from .junction_bank_transitions import TransitionBoundaryData
            return TransitionBoundaryData(result.xy, result.heights, result.gradients, result.rules,
                result.ports, result.datum, result.ceiling, result.hole, tuple(self.transitions))
        return result


def apply(context, group, loop, data, *, maximum_slope=None):
    """Explicit opt-in; callers retain the unchanged default height policy."""
    return SourceCurbConditions(context, group, maximum_slope=maximum_slope).apply(loop, data)
