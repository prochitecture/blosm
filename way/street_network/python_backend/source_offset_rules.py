"""Faithful jets on declared source lead maps and finite curve pieces.

Only prediction pieces are introduced. Physical curves, chords, mouth planes,
source profiles and original native station domains remain caller-owned data.
"""
import numpy as np
from .model import End,Intersection,Key
from .routes import Route
from .paths import Path
from .source_identity import matches_current_section


def whole_original_curve(context,group,sources,curve):
    """Actual unmodified original layout plus complete exposed source range."""
    from .compounds import rim_label
    from .compound_curbs import SamplingBudgetCurb
    from .ownership import footprint_key
    if (len(group.members)!=1 or group.claims or group.physical_layout is not None
            or len(group.footprints)!=1 or len(group.loops)!=1 or group.loops[0].hole):
        return False
    owner=group.members[0];source=sources.get(owner)
    node=context.network.records().get(owner)
    if (source is None or type(node) is not Intersection
            or source.layout is not context.products.get(footprint_key(owner))
            or source.region.boundary_only):return False
    layout=source.layout
    if ({End(a.source.id,a.end) for a in layout.approaches}!=set(node.ends)
            or len(layout.approaches)!=len(node.ends)):return False
    for a in layout.approaches:
        if (not matches_current_section(context,a.source)
                or hasattr(a,'field') or a.profile!=context[Key('profile',a.source.id)].oriented(a.end)
                or a.width!=a.source.width
                or a.offset!=a.source.offset*(1. if a.end==0 else -1.)):return False
        pieces=a.path.pieces if type(a.path) is Route else ()
        if type(a.path) is not Route or any(piece.path is not context[Key('path',piece.section.id)] for piece in pieces):
            return False
    for corner in layout.corners:
        end=End(corner.first.source.id,corner.first.end)
        if end not in node.ends:return False
        radius,mode=node.corner_radius(end)
        if (corner.radius_mode!=mode or corner.requested_radius!=radius
                or mode=='STRICT' and corner.radius!=radius):return False
    matched=[]
    for index,boundary in enumerate(layout.boundaries):
        while isinstance(boundary,SamplingBudgetCurb):boundary=boundary.curve
        if boundary is curve:matched.append(index)
    if len(matched)!=1:return False
    intervals=[]
    for span in group.loops[0].spans:
        if span.owner!=owner:return False
        label=rim_label(span,sources)
        if label.kind=='curb' and label.index==matched[0]:
            intervals.append((float(label.first),float(label.last)))
    intervals.sort()
    return bool(intervals and intervals[0][0]==0. and intervals[-1][1]==1.
        and all(a<b for a,b in intervals)
        and all(a[1]==b[0] for a,b in zip(intervals,intervals[1:])))


class SourceCurvePieces:
    def __init__(self, first, last, pieces):
        self.first, self.last = float(first), float(last)
        self.records = []
        for begin, end, rule in pieces:
            declared = getattr(rule, 'polynomial_breaks', None)
            bounds = (0., *(declared() if declared is not None else ()), 1.)
            for low, high in zip(bounds, bounds[1:]):
                a = (begin+low*(end-begin)-first)/(last-first)
                b = (begin+high*(end-begin)-first)/(last-first)
                self.records.append((a, b, rule, low, high))
        if (not self.records or len(self.records)>20000
                or self.records[0][0] != 0. or self.records[-1][1] != 1.
                or any(a[1] != b[0] for a,b in zip(self.records,self.records[1:]))):
            raise ValueError('Source curb pieces must retain one complete finite parameter interval')
        self._breaks = tuple(row[1] for row in self.records[:-1])
        self.ends = np.array([row[1] for row in self.records])
        self.corner_fragment = all(getattr(rule,'corner_fragment',False) for _,_,rule in pieces)
        self.source_anchor = all(getattr(rule,'source_anchor',True) for _,_,rule in pieces)
        self.z0 = float(np.asarray(self.sample(0.)[0]))
        self.z1 = float(np.asarray(self.sample(1.)[0]))

    def polynomial_breaks(self):
        return self._breaks

    def sample_piece(self,t,first,last):
        selected=[record for record in self.records if record[:2]==(first,last)]
        if len(selected)!=1:
            raise ValueError('Source jet piece must retain its declared original boundaries')
        _,_,rule,low,high=selected[0]
        sampler=getattr(rule,'sample_piece',None)
        return (sampler(t,low,high) if sampler is not None
                else rule.sample(low+np.asarray(t)*(high-low)))

    def sample(self,t):
        at=np.asarray(t);flat=at.reshape(-1)
        indices=np.clip(np.searchsorted(self.ends,flat,side='right'),0,len(self.records)-1)
        height=np.empty(len(flat));gradient=np.empty((len(flat),2))
        for index in np.unique(indices):
            selected=indices==index;first,last,*_=self.records[int(index)]
            local=(flat[selected]-first)/(last-first)
            height[selected],gradient[selected]=self.sample_piece(np.clip(local,0.,1.),first,last)
        return height.reshape(at.shape),gradient.reshape(at.shape+(2,))


class SourceOffsetRule:
    corner_fragment=False
    source_anchor=True

    def __init__(self,context,curve,first,last):
        self.context,self.curve=context,curve
        self.first,self.last=float(first),float(last)
        self.q=2.*(curve.offset-curve.approach.offset)/curve.approach.width
        if not 0.<=first<last<=1.+1e-9 or abs(abs(self.q)-1.)>1e-8:
            raise ValueError('Source offset rule needs a finite actual native bank')
        parameters=np.asarray(curve.parameters);cumulative=np.asarray(curve.cumulative)
        if (len(parameters)!=len(cumulative) or len(parameters)<2 or cumulative[0]!=0.
                or np.any(np.diff(cumulative)<=0.)
                or not (np.all(np.diff(parameters)>0.) or np.all(np.diff(parameters)<0.))):
            raise ValueError('Source lead needs its actual strictly monotone parameter map')
        knots=list(cumulative[1:-1]/curve.length)
        station_knots=[(float(value/curve.length),float(station))
                       for value,station in zip(cumulative,parameters)]
        path=curve.approach.path
        self.native=(type(path) in (Path,Route) and not hasattr(curve.approach,'field')
                     and not getattr(path,'mouth_piece',None)
                     and not getattr(path,'absorbed_sites',()))
        stations=[]
        from .paired_source_jet_domains import PairedSourceJetDomain
        self.paired_domain = None if self.native else PairedSourceJetDomain.create(context, curve.approach)
        if self.native:
            if type(path) is Route:
                for index,piece in enumerate(path.pieces):
                    values=piece.path.station(piece.path.knots)
                    stations.extend(path.stations[index]+(values if piece.forward else piece.path.length-values))
            else:
                values=path.station(path.knots)
                stations.extend(values if curve.approach.end==0 else path.length-values)
        elif self.paired_domain is not None:
            stations.extend(self.paired_domain.outward_knots())
        low,high=sorted((float(parameters[0]),float(parameters[-1])))
        for station in stations:
            if low<station<high:
                values,metric=(parameters,cumulative) if parameters[0]<parameters[-1] else (parameters[::-1],cumulative[::-1])
                value=float(np.interp(station,values,metric)/curve.length)
                knots.append(value);station_knots.append((value,float(station)))
        # These are the curve definition's own map seams, not the native
        # Path's adaptive integration lookup entries.
        self._breaks=tuple(sorted({float((value-first)/(last-first)) for value in knots if first<value<last}))
        self.ends=np.array((*self._breaks,1.))
        self._station_anchors={0.:float(self._stations(0.)),1.:float(self._stations(1.))}
        for value,station in station_knots:
            if value==first:self._station_anchors[0.]=station
            if value==last:self._station_anchors[1.]=station
            if first<value<last:self._station_anchors[float((value-first)/(last-first))]=station

    def polynomial_breaks(self):return self._breaks

    def _stations(self,t):
        return np.interp((self.first+np.asarray(t)*(self.last-self.first))*self.curve.length,
                         self.curve.cumulative,self.curve.parameters)

    def _gradient(self,station,interval):
        from .source_curb_conditions import _native_gradient
        approach=self.curve.approach;path=approach.path
        if not self.native:
            paired = getattr(self, 'paired_domain', None)
            if paired is not None and min(interval) >= 0. and max(interval) <= approach.length:
                return paired.gradient(station, self.q, interval)
            return np.asarray([approach.physical_gradient(float(s),self.q) for s in np.asarray(station).reshape(-1)])
        midpoint=.5*sum(interval)
        if midpoint<0.:
            if max(interval)>0.:
                raise ValueError('Source offset auxiliary jet crosses an undeclared original zero station')
            # Approach.frame defines this original auxiliary construction
            # exactly as p0+s*T0+d*N0. Its physical covector is constant; no
            # negative station becomes a native body coordinate.
            _,tangent,normal=approach.frame(0.)
            basis=np.stack((tangent,.5*approach.width*normal))
            gradient=np.linalg.solve(basis,np.array((0.,float(approach.profile.derivative(self.q)))))
            return np.broadcast_to(gradient,np.asarray(station).shape+(2,)).reshape(-1,2).copy()
        if min(interval)<0.:
            raise ValueError('Source offset native jet crosses an undeclared original zero station')
        if type(path) is Route:
            index=int(np.clip(np.searchsorted(path.stations,midpoint,side='right')-1,0,len(path.pieces)-1))
            piece=path.pieces[index];native=piece.path
            at=np.asarray(station)-path.stations[index]
            ends=np.asarray(interval)-path.stations[index]
            forward=piece.forward
        else:
            native=path;at=np.asarray(station);ends=np.asarray(interval)
            forward=approach.end==0
        if not forward:at=native.length-at;ends=native.length-ends
        profile=self.context[Key('profile',native.section.id)]
        side=self.q*(1. if forward else -1.)
        return _native_gradient(native,profile,at,side,interval=ends)

    def sample_piece(self,t,first,last):
        boundaries=(0.,*self._breaks,1.)
        if (first not in boundaries or last not in boundaries or not first<last
                or any(first<value<last for value in self._breaks)):
            raise ValueError('Source offset jet cannot skip its original parameter seams')
        local=first+np.asarray(t)*(last-first)
        station=self._stations(local)
        interval=self._station_anchors[first],self._station_anchors[last]
        values=np.asarray(t)
        station=np.where(values==0.,interval[0],np.where(values==1.,interval[1],station))
        gradient=self._gradient(station,interval).reshape(np.asarray(t).shape+(2,))
        profile=self.curve.approach.profile
        height=np.full(np.asarray(t).shape,profile.datum+float(profile.value(self.q)))
        return height,gradient

    def sample(self,t):
        at=np.asarray(t);flat=at.reshape(-1)
        indices=np.clip(np.searchsorted(self.ends,flat,side='right'),0,len(self.ends)-1)
        boundaries=(0.,*self._breaks,1.);height=np.empty(len(flat));gradient=np.empty((len(flat),2))
        for index in np.unique(indices):
            mask=indices==index;first,last=boundaries[int(index):int(index)+2]
            local=np.clip((flat[mask]-first)/(last-first),0.,1.)
            height[mask],gradient[mask]=self.sample_piece(local,first,last)
        return height.reshape(at.shape),gradient.reshape(at.shape+(2,))
