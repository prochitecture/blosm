"""One finite physical port for overlapping ground roads at a scope boundary.

Two original ground incidences can share a junction while their short retained
tails still overlap at distinct, recorded scope ends. There is then no finite
internal gore. Their *finite union* supplies the outgoing XY boundary, and one
explicit reconciled crown replaces the incompatible overlapping full crowns.
Neither source topology nor either original scope incidence is changed.
"""
from dataclasses import dataclass
import math
import numpy as np

from .model import End, Key
from .paths import right, unit
from .regions import Region, RegionUnion
from .source_mouths import SourceMouth
from .finite_containment import contains_segment as _contains_segment

MAXIMUM_LEAD_SAMPLES=20000
MAXIMUM_RIM_SEGMENTS=8192


@dataclass(frozen=True)
class FiniteScopeAlias:
    junction_end: End
    scope_end: End
    scope_station: float
    length: float


@dataclass(frozen=True)
class FiniteScopeInventory:
    pieces: tuple


class ReconciledScopeCrown:
    """Monotone crown halves with a zero jet at each actual outline bend.

    A bent finite mouth has no unique transverse tangent at a bend. Its new
    physical crown therefore has a stationary jet there, preserving a common
    height/gradient while never pretending to retain the old crossing crowns.
    """
    half=0
    peak=0.

    def __init__(self,rise,datum,parameters,edge_slopes):
        self.rise,self.datum=rise,datum
        self.knots=np.asarray(sorted(set((0.,.5,1.,*parameters))))
        self.values=rise*(1.-(2.*self.knots-1.)**2)
        self.end_segments = []
        breaks = set(self.knots)
        for side, slope in enumerate(edge_slopes):
            i = 0 if side == 0 else len(self.knots)-2
            length = self.knots[i+1]-self.knots[i]
            delta = abs(self.values[i+1]-self.values[i])
            if not np.isfinite(slope) or slope < -1e-9 or (delta == 0. and slope > 0.):
                raise ValueError('Finite scope crown cannot match an outward-rising source bank')
            slope = max(0., slope)
            # Integral of slope*(1-u/support)^4 is slope*support/5.
            # Reserve at least a tenth of the height for the ordinary smooth
            # rise across this segment. Halving prevents a support breakpoint
            # from producing a tiny sliver next to an existing outline bend.
            support = length
            if slope and slope*length > 4.5*delta:
                halves = math.ceil(math.log2(slope*length/(4.5*delta)))
                support = math.ldexp(length, -halves)
            self.end_segments.append((i, slope, support))
            if support < length:
                breaks.add(support if side == 0 else 1.-support)
        self.breakpoints = sorted(breaks)

    def _piece(self,q):
        u=.5*(np.asarray(q)+1.)
        index=np.clip(np.searchsorted(self.knots,u,side='right')-1,0,len(self.knots)-2)
        length=self.knots[index+1]-self.knots[index]
        return index,np.clip((u-self.knots[index])/length,0.,1.),length

    def value(self,q):
        i,t,length=self._piece(q)
        weight=t**3*(10.+t*(-15.+6.*t))
        value = self.values[i]+weight*(self.values[i+1]-self.values[i])
        for side,(segment,slope,support) in enumerate(self.end_segments):
            inward = t if side == 0 else 1.-t
            blend = weight if side == 0 else 1.-weight
            decay = np.maximum(0.,1.-inward*length/support)
            correction = slope*support/5.*(1.-decay**5-blend)
            value = value + np.where(i==segment, correction, 0.)
        return value

    def derivative(self,q):
        i,t,length=self._piece(q)
        blend = 30.*t*t*(1.-t)**2/length
        derivative = (self.values[i+1]-self.values[i])*blend
        for side,(segment,slope,support) in enumerate(self.end_segments):
            inward = t if side == 0 else 1.-t
            decay = np.maximum(0.,1.-inward*length/support)
            correction = slope*(decay**4 - support/5.*blend)
            derivative = derivative + np.where(i==segment, correction if side==0 else -correction, 0.)
        return .5*derivative

    def coordinate(self,side,v):
        return side*(1.-np.asarray(v))

    def intervals(self,side,base):
        return base


class FiniteScopePort:
    """A complete crown on the actual piecewise linear distal union boundary.

    This is deliberately not a Route and cannot extrapolate a source endpoint.
    The port has its own physical profile; aliases identify the consumed source
    intervals rather than claiming that two intersecting old crowns survived.
    """
    def __init__(self, first, second, points):
        self.original_approaches = (first, second)
        self.source, self.end = first.source, first.end
        self.path = FiniteScopeInventory(first.path.pieces + second.path.pieces)
        self.points = np.asarray(points, float)
        lengths = np.linalg.norm(np.diff(self.points, axis=0), axis=1)
        if np.any(lengths <= 1e-7):
            raise ValueError('Finite scope port has a collapsed physical boundary segment')
        self.stations = np.r_[0., np.cumsum(lengths)]
        self.width, self.offset, self.length = float(self.stations[-1]), 0., 0.
        self.parameters = self.stations / self.width
        self.origin = self.point(0., 0.)
        self.direction = unit(first.direction + second.direction)
        # The two exterior endpoints lie on the actual second and first
        # source banks, in that order. Match the full physical gradients
        # there. Interior bends of the finite mouth still have zero gradient.
        end_tangents = unit(np.diff(self.points, axis=0))[[0,-1]]
        gradients = np.asarray((second.physical_gradient(second.length,-1.),
                                first.physical_gradient(first.length,1.)))
        along = np.sum(gradients*end_tangents,axis=1)
        self.end_normal_gradients = gradients - along[:,None]*end_tangents
        self.profile = ReconciledScopeCrown(max(first.profile.rise, second.profile.rise),
            first.profile.datum,self.parameters,self.width*np.array((along[0],-along[1])))
        self.aliases = tuple(FiniteScopeAlias(End(a.source.id, a.end), End(a.source.id, 1-a.end),
            a.length if a.end == 0 else 0., a.length) for a in (first, second))

    def _parameter(self, station, distance):
        if np.any(abs(np.asarray(station)) > 1e-10):
            raise ValueError('A finite merged scope port cannot be extrapolated')
        u = .5 + np.asarray(distance) / self.width
        if np.any((u < -1e-10) | (u > 1.+1e-10)):
            raise ValueError('A finite merged scope port cannot leave its finite outline')
        return np.clip(u, 0., 1.)

    def point(self, station, distance):
        u = self._parameter(station, distance)
        return np.stack([np.interp(u, self.parameters, self.points[:, axis]) for axis in range(2)], axis=-1)

    def physical_gradient(self, station, q):
        u = self._parameter(station, .5*self.width*np.asarray(q))
        index = np.clip(np.searchsorted(self.parameters, u, side='right')-1, 0, len(self.points)-2)
        tangent = unit(np.diff(self.points, axis=0))[index]
        gradient = self.profile.derivative(q)[..., None]*tangent/(.5*self.width)
        for side,(_,_,support) in enumerate(self.profile.end_segments):
            inward = u if side == 0 else 1.-u
            decay = np.maximum(0.,1.-inward/support)
            gradient = gradient + decay[...,None]**4 * self.end_normal_gradients[side]
        return gradient

    def mouth_parameters(self, radial, spacing):
        # Every actual finite-union corner and the new complete crest survive.
        # Even intervals per physical segment preserve quad boundary parity.
        knots = self.profile.breakpoints
        result = []
        for first, last in zip(knots, knots[1:]):
            count = max(2, int(math.ceil((last-first)*self.width/spacing)))
            count += count % 2
            result.extend(np.linspace(first, last, count+1)[:-1])
        return np.r_[result, 1.][::-1]

    def source_mouths(self, station):
        self._parameter(station, 0.)
        # A physical half belongs to each source alias. These are identities
        # and parameter inventories, not independent incoming full crowns.
        return tuple(SourceMouth(alias.junction_end, alias.scope_station,
            float(role == 0), .5, float(role == 0), .5)
            for role, alias in enumerate(self.aliases))

    def report(self):
        return dict(mode='finite_ground_scope_union', profile='reconciled_complete_crown',
            datum=self.profile.datum, rise=self.profile.rise, outline=self.points.tolist(),
            source_aliases=[dict(section=a.junction_end.section_id, junction_end=a.junction_end.end,
                scope_end=a.scope_end.end, scope_station=a.scope_station, consumed_interval=[0., a.length])
                for a in self.aliases], endpoint_extrapolation=False,
            source_containment='complete_finite_linear_footprints')


def _straight_finite(context, approach):
    from .profiles import memberships
    from .routes import Route
    a = approach
    if (not isinstance(a.path, Route) or len(a.path.pieces) != 1 or a.profile.half
            or a.source.id in memberships(context.network)
            or a.source.selection_ends or 1-a.end not in a.source.scope_ends
            or abs(a.offset) > 1e-9 or getattr(a.path, 'absorbed_sites', ())
            or getattr(a.path, 'minimum_mouth', 0.) > 1e-9
            or a.length/.2 > MAXIMUM_LEAD_SAMPLES):
        return False
    points, tangents, normals = a.frame(np.linspace(0., a.length, 9))
    expected = a.origin + np.linspace(0., a.length, 9)[:, None]*a.direction
    return bool(np.max(np.linalg.norm(points-expected, axis=1)) < 1e-8
        and np.max(np.linalg.norm(tangents-a.direction, axis=1)) < 1e-8
        and np.max(np.linalg.norm(normals-right(a.direction), axis=1)) < 1e-8)


def _distal_union(first, second):
    """Exact finite straight-strip union; its distal chain has no extension."""
    regions=[]
    for a in (first, second):
        points=np.array((a.point(0., .5*a.width), a.point(a.length, .5*a.width),
            a.point(a.length, -.5*a.width), a.point(0., -.5*a.width)))
        regions.append(Region(Key('scope_strip', a.source.id), points))
    loops=RegionUnion().build(regions)
    if len(loops) != 1 or loops[0].hole:
        return None
    xy=loops[0].xy
    # first is the clockwise arm, second the counterclockwise arm.
    endpoints=(second.point(second.length, -.5*second.width), first.point(first.length, .5*first.width))
    indices=[]
    for point in endpoints:
        index=int(np.argmin(np.linalg.norm(xy-point, axis=1)))
        if np.linalg.norm(xy[index]-point)>1e-7:
            return None
        indices.append(index)
    candidates=[]
    axis=unit(first.direction+second.direction)
    for step in (1,-1):
        ids=[indices[0]]
        while ids[-1] != indices[1]:
            ids.append((ids[-1]+step)%len(xy))
            if len(ids)>len(xy):raise ValueError('Finite scope union failed to close')
        chain=xy[ids]
        candidates.append((float(np.min((chain-first.origin)@axis)), chain))
    score, points=max(candidates, key=lambda item:item[0])
    # The distal cap must lie strictly beyond the shared source junction.
    if score<=max(first.width,second.width):
        return None
    return points


def _overlap_to_terminal(first,second):
    """Exact proximity bound for a pair already certified to be straight.

    Distance to a finite line segment is convex on the first line interval.
    The shared start has zero distance, so its opposite endpoint is the only
    remaining maximum. The exact finite strip union still supplies the rim;
    this bound alone never authorizes a surface or endpoint extension.
    """
    point=first.frame(first.length)[0]
    station=float(np.clip((point-second.origin)@second.direction,0.,second.length))
    distance=float(np.linalg.norm(point-second.frame(station)[0]))
    return distance<.5*(first.width+second.width)-1e-7


def from_failure(context, junction, approaches, failure):
    """Bounded fallback for one failed pair with *recorded* finite scope ends."""
    from .finite_scope_spurs import from_failure as scope_spur
    expanded=scope_spur(context,junction,approaches,failure)
    if expanded is not None:return expanded
    from .intersections import CornerConstructionError, BoundaryCurve, JunctionLayout, Mouth, half_intervals
    from .corner_radii import CornerRadiusPolicy
    if not isinstance(failure, CornerConstructionError) or len(failure.approaches)!=2:
        return None
    first,second=failure.approaches
    if (junction.corner_radius(End(first.source.id,first.end))[1] != 'STYLE'
            or not all(_straight_finite(context,a) for a in (first,second))
            or np.linalg.norm(first.origin-second.origin)>1e-8
            or abs(first.profile.datum-second.profile.datum)>1e-9
            or first.direction@second.direction<=0.
            or not _overlap_to_terminal(first,second) or not _overlap_to_terminal(second,first)):
        return None
    ordered=sorted(approaches,key=lambda a:math.atan2(a.direction[1],a.direction[0])%math.tau)
    index=ordered.index(first)
    if ordered[(index+1)%len(ordered)] is not second:
        return None
    ordered=ordered[index:]+ordered[:index]
    outline=_distal_union(first,second)
    if outline is None:return None
    port=FiniteScopePort(first,second,outline)
    corners=[]
    pairs=[(second,ordered[2]), *list(zip(ordered[2:], ordered[3:])), (ordered[-1],first)]
    settings=context.network.settings
    for a,b in pairs:
        radius,mode=junction.corner_radius(End(a.source.id,a.end))
        corners.append(CornerRadiusPolicy.solve(a,b,radius,mode,
            math.radians(settings.continuation_angle_degrees)))
    external=ordered[2:]
    distances=[max(0.,corners[i].stations[1],corners[i+1].stations[0])+min(.25,a.width*.05)
        for i,a in enumerate(external)]
    if any(distance>=a.length-1e-6 for a,distance in zip(external,distances)):
        return None
    if (first.length+second.length+2.*sum(distances))/.2 > MAXIMUM_LEAD_SAMPLES:
        return None
    # Keep all original external circular contacts. Only their finite leads
    # reach the actual far source ends instead of an impossible internal gore.
    boundaries=[BoundaryCurve(corners[0],(second.length,distances[0]),settings.curve_tolerance)]
    boundaries.extend(BoundaryCurve(corners[i+1],(distances[i],distances[i+1]),settings.curve_tolerance)
        for i in range(len(external)-1))
    boundaries.append(BoundaryCurve(corners[-1],(distances[-1],first.length),settings.curve_tolerance))
    radial=half_intervals(settings)
    rows=4*int(math.ceil(max(2*radial,*[b.length/min(settings.spacing,
        2.*math.sqrt(max(1e-12,2.*b.corner.radius*settings.curve_tolerance-settings.curve_tolerance**2)))
        for b in boundaries])/4.))
    if (rows+2*radial)*(len(external)+1)+16 > MAXIMUM_RIM_SEGMENTS:
        raise ValueError('Finite merged scope port exceeds its physical rim sample budget')
    mouths={}
    for a,distance in zip(external,distances):
        for identity,station in a.path.trims(distance).items():mouths[identity]=Mouth(float(station),2*radial)
    for alias in port.aliases:
        mouths[alias.junction_end.section_id,alias.junction_end.end]=Mouth(alias.scope_station,2*radial)
    return JunctionLayout((port,*external),tuple(corners),tuple(boundaries),(0.,*distances),mouths,
        rows,radial,np.zeros(2),[0.]*(1+len(external)),mesh_strategy='boundary_chart',finite_scope_ports=(port,))


def groups(context, layouts):
    """Qualify physical finite scope footprints before normal surface dispatch."""
    from .compounds import CompoundGroup, SharedClaim, sample_footprint
    result=[]
    for owner,layout in layouts.items():
        ports=getattr(layout,'finite_scope_ports',())
        if not ports:continue
        footprint=sample_footprint(owner,layout)
        support=tuple(port.support_footprint(context,owner,layout) for port in ports
                      if hasattr(port,'support_footprint'))
        footprints=(footprint,*support)
        loops=RegionUnion().build(tuple(f.region for f in footprints))
        if len(loops)!=1 or loops[0].hole:
            raise ValueError('Finite merged scope port needs one simple physical exterior')
        from .compounds import rim_label
        sources={f.region.owner:f for f in footprints}
        if any(rim_label(span,sources).kind=='internal' for loop in loops for span in loop.spans):
            raise ValueError('Finite scope support exposed an internal original source cross-section')
        claims=[]
        for port in ports:
            for approach,alias in zip(port.original_approaches,port.aliases):
                # These straight original intervals end *on* the explicit
                # scope port, so strict interior containment is inappropriate.
                if not _contains_interval(approach,loops[0].xy):
                    raise ValueError('Finite scope union lost an original source interval')
                rectangle=np.array((approach.point(0.,.5*approach.width),
                    approach.point(approach.length,.5*approach.width),
                    approach.point(approach.length,-.5*approach.width),
                    approach.point(0.,-.5*approach.width)))
                if not all(_contains_segment(first,last,loops[0].xy)
                           for first,last in zip(rectangle,np.roll(rectangle,-1,axis=0))):
                    raise ValueError('Finite scope union lost an original finite source footprint')
                claims.append(SharedClaim(alias.junction_end.section_id,owner,owner,alias.length,0.))
            validate=getattr(port,'validate_group',None)
            if validate is not None:validate(context,loops[0])
            additional=getattr(port,'additional_claims',None)
            if additional is not None:claims.extend(additional(owner))
        result.append(CompoundGroup((owner,),tuple(claims),footprints,loops,physical_layout=layout))
    return tuple(result)


def _contains_interval(approach,polygon):
    """Certify every interval of one straight source against a finite polygon.

    Split at all boundary intersections and classify each open interval. This
    catches an arbitrarily short exit/re-entry, independently of sample spacing.
    Original terminal endpoints may lie on the explicit open scope boundary.
    """
    first,last=approach.frame(np.array((0.,approach.length)))[0]
    return _contains_segment(first,last,polygon)


