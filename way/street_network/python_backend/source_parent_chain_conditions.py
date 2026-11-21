"""Actual source endpoint conditions for a typed finite parent-chain fillet."""
import numpy as np
from .model import Key,End
from .source_parent_chain_fillets import PreparedParentChainFillet,radius_candidates
from .compound_dimensional_corners import DimensionalCompositeCurb
from .source_owned_connector_fillets import NativeCurbContact,_policy
from .source_curb_coordinates import CurbCoordinates
from .source_bank_reconciliation import owned_ranges
from .source_circular_intervals import PreparedCircularInterval


def endpoint_slopes(curve,gradients):
    """A zero-along-arc blend must retain each actual inward endpoint jet."""
    from .source_curb_conditions import SourceCurbConditionUnavailable
    values=np.asarray(gradients);tangents=curve.tangent(np.array((0.,1.)))
    if values.shape!=(2,2) or not np.isfinite(values).all():
        raise SourceCurbConditionUnavailable('Parent-chain endpoints need finite actual physical source jets')
    inward=np.column_stack((-tangents[:,1],tangents[:,0]))
    along=np.sum(values*tangents,axis=1);normal=np.sum(values*inward,axis=1)
    if np.max(abs(along))>1e-8 or np.min(normal)<-1e-10:
        raise SourceCurbConditionUnavailable('Parent-chain endpoint jet is incompatible with its true arc tangent/inward normal')
    return tuple(map(float,np.linalg.norm(values,axis=1)))


def build_rule(resolver,source,index,curve,first,last):
    """Validate original authorities, then blend their genuine endpoint jets.

    The core's existing semantic rule supplies the first condition. The second
    is the actual native carriageway bank, with its own profile and physical
    frame. The generated arc remains unanchored; it creates no source End.
    """
    from .source_curb_conditions import SourceCurbConditionUnavailable,_SourceRule,_native_gradient
    from .junction_bank_transitions import original_corner
    from .retained_arc_fillets import OriginalAnalyticalCircleCurve
    def refuse(message):raise SourceCurbConditionUnavailable(message)
    if (type(curve) is not PreparedParentChainFillet or index!=0
            or getattr(source.layout,'generated_junction_area',None) is not curve
            or getattr(source.layout,'source_parent_chain_fillet',None) is not curve
            or source.layout.boundaries[index] is not curve or source.region.owner!=curve.owner):
        refuse('Parent-chain rule lost its typed generated source identity')
    if not np.isfinite((first,last)).all() or not 0.<=first<last<=1.:
        refuse('Parent-chain rule needs one finite increasing generated arc interval')
    core=curve.core_contact;parent=curve.parent_circle;native=curve.native_contact
    if type(native) is not NativeCurbContact or not native.owned_intervals:
        refuse('Parent-chain rule needs a sealed actual native contact')
    source_core=resolver.sources.get(core.owner);source_circle=resolver.sources.get(parent.owner)
    if source_core is None or source_circle is None:
        refuse('Parent-chain rule lost an original physical source owner')
    if (type(core.curve) is not DimensionalCompositeCurb
            or not 0<=core.index<len(source_core.layout.boundaries)
            or source_core.layout.boundaries[core.index] is not core.curve
            or core.curve.corner is not core.original_corner
            or not np.isfinite((core.parameter,*core.visible_parameters)).all()
            or not 0.<=core.visible_parameters[0]<core.parameter<core.visible_parameters[1]<=1.):
        refuse('Parent-chain core changed its actual certified curve or finite original parameter')
    if (not 0<=parent.index<len(source_circle.layout.boundaries)
            or type(parent.interval) is not PreparedCircularInterval
            or parent.interval is not getattr(parent.curve,'original_circle_interval',None)
            or parent.interval.original_corner is not parent.corner):
        refuse('Parent-chain circle lost its original producer interval authority')
    bound=source_circle.layout.boundaries[parent.index]
    if (bound is not parent.curve and not (type(bound) is OriginalAnalyticalCircleCurve
            and bound.original_prepared_curve is parent.curve and bound.corner is parent.corner)):
        refuse('Parent-chain circle changed its actual original prepared source')
    parent.interval.validate(parent.curve,parent.corner)
    policies=[]
    for original,site,corner in ((source_core,core.original_site,core.original_corner),
                                  (source_circle,parent.original_site,parent.corner)):
        node=next((item for item in resolver.context.network.intersections if item.key==site),None)
        if (site not in resolver.group.members or getattr(original.layout,'original_junction',None)!=site
                or node is None or not original_corner(resolver,node,corner)):
            refuse('Parent-chain rule changed an original junction/profile/radius authority')
        endpoint=End(corner.first.source.id,corner.first.end);radius,mode=node.corner_radius(endpoint)
        policies.append((site,endpoint,radius,mode))
    policy=_policy(policies)
    allowed=radius_candidates(curve.preferred_radius,curve.radius_mode)
    if (policy is None or policy!=(curve.preferred_radius,curve.radius_mode)
            or not curve.attempted_radii or curve.attempted_radii[0]!=curve.preferred_radius
            or curve.attempted_radii[-1]!=curve.radius or len(curve.attempted_radii)>17
            or tuple(curve.attempted_radii)!=allowed[:len(curve.attempted_radii)]
            or curve.radius_mode=='STRICT' and curve.radius!=curve.preferred_radius):
        refuse('Parent-chain rule changed its original STYLE/STRICT radius policy')
    sid=native.path.section.id;ranges=owned_ranges(resolver.context,resolver.group).get(sid,())
    if (native.path is not resolver.context[Key('path',sid)] or native.side not in (-1,1)
            or ranges!=native.owned_intervals
            or not any(a<=native.station<=b for a,b in ranges)):
        refuse('Parent-chain rule changed its sealed actual native ownership')
    tolerance=max(1e-10,32.*np.spacing(max(1.,float(np.max(abs(np.asarray(curve.center)))))))
    angles=curve.start_angle+curve.sweep*np.array((0.,1.))
    endpoints=np.asarray(curve.center)+curve.radius*np.column_stack((np.cos(angles),np.sin(angles)))
    if np.max(np.linalg.norm(endpoints-np.asarray((curve.first,curve.last)),axis=1))>tolerance:
        refuse('Parent-chain rule changed its analytical circle/native endpoint agreement')
    if np.linalg.norm(core.curve(core.parameter)-np.asarray(curve.first))>tolerance:
        refuse('Parent-chain rule lost its exact certified core contact')
    found=[]
    for native_source in resolver.sources.values():
        if native_source.region.owner.kind!='cycle_strip' or native_source.region.owner.id!=sid:continue
        for bank in native_source.layout.boundaries:
            data=getattr(bank,'source_coordinates',None)
            if type(data) is not CurbCoordinates or data.path is not native.path or data.approach is not None:continue
            data.validate(bank.points,resolver.context)
            for at in np.flatnonzero((data.stations==native.station)&(data.sides==native.side)):
                found.append(bank.points[at])
    if not found or any(np.linalg.norm(point-np.asarray(curve.last))>tolerance for point in found):
        refuse('Parent-chain rule lost its exact original native contact knot')
    expected=native.path.point(native.station,native.path.section.offset+.5*native.path.section.width*native.side)
    if np.linalg.norm(expected-np.asarray(curve.last))>tolerance:
        refuse('Parent-chain native contact changed its actual profile/frame')
    labels=[(i,label) for i,label in enumerate(source.segments) if label.kind=='curb' and label.index==0]
    values=np.asarray([label.first for i,label in labels]+[labels[-1][1].last]) if labels else np.array(())
    points=np.asarray([source.region.xy[i] for i,label in labels]+[source.region.xy[(labels[-1][0]+1)%len(source.region.xy)]]) if labels else np.empty((0,2))
    if (not labels or tuple(values)!=curve.source_reference_parameters
            or np.max(np.linalg.norm(curve(values)-points,axis=1))>tolerance):
        refuse('Parent-chain generated arc changed its original prepared knots')
    rule=resolver.source_rule(core.owner,core.index,core.parameter,min(1.,core.parameter+1e-5))
    z0,g0=rule.sample(0.);profile=resolver.context[Key('profile',sid)]
    z1=profile.datum+float(profile.value(native.side));g1=_native_gradient(native.path,profile,native.station,native.side)
    slopes=endpoint_slopes(curve,(g0,g1))
    result=_SourceRule(curve,first,last,float(z0),z1,*slopes,corner_fragment=True)
    result.source_anchor=False
    return result
