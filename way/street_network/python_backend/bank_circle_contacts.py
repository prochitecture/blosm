"""Finite contacts between an exposed Section bank and an original circular corner."""
from dataclasses import dataclass
import math
import numpy as np

from .compound_curbs import SamplingBudgetCurb
from .generated_junction_curbs import GeneratedJunctionArea
from .retained_arc_fillets import CircularBoundaryGeometry, OriginalAnalyticalCircleCurve
from .sampling import AnalyticalContactError
from .source_curb_coordinates import CurbCoordinates


class OriginalCornerArc(OriginalAnalyticalCircleCurve):
    """Analytical tangent on the circular portion of an original prepared curb."""

    def derivative(self, parameter):
        t=np.asarray(parameter)
        phase=self.circle_parameter(t)
        scale=self.corner.radius*self.corner.sweep
        if self.circle_interval is not None:
            knots=np.asarray(self.circle_interval.prepared_parameters)
            rates=np.diff(self.circle_interval.circle_parameters)/np.diff(knots)
            indices=np.clip(np.searchsorted(knots,t,side='right')-1,0,len(rates)-1)
            scale=np.asarray(scale*rates[indices])[...,None]
        angle = self.corner.start_angle + self.corner.sweep * phase
        return scale * np.stack((-np.sin(angle), np.cos(angle)), axis=-1)


class ContactSampling(SamplingBudgetCurb):
    def __init__(self, curve, tolerance, contact):
        super().__init__(curve, tolerance)
        self.contact = contact

    def certified_parameters(self, sampler, first, last):
        parameters = super().certified_parameters(sampler, first, last)
        if first < self.contact < last:
            parameters = np.union1d(parameters, [self.contact])
        return parameters


class CircularJoin(GeneratedJunctionArea):
    """An exact, sampled circle between two finite, identified source contacts."""
    leads = ()
    continuation = None

    @property
    def length(self):
        return self.corner.radius * abs(self.corner.sweep)

    @property
    def lengths(self):
        return np.asarray((0., 0., self.length, self.length))

    @property
    def curvature_bound(self):
        return 1. / self.corner.radius

    @property
    def contact_reference_points(self):
        return self(self.source_reference_parameters)

    def __call__(self, t):
        t = np.asarray(t)
        angle = self.corner.start_angle + self.corner.sweep*t
        xy = np.asarray(self.corner.center) + self.corner.radius*np.stack((np.cos(angle), np.sin(angle)), axis=-1)
        return np.where((t == 0.)[..., None], self.first, np.where((t == 1.)[..., None], self.last, xy))

    def tangent(self, t):
        angle = self.corner.start_angle + self.corner.sweep*np.asarray(t)
        return np.sign(self.corner.sweep)*np.stack((-np.sin(angle), np.cos(angle)), axis=-1)

    def certified_parameters(self, sampler, first, last):
        error = min(self.source_chord_tolerance, sampler.tolerance)
        step = min(sampler.spacing, math.sqrt(8.*error*self.corner.radius))
        knots = np.asarray(self.source_reference_parameters)
        return np.union1d(sampler.regular(first, last, step/self.length), knots[(first < knots) & (knots < last)])


@dataclass(frozen=True)
class BankCircleFillet(CircularJoin):
    owner: object
    corner: CircularBoundaryGeometry
    first: tuple
    last: tuple
    native_contact: object
    native_owner: object
    native_index: int
    native_parameter: float
    parent_owner: object
    parent_index: int
    parent_prepared: object
    parent_corner: object
    parent_parameter: float
    parent_interval: tuple
    bank_interval: tuple
    original_sectors: tuple
    closure_ancestry: tuple
    source_reference_parameters: tuple
    source_chord_tolerance: float

    def contact(self, source_owner, index, curve):
        if (source_owner, index) == (self.parent_owner, self.parent_index):
            parent = curve.curve if isinstance(curve, ContactSampling) else curve
            if (not isinstance(parent, OriginalAnalyticalCircleCurve)
                    or parent.original_prepared_curve is not self.parent_prepared
                    or parent.corner is not self.parent_corner):
                raise AnalyticalContactError('Bank fillet lost its original circular parent')
            return 1., self.parent_parameter
        if (source_owner, index) == (self.native_owner, self.native_index):
            coordinates = getattr(curve, 'source_coordinates', None)
            bank = self.native_contact
            if (not isinstance(coordinates, CurbCoordinates) or coordinates.path is not bank.path
                    or not np.any((coordinates.stations == bank.station) & (coordinates.sides == bank.side))):
                raise AnalyticalContactError('Bank fillet lost its original Section contact')
            coordinates.validate(curve.points)
            return 0., self.native_parameter
        raise AnalyticalContactError('Bank fillet meets a boundary outside its original source pair')

    def report(self):
        return dict(owner=self.owner.text(), native_section=self.native_contact.path.section.id,
                    native_side=self.native_contact.side, native_station=self.native_contact.station,
                    parent_owner=self.parent_owner.text(), parent_index=self.parent_index,
                    parent_contact_parameter=self.parent_parameter, parent_radius=self.parent_corner.radius,
                    radius=self.corner.radius, center=self.corner.center, sweep=self.corner.sweep,
                    actual_original_sectors=[(end.section_id,end.end,radius,mode)
                                             for end,radius,mode in self.original_sectors])


def shared_contact(sampler, curves, spans, labels, layouts):
    fillets = [i for i, curve in enumerate(curves) if isinstance(curve, CircularJoin)]
    if len(fillets) != 1:
        raise AnalyticalContactError('One bank/circle fillet must meet one original source')
    i, = fillets
    j = 1 - i
    original = np.array((labels[0].last, labels[1].first))
    parameters = original.copy()
    parameters[i], parameters[j] = curves[i].contact(spans[j].owner, labels[j].index, curves[j])
    reach = np.array([2./layouts[span.owner].rows for span in spans])
    if np.any(abs(parameters-original) > reach) or np.any((parameters < 0.) | (parameters > 1.)):
        raise AnalyticalContactError('Bank/circle fillet contact leaves its sampled neighborhood')
    xy = np.asarray([curve(t) for curve, t in zip(curves, parameters)])
    precision = max(1e-11, 32.*np.spacing(max(1., float(np.max(abs(xy))))))
    if np.linalg.norm(xy[0]-xy[1]) > precision:
        raise AnalyticalContactError('Bank and circular parent disagree at their declared contact')
    sampler._check_contact_budget(curves, parameters, spans, layouts)
    return parameters


def height_rule(resolver, source, index, curve, first, last):
    """Recover the actual source profiles and inward crossfalls at both contacts."""
    from .displaced_corner_authority import original_parent_corner
    from .model import Key
    from .paths import right
    from .source_bank_reconciliation import owned_ranges
    from .source_curb_conditions import SourceCurbConditionUnavailable, _SourceRule, _corner_slopes, _native_gradient
    from .source_parent_chain_conditions import endpoint_slopes
    from .source_parent_chain_fillets import radius_candidates

    if (index != 0 or source.region.owner != curve.owner
            or getattr(source.layout,'generated_junction_area',None) is not curve):
        raise SourceCurbConditionUnavailable('Bank/circle height rule lost its generated boundary')
    context = resolver.context
    bank = curve.native_contact
    sid = bank.path.section.id
    ranges = owned_ranges(context, resolver.group).get(sid, ())
    if (context[Key('path',sid)] is not bank.path or bank.owned_intervals != ranges
            or not any(a <= bank.station <= b for a,b in ranges)
            or not curve.bank_interval[0] < bank.station < curve.bank_interval[1]):
        raise SourceCurbConditionUnavailable('Circular join lost its actual finite Section bank')
    for owner, idx in ((curve.native_owner,curve.native_index),(curve.parent_owner,curve.parent_index)):
        curve.contact(owner, idx, resolver.sources[owner].layout.boundaries[idx])
    parent = resolver.sources[curve.parent_owner]
    original = curve.parent_corner
    site = getattr(parent.layout,'original_junction',None)
    node = next((v for v in context.network.intersections if v.key == site), None)
    phase = curve.parent_parameter
    if (node is None or site not in resolver.group.members or not original_parent_corner(resolver,node,original)
            or not curve.parent_interval[0] < phase < curve.parent_interval[1]
            or not curve.original_sectors or any(end not in node.ends or mode != 'STYLE'
                or node.corner_radius(end) != (radius,mode) for end,radius,mode in curve.original_sectors)
            or curve.corner.radius not in radius_candidates(min(r for _,r,_ in curve.original_sectors),'STYLE')):
        raise SourceCurbConditionUnavailable('Circular join lost its original radius or profile authority')
    station_xy = bank.path.point(bank.station,bank.path.section.offset+.5*bank.path.section.width*bank.side)
    circle_xy = parent.layout.boundaries[curve.parent_index](phase)
    precision = max(1e-10, 32.*np.spacing(max(1.,float(np.max(np.abs(curve.corner.center))))))
    if max(np.linalg.norm(station_xy-curve.first),np.linalg.norm(circle_xy-curve.last)) > precision:
        raise SourceCurbConditionUnavailable('Circular join no longer meets its two source contacts')
    levels = [a.profile.datum+float(a.profile.value(side))
              for a,side in ((original.first,-1.),(original.second,1.))]
    if abs(levels[0]-levels[1]) > 1e-10:
        raise SourceCurbConditionUnavailable('A mixed-height parent circle needs a longitudinal contact grade')
    a,b = _corner_slopes(original)
    slope = a+(b-a)*phase**3*(10.+phase*(-15.+6.*phase))
    angle = original.start_angle+original.sweep*phase
    tangent = np.sign(original.sweep)*np.array((-math.sin(angle),math.cos(angle)))
    profile = context[Key('profile',sid)]
    native_gradient = _native_gradient(bank.path,profile,bank.station,bank.side,interval=curve.bank_interval)
    slopes = endpoint_slopes(curve,(native_gradient,-slope*right(tangent)))
    z = profile.datum+float(profile.value(bank.side))
    rule = _SourceRule(curve,first,last,z,levels[0],*slopes,corner_fragment=True)
    rule.source_anchor = False
    return rule
