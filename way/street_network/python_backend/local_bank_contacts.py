"""Fit a circular join on the exposed parts of two original Section banks."""
import math
import numpy as np

from .corner_radii import CornerContactError, CornerRadiusPolicy
from .curb_contact import CircularCurbContact
from .intersections import Approach, Corner
from .model import Key
from .source_bank_reconciliation import SourceBankDiscovery, visible_bank_stations
from .source_curb_rounding import FiniteSourceBank, SourceBankFillet
from .source_parent_chain_fillets import ParentChainUnavailable, bank_tangent, offset_solution


class DirectedBank:
    def __init__(self, bank, domain):
        self.path = bank.path
        self.first, self.last = domain
        self.offset = bank.path.section.offset+.5*bank.path.section.width*bank.side

    def station(self, t):
        return self.first+(self.last-self.first)*np.asarray(t)

    def __call__(self, t):
        return self.path.point(self.station(t),self.offset)

    def derivative(self, t):
        # The offset solver needs the direction of this regular locus.
        # An Exact road's miter vector can rotate along one affine segment.
        station = self.station(t)
        path = self.path
        if path.exact:
            i = int(np.clip(np.searchsorted(path.knots,.5*(self.first+self.last),side='right')-1,0,len(path.chords)-1))
            vector = path.directions[i] + self.offset*(path.normals[i+1]-path.normals[i])/path.chords[i]
            tangent = np.broadcast_to(vector,station.shape+(2,))
        else:
            tangent = path.frame(station)[1]
        return tangent*(self.last-self.first)


def _outward(arm, values):
    values = np.asarray(values)
    return values if arm.end == 0 else arm.length-values


def fit(context, group, proposal):
    discovery = SourceBankDiscovery(context,group)
    banks = proposal.first, proposal.last
    domains = tuple(visible_bank_stations(discovery,proposal,incoming=incoming) for incoming in (True,False))
    if any(domain is None for domain in domains):
        raise CornerContactError('Local fillet needs two continuous visible bank intervals')
    incoming, outgoing = (discovery.tangent(bank) for bank in banks)
    turn = math.atan2(float(incoming[0]*outgoing[1]-incoming[1]*outgoing[0]),float(incoming@outgoing))
    if not -math.pi < turn < -math.radians(1e-4):
        raise CornerContactError('The two original banks have no concave local turn')
    return _fit_on_intervals(context, proposal, domains, turn)


def _fit_on_intervals(context, proposal, domains, turn):
    banks = proposal.first, proposal.last
    arms = tuple(Approach(bank.path,end,context[Key('profile',bank.path.section.id)])
                 for bank,end in zip(banks,(int(banks[0].side == 1),int(banks[1].side == -1))))
    owned = tuple(FiniteSourceBank(bank.path,bank.side,bank.intervals) for bank in banks)
    curvature = []
    for bank, domain in zip(banks,domains):
        low,high = sorted(domain)
        if bank.path.exact and np.any((low < bank.path.knots) & (bank.path.knots < high)):
            raise CornerContactError('Local bank contact crosses a miter in the original polyline')
        curvature.append(0. if bank.path.exact else bank.path.visible_curvature_bound(low,high))
    # The source banks intersect before rounding. Restrict the crossing solve
    # to their actual owned intervals, seeded by the observed rim switch.
    bounds = [sorted(_outward(arm,next(interval for interval in bank.intervals
                    if interval[0] <= bank.station <= interval[1]))) for bank,arm in zip(banks,arms)]
    solver = CircularCurbContact(*arms,0.,1.)
    solver.tolerance = max(1e-11,8*np.finfo(float).eps*max(1.,float(np.max(np.abs(proposal.point)))))
    low,high = np.asarray(bounds).T
    crossing,xy,points = solver.solve([float(_outward(arm,bank.station)) for arm,bank in zip(arms,banks)],low,high)
    if (np.linalg.norm(points[0]-points[1]) > 1e-8
            or np.linalg.norm(xy-np.asarray(proposal.point)) > 2.*context.network.settings.curve_tolerance):
        raise CornerContactError('Solved bank crossing leaves the observed local source contact')
    first = DirectedBank(banks[0],domains[0])
    last = tuple(_outward(arms[1],domains[1]))
    offsets = [abs(bank.path.section.offset+.5*bank.path.section.width*bank.side) for bank in banks]
    attempts = []
    for radius in CornerRadiusPolicy.candidates(proposal.radius,proposal.mode):
        attempts.append(radius)
        if any((offset+radius)*k >= 1.-1e-10 for offset,k in zip(offsets,curvature)):
            continue
        try:
            result = offset_solution(first,(0.,1.),arms[1],last,radius)
        except ParentChainUnavailable:
            continue
        station = float(_outward(arms[0],first.station(result['poly_contact']))),result['bank_outward_contact']
        if not all(bank.contains(*_outward(arm,(s,cut))) for bank,arm,s,cut in zip(owned,arms,station,crossing)):
            continue
        tangents = np.stack((first.derivative(result['poly_contact']),bank_tangent(arms[1],station[1],last)))
        tangents /= np.linalg.norm(tangents,axis=1)[:,None]
        radial = np.asarray([result['first'],result['last']])-result['center']
        if np.max(np.abs(np.sum(radial*tangents,axis=1))) > 1e-8:
            continue
        corner = Corner(*arms, station,np.asarray(result['center']),radius,result['start_angle'],result['sweep'],
                        math.pi+turn,requested_radius=proposal.radius,radius_mode=proposal.mode,
                        radius_attempts=tuple(attempts),fit_policy='FINITE_NATIVE_LOCAL_OFFSET')
        return SourceBankFillet(corner,tuple(map(float,crossing)),owned)
    raise CornerContactError('No permitted radius meets both visible native banks with the local turn')
