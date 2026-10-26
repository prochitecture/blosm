"""Bounded Exact-seam fitting with certified chord error and explicit clearance."""
import math
import numpy as np
from .exact_transitions import ExactTransitionGeometry, restricted
from .model import End


class ExactRowSampler:
    """Uniform rows per mandatory piece, certified by Bernstein convex hulls."""
    def __init__(self, curbs, tolerance, max_intervals):
        if not math.isfinite(tolerance) or tolerance <= 0.:
            raise ValueError('Exact fitted rows need a positive finite curve tolerance')
        self.curbs, self.tolerance, self.max_intervals = tuple(curbs), tolerance, max_intervals
        self.knots = self.curbs[0].knots
        if any(not np.array_equal(c.knots, self.knots) for c in self.curbs):
            raise ValueError('Exact fitted curbs need identical mandatory knots')
        if 2*(len(self.knots)-1) > max_intervals:
            raise ValueError('Exact fitted rows exceed their allocation budget')
        counts = []
        for i in range(len(self.knots)-1):
            pieces = tuple(c.segments[i] for c in self.curbs)
            second = max(p.maximum_second for p in pieces)
            upper = max(2, int(math.ceil(math.sqrt(second/(8.*tolerance)))))
            upper += upper % 2
            limit = max_intervals-max_intervals % 2
            upper = min(upper, limit)
            # The second-derivative estimate seeds the search. A control-hull
            # certificate can be more conservative than the true curve error.
            while self.error(pieces, upper) > tolerance:
                if upper >= limit:
                    raise ValueError('Exact fitted curve tolerance exceeds its sampling budget')
                upper = min(2*upper,limit)
            low, high = 1, upper//2
            while low < high:
                mid = (low+high)//2
                if self.error(pieces, 2*mid) <= tolerance:
                    high = mid
                else:
                    low = mid+1
            counts.append(2*low)
        if sum(counts) > max_intervals:
            raise ValueError('Exact fitted rows exceed their allocation budget')
        self.counts = tuple(counts)

    @staticmethod
    def error(pieces, count):
        maximum = 0.
        parameters = np.linspace(0., 1., count+1)
        for piece in pieces:
            # Covers the small fixed-degree de Casteljau arithmetic as well as
            # the Bernstein residual hull; world translation is not ignored.
            roundoff = 64.*np.finfo(float).eps*max(1., float(abs(piece.controls).max()))
            for a,b in zip(parameters, parameters[1:]):
                controls = restricted(piece.controls, float(a), float(b))
                linear = np.linspace(controls[0], controls[-1], len(controls))
                maximum = max(maximum, float(np.linalg.norm(controls-linear,axis=1).max())+roundoff)
        return maximum

    def candidates(self, resolution):
        seen = set()
        target = float(resolution)
        for scale in (target,max(1.,.75*target),max(1.,.5*target),1.,1.25*target,1.5*target,2.*target):
            counts = tuple(2*int(math.ceil(n*scale/2.)) for n in self.counts)
            if counts in seen or sum(counts) > self.max_intervals:
                continue
            seen.add(counts)
            bound = max(self.error(tuple(c.segments[i] for c in self.curbs),n) for i,n in enumerate(counts))
            if bound > self.tolerance:
                continue
            parameters = np.r_[np.concatenate([np.linspace(a,b,n+1)[:-1]
                for a,b,n in zip(self.knots,self.knots[1:],counts)]),1.]
            parameters.flags.writeable = False
            yield parameters, bound


class ExactTransitionFitter:
    """Fallback after the established transition branch cannot certify a mesh."""
    max_candidates = 64

    @staticmethod
    def free_limits(network, approaches, clearance):
        from .paths import Path
        from .routes import Route
        incidents = {end for site in (*network.intersections,*network.connections) for end in site.ends}
        bundled = {t.section_id for b in network.bundles for m in b.members for t in m.path}
        limits = np.array(clearance.limits,copy=True)
        records = []
        for i,approach in enumerate(approaches):
            path = approach.path
            if isinstance(path,Route):
                if not path.terminal_free or getattr(path,'mouth_piece',None) is not None:
                    continue
                ids = tuple(p.section.id for p in path.pieces)
                last = path.pieces[-1]
                far = End(last.section.id,int(last.forward))
            elif isinstance(path,Path):
                ids = (path.section.id,)
                far = End(path.section.id,1-approach.end)
                if far.end in path.section.selection_ends:
                    continue
            else:
                continue  # No explicit free-terminal contract on this path view.
            if (far in incidents or bundled.intersection(ids)
                    or any(c['section']==approach.source.id for c in clearance.constraints)):
                continue
            gap = min(.25,.05*approach.width)
            cap = min(.95*approach.length,approach.length-gap)
            if cap > limits[i]:
                limits[i] = cap
                records.append(dict(approach=i,section=approach.source.id,
                    free_end=dict(section_id=far.section_id,end=far.end),maximum_distance=float(cap),
                    remaining_full_profile_length=float(approach.length-cap)))
        limits.flags.writeable = False
        return limits, records

    @classmethod
    def fit(cls, context, site, approaches, clearance, desired, across, failures):
        from .transitions import SectionTransition, TransitionClearance, _certify_orientation
        limits, free = cls.free_limits(context.network,approaches,clearance)
        settings = context.network.settings
        seen, attempts = set(), 0
        for cap in (clearance.limits,limits):
            effective = TransitionClearance(cap,clearance.constraints)
            for factor in (1.,1.5,2.,3.,4.,6.,8.):
                distances = np.minimum(.5*desired*factor,cap)
                identity = tuple(distances)
                if identity in seen:
                    continue
                seen.add(identity)
                try:
                    curbs, geometry = ExactTransitionGeometry.build(approaches,distances,coherent_collinear=True)
                    orientation = _certify_orientation(*curbs)
                    sampler = ExactRowSampler(curbs,settings.curve_tolerance,50000//across)
                    for parameters,bound in sampler.candidates(settings.resolution):
                        attempts += 1
                        if attempts > cls.max_candidates:
                            raise ValueError('Exact transition fitting exhausted its candidate budget')
                        certification = dict(orientation,**geometry,sampling_policy='certified_curve_error',
                            certified_curve_error=bound,
                            sampling_piece_intervals=[int(np.count_nonzero((parameters>=a)&(parameters<=b))-1)
                                for a,b in zip(curbs[0].knots,curbs[0].knots[1:])],
                            extended_free_ends=[r for r in free if distances[r['approach']]>clearance.limits[r['approach']]],
                            original_trim_distances=list(map(float,clearance.limits)))
                        result = SectionTransition(site.key,tuple(approaches),tuple(map(float,distances)),
                            tuple(curbs),len(parameters)-1,across,parameters,certification)
                        try:
                            return SectionTransition.certify(context,result,effective)
                        except ValueError as error:
                            failures.append(str(error))
                except ValueError as error:
                    failures.append(str(error))
                if attempts >= cls.max_candidates:
                    return None
        return None
