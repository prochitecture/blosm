"""Finite original cross-sections constrain near-through circular contacts."""
import math
import numpy as np
from .corner_patch import corner_patch, _stations
from .model import Key
from .regions import Region, RegionUnion, signed_area


def _canonical_root_stations(samples,extent,contact):
    fixed=[]
    for value in sorted({0.,extent,float(np.clip(contact,0.,extent))}):
        if not fixed or value-fixed[-1]>1e-9:fixed.append(value)
        elif value==extent:fixed[-1]=value
    result=[]
    for value in sorted(np.r_[np.clip(samples,0.,extent),fixed]):
        nearest=min(fixed,key=lambda compulsory:abs(compulsory-value))
        if abs(nearest-value)<=1e-9:value=nearest
        if not result or value-result[-1]>1e-9:result.append(float(value))
    return np.asarray(result)


def complete_source_roots(corner, spacing, tolerance, *, support_approaches=(), audit_all_contacts=False,
                          audit_facing_halves=False,on_refusal=None):
    """Keep the existing candidate when its local source roots stay covered.

    Tiny positive contacts on curved, nearly opposite roads can put a sliver
    of the original full root cross-section outside the prescribed circle.
    The local patch and finite authored strips at the same original node prove feasibility;
    fitting can then try the next already-permitted STYLE candidate. Negative
    contacts use the existing complete source-curb leads.

    The caller supplies support_approaches from one declared original junction
    incidence set. The origin check is a consistency guard, not a proximity
    search or a way to connect unrelated source roads.
    """
    if type(audit_all_contacts) is not bool:
        raise ValueError('Source-root contact auditing needs an explicit boolean policy')
    if type(audit_facing_halves) is not bool or (audit_facing_halves and not audit_all_contacts):
        raise ValueError('Candidate facing-half auditing requires forced contact auditing')
    if on_refusal is not None and not callable(on_refusal):
        raise ValueError('Source-root refusal diagnostics need an optional callable')
    if corner.center is None:
        return True
    if (not audit_all_contacts and (min(corner.stations) < 0. or
            min(corner.stations) >= .1*min(corner.first.width,corner.second.width) or
            abs(corner.gap-math.pi) > math.radians(15.))):
        return True
    patch=corner_patch(Key('root_patch',0),corner,spacing,tolerance)
    regions=[patch.region];caps={}
    targets=[(corner.first,corner.stations[0]),(corner.second,corner.stations[1])]
    seen={(a.source.id,a.end) for a,_ in targets}
    for approach in support_approaches:
        if np.linalg.norm(approach.origin-corner.first.origin)>1e-7:
            raise ValueError('Corner source-root support must use the same original physical junction')
        token=(approach.source.id,approach.end)
        if token not in seen:
            seen.add(token);targets.append((approach,0.))
    for index,(approach,contact) in enumerate(targets):
        extent=min(approach.length,max(2*approach.width,contact+.25))
        at=_stations(approach,extent,min(.5,spacing),tolerance,8192)
        # Negative stations belong only to the bounded auxiliary corner lead;
        # an original road strip always starts at its real source cross-section.
        at=_canonical_root_stations(at,extent,contact)
        points=np.vstack((approach.point(at,approach.offset+.5*approach.width),
                          approach.point(at,approach.offset-.5*approach.width)[::-1]))
        cap=points[[-1,0]]
        if signed_area(points)<0.:points=points[::-1]
        owner=Key('root_strip',index)
        regions.append(Region(owner,points))
        if index<2:caps[owner]=cap
    loops=RegionUnion().build(regions)
    by_owner={region.owner:region for region in regions}
    for loop in loops:
        for span in loop.spans:
            if span.owner not in caps:continue
            points=by_owner[span.owner].xy
            a,b=points[span.segment],points[(span.segment+1)%len(points)]
            cap=caps[span.owner]
            if ((np.linalg.norm(a-cap[0])<=1e-8 and np.linalg.norm(b-cap[1])<=1e-8) or
                    (np.linalg.norm(a-cap[1])<=1e-8 and np.linalg.norm(b-cap[0])<=1e-8)):
                if np.linalg.norm(b-a)*(span.last-span.first)>1e-7:
                    if not audit_facing_halves:
                        if on_refusal is not None:
                            direction=cap[1]-cap[0];squared=float(direction@direction)
                            endpoints=np.array((a+span.first*(b-a),a+span.last*(b-a)))
                            parameters=(endpoints-cap[0])@direction/squared
                            approach=corner.first if span.owner.id==0 else corner.second
                            on_refusal(dict(reason='EXPOSED_ORIGINAL_SOURCE_ROOT',section=approach.source.id,
                                end=approach.end,side=None,interval=tuple(map(float,sorted(parameters))),
                                exposed_length=float(np.linalg.norm(b-a)*(span.last-span.first))))
                        return False
                    # Each prescribed corner owns first.LEFT and second.RIGHT.
                    # Other root halves belong to other local corner sectors;
                    # the final full-source union still certifies all of them.
                    direction=cap[1]-cap[0]
                    squared=float(direction@direction)
                    endpoints=np.array((a+span.first*(b-a),a+span.last*(b-a)))
                    parameters=(endpoints-cap[0])@direction/squared
                    low,high=((0.,.5) if span.owner.id==0 else (.5,1.))
                    overlap=min(high,float(parameters.max()))-max(low,float(parameters.min()))
                    if overlap*np.sqrt(squared)>1e-7:
                        if on_refusal is not None:
                            approach=corner.first if span.owner.id==0 else corner.second
                            on_refusal(dict(reason='EXPOSED_ORIGINAL_FACING_ROOT_HALF',section=approach.source.id,
                                end=approach.end,side=-1 if span.owner.id==0 else 1,
                                interval=(max(low,float(parameters.min())),min(high,float(parameters.max()))),
                                exposed_length=float(overlap*np.sqrt(squared))))
                        return False
    return True


def raise_source_root_refusal(record):
    """Report the actual failed root contract, rather than invent a reach failure."""
    from .corner_radii import CornerContactError
    side={-1:'LEFT',1:'RIGHT',None:'FULL'}[record['side']]
    first,last=record['interval']
    raise CornerContactError(f"Original Source {record['section']}/end{record['end']} retains an exposed {side} root "
        f"interval [{first:.17g}, {last:.17g}] of length {record['exposed_length']:.17g} m")
