"""Original adjacent member ports retained through a side-junction compound."""
from dataclasses import dataclass
from .model import End,Key


@dataclass(frozen=True)
class RetainedBundleThroughPair:
    bundle: int
    role: int
    cuts: tuple
    junction: object
    sections: tuple = ()
    bundle_end: None = None

    def report(self):
        return dict(bundle=self.bundle,member=self.role,sections=[],
            policy='RETAINED_ORIGINAL_MEMBER_THROUGH_PAIR',original_junction=self.junction.text(),
            complete_exterior_half_profiles=[dict(section=e.section_id,end=e.end) for e in self.cuts],
            bundle_body_consumed=False,endpoint_extrapolation=False)

    @property
    def token(self):
        return (self.bundle,self.role,('through',self.junction.kind,self.junction.id))


def prepare(context,bundle,consumed,approaches,members):
    """Prove both ports of an unconsumed complete original member passage.

    Pair only consecutive declared traversals meeting at the same owned
    physical junction. This grants no Bundle body interval and no new source
    connection. Final complete-mouth, curb-cap and surface checks still apply.
    """
    if len(bundle.members)!=2 or not all(member.complete for member in bundle.members):
        return ()
    sites={end:site for site in context.network.intersections for end in site.ends if len(site.ends)>2}
    by_end={End(a.source.id,a.end):a for a in approaches}
    result=[]
    for role,member in enumerate(bundle.members):
        for before,after in zip(member.path,member.path[1:]):
            # Ownership of another member interval is independent of this
            # surviving original passage. Never turn a consumed traversal
            # itself into an exposed through port.
            if before.section_id in consumed or after.section_id in consumed:
                continue
            ends=(End(before.section_id,1 if before.forward else 0),
                  End(after.section_id,0 if after.forward else 1))
            site=sites.get(ends[0])
            if (site is None or site.key not in members or sites.get(ends[1]) is not site or
                    any(end not in by_end for end in ends)):
                continue
            valid=True
            for end in ends:
                approach=by_end[end]
                expected=context[Key('profile',end.section_id)].oriented(end.end)
                if (not expected.half or approach.profile!=expected or
                        getattr(approach.path,'absorbed_sites',()) or hasattr(approach.path,'mouth_piece')):
                    valid=False;break
                first=approach.path.pieces[0]
                if (End(first.section.id,0 if first.forward else 1)!=end or
                        any(piece.section.id in consumed for piece in approach.path.pieces) or
                        abs(approach.width-first.section.width)>1e-8 or
                        abs(approach.offset-first.section.offset*(1 if first.forward else -1))>1e-8):
                    valid=False;break
            if valid:result.append(RetainedBundleThroughPair(bundle.id,role,ends,site.key))
    return tuple(result)
