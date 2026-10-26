"""One-sided finite strips opposite an aligned, merged free Bundle end.

Source endpoints and profiles are unchanged. A free overlapping end supplies
the common chart direction; separate profile strips are created only at the
other end. The original full boundary remains the actual pavement footprint.
"""
import numpy as np
from .bundle_axis import AxialChartUnavailable
from .bundle_caps import FiniteAxialPlan
from .bundle_ends import member_ends
from .model import Key
from .paths import unit


class MergedFreeCapPlan(FiniteAxialPlan):
    coordinated_fibres = True
    @classmethod
    def build(cls, context, bundle):
        sites={e for j in (*context.network.intersections,*context.network.connections) for e in j.ends}
        errors=[]
        for end in (0,1):
            if any(member_ends(m)[end] in sites for m in bundle.members):
                continue
            corners,directions=[],[]
            for member in bundle.members:
                visit=member.path[0 if end==0 else -1]
                path=context[Key('path',visit.section_id)]
                source=path.section
                at=0. if (end==0)==visit.forward else path.length
                _,tangent,_=path.frame(at)
                directions.append(tangent*(1 if visit.forward else -1))
                corners.append(path.point(at,source.offset+.5*source.width*np.array((-1.,1.))))
            combined=np.asarray(directions).sum(axis=0)
            if np.linalg.norm(combined)<1e-8:
                continue
            axis=unit(combined)
            if any(np.linalg.norm(v-axis)>1e-6 for v in directions):
                continue
            points=np.vstack(corners)
            origin=points[0]
            if np.ptp((points-origin)@axis)>1e-6:
                continue
            left=np.array((-axis[1],axis[0]))
            values=[(p-origin)@left for p in corners]
            # Preserve the Bundle's declared lateral order and select only an
            # end whose inner curbs already touch or overlap.
            if values[0].mean()<=values[1].mean() or values[0].min()-values[1].max()>1e-8:
                continue
            try:
                return super().build(context,bundle,merged_end=end,axis=axis)
            except AxialChartUnavailable as error:
                errors.append(str(error))
        if errors:
            raise AxialChartUnavailable('No qualified one-sided merged free cap: '+'; '.join(errors))
        return None

    def __init__(self, context, bundle, *, margin_factor, merged_end, axis):
        self.merged_end=merged_end
        self.merged_axis=np.asarray(axis,dtype=float)
        super().__init__(context,bundle,margin_factor=margin_factor)

    def chart_axis(self, context, directions):
        return self.merged_axis.copy()

    def strip_ends(self):
        return (1-self.merged_end,)

    def pavement(self, context):
        part=super().pavement(context)
        part.metadata.update(merged_free_end=self.merged_end,finite_strip_ends=list(self.strip_ends()))
        return part
