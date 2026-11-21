"""Final prepared-owner exclusion for geometry-only junction infills."""
from copy import copy
from dataclasses import replace
import numpy as np
from .compounds import sample_footprint
from .generated_junction_curbs import is_generated_junction_area
from .source_bank_reconciliation import SourceBankReconciliationError
from .source_mouth_clearance import _positive_overlap


class PreparedNeighborGuard:
    """Use eventual physical owners, never substitute their raw input layout."""
    def __init__(self,context,groups,layouts,*,maximum_pairs=2000000):
        if type(maximum_pairs) is not int or not 1<=maximum_pairs<=2000000:
            raise ValueError('Generated curb neighbors need a finite pair allowance')
        self.context,self.groups,self.layouts=context,tuple(groups),layouts
        self.maximum_pairs=maximum_pairs
        self.area,self.area_owners,self.mouths=[],[],[]
        occupied=set()
        from .ordinary_mouth_clearance import _owner_mouths
        for group in self.groups:
            if occupied.intersection(group.members):
                raise SourceBankReconciliationError('Prepared physical neighbor inventory repeats an original owner')
            occupied.update(group.members)
            owner=tuple(group.members)
            self.mouths.append((owner,_owner_mouths(context,group).mouths))
            for source in group.area_footprints:
                self.area.append(source.region.xy);self.area_owners.append(owner)
        for key,layout in sorted(layouts.items()):
            if key in occupied:continue
            if not hasattr(layout,'mouths') or not hasattr(layout,'boundaries'):
                raise SourceBankReconciliationError('Generated curb guard needs every eventual physical neighbor owner')
            footprint=sample_footprint(key,layout)
            self.area.append(footprint.region.xy);self.area_owners.append((key,))
            self.mouths.append(((key,),layout.mouths))
        self.low=np.array([area.min(0) for area in self.area]).reshape(-1,2)
        self.high=np.array([area.max(0) for area in self.area]).reshape(-1,2)

    def validate(self):
        result=[];pairs=0
        from .ordinary_mouth_clearance import _owner_mouths
        band=2.*min(.001,.05*self.context.network.settings.curve_tolerance)+1e-7
        for group in self.groups:
            generated=tuple(source for source in group.area_footprints if is_generated_junction_area(source))
            if not generated:
                result.append(group);continue
            current=tuple(group.members);mouths=_owner_mouths(self.context,group).mouths
            for owner,neighbor in self.mouths:
                if owner==current:continue
                for (sid,end),mouth in mouths.items():
                    opposite=neighbor.get((sid,1-end))
                    if opposite is not None and (mouth.station>opposite.station+1e-7 if end==0
                                                else mouth.station<opposite.station-1e-7):
                        raise SourceBankReconciliationError('Generated curb conflicts with an eventual prepared neighbor mouth')
            examined=0
            for source in generated:
                polygon=source.region.xy
                selected=np.flatnonzero(np.all(polygon.max(0)>=self.low,axis=1)
                    &np.all(polygon.min(0)<=self.high,axis=1))
                for index in selected:
                    if self.area_owners[index]==current:continue
                    pairs+=1;examined+=1
                    if pairs>self.maximum_pairs:
                        raise SourceBankReconciliationError('Generated curb neighbor exclusion exceeds its finite pair budget')
                    if _positive_overlap(polygon,self.area[index],tolerance=band):
                        neighbor=','.join(key.text() for key in self.area_owners[index])
                        raise SourceBankReconciliationError('Generated curb overlaps prepared physical owner '+neighbor)
            proof=dict(policy='EVENTUAL_PREPARED_PHYSICAL_OWNER_EXCLUSION',
                original_members=[key.text() for key in group.members],generated_areas=len(generated),
                owner_inventory_count=len(self.mouths),area_operands=len(self.area),
                tested_pairs=examined,numerical_source_error_band=band,
                original_source_claims_and_mouths_unchanged=True)
            layout=copy(group.physical_layout)
            layout.generated_curb_prepared_neighbor_guard=proof
            from .source_root_templates import NativeRootEnvelope
            if any(isinstance(getattr(source.layout,'generated_junction_area',None),NativeRootEnvelope)
                   for source in generated):
                layout.source_root_prepared_neighbor_guard=proof
            record=getattr(layout,'source_bank_reconciliation',None)
            if record:
                layout.source_bank_reconciliation=dict(record,prepared_owner_area_guard=proof)
            arc=getattr(layout,'retained_arc_reconciliation',None)
            if arc:
                layout.retained_arc_reconciliation=dict(arc,prepared_owner_area_guard=proof)
            chain=getattr(layout,'source_parent_chain_reconciliation',None)
            if chain:
                layout.source_parent_chain_reconciliation=dict(chain,prepared_owner_area_guard=proof)
            gaps=getattr(layout,'residual_road_gap_reconciliation',None)
            if gaps:
                layout.residual_road_gap_reconciliation=dict(gaps,prepared_owner_area_guard=proof)
            result.append(replace(group,physical_layout=layout))
        return tuple(result)


def validate(context,groups,layouts):
    groups=tuple(groups)
    if not any(is_generated_junction_area(source) for group in groups for source in group.area_footprints):
        return groups
    return PreparedNeighborGuard(context,groups,layouts).validate()
