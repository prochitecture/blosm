"""Producer-declared circle pieces inside unchanged prepared curb curves."""
from dataclasses import dataclass
import numpy as np

SOURCE_CIRCLE_INTERVAL_POLICY='1-original-prepared-circle-subintervals'


def _finite_map(first,last):
    return (len(first)>=2 and len(first)==len(last)
        and np.isfinite(np.asarray((first,last))).all()
        and all(0.<=value<=1. for value in first)
        and all(b>a for a,b in zip(first,first[1:])))


def _bounded_interp(value,xp,fp):
    values=np.asarray(value)
    if not np.isfinite(values).all() or np.any(values<xp[0]) or np.any(values>xp[-1]):
        raise ValueError('Circle piece parameter leaves its explicitly declared finite interval')
    return np.interp(values,xp,fp)


@dataclass(frozen=True)
class PreparedAuxiliaryLead:
    original_approach: object
    bank_side: int
    prepared_parameters: tuple
    outward_stations: tuple
    vertex_indices: tuple
    reference_points: tuple

    @property
    def prepared_bounds(self):return self.prepared_parameters[0],self.prepared_parameters[-1]

    def source_station(self,t):return _bounded_interp(t,self.prepared_parameters,self.outward_stations)

    def validate(self,curve,corner):
        expected=corner.first if self.bank_side==-1 else corner.second
        ids=np.asarray(self.vertex_indices);points=np.asarray(self.reference_points)
        if (type(self.bank_side) is not int or self.bank_side not in (-1,1)
                or self.original_approach is not expected
                or not _finite_map(self.prepared_parameters,self.outward_stations)
                or ids.dtype.kind not in 'iu' or np.any(ids<0) or np.any(ids>=len(curve.points))
                or len(ids)!=len(self.prepared_parameters)
                or not np.array_equal(curve.points[ids],points)
                or not np.array_equal(curve.stations[ids]/curve.length,self.prepared_parameters)):
            raise ValueError('Original auxiliary lead changed its declared source/vertex ancestry')
        arm=self.original_approach
        actual=arm.point(np.asarray(self.outward_stations),arm.offset+.5*arm.width*self.bank_side)
        if np.max(np.linalg.norm(actual-points,axis=1))>1e-8:
            raise ValueError('Declared auxiliary lead no longer evaluates on its original construction')
        return True


@dataclass(frozen=True)
class PreparedCircularInterval:
    original_corner: object
    prepared_parameters: tuple
    circle_parameters: tuple
    vertex_indices: tuple
    reference_points: tuple
    auxiliary_leads: tuple
    maximum_chord_error: float

    def __post_init__(self):
        phase=np.asarray(self.circle_parameters)
        if (not _finite_map(self.prepared_parameters,self.circle_parameters)
                or not np.isin(phase[[0,-1]],(0.,1.)).all()
                or phase[0]==phase[-1] or not (np.all(np.diff(phase)>0.) or np.all(np.diff(phase)<0.))
                or len(self.vertex_indices)!=len(phase) or len(self.reference_points)!=len(phase)
                or not np.isfinite(self.maximum_chord_error) or self.maximum_chord_error<=0.):
            raise ValueError('Prepared circle piece needs exact finite original arc-knot ancestry')

    @property
    def prepared_bounds(self):return self.prepared_parameters[0],self.prepared_parameters[-1]

    @property
    def lead_domains(self):return tuple(lead.prepared_bounds for lead in self.auxiliary_leads)

    def circle_parameter(self,t):return _bounded_interp(t,self.prepared_parameters,self.circle_parameters)

    def source_parameter(self,t):
        xp,fp=self.circle_parameters,self.prepared_parameters
        if xp[-1]<xp[0]:xp,fp=xp[::-1],fp[::-1]
        return _bounded_interp(t,xp,fp)

    def validate(self,curve,corner=None):
        if corner is not None and self.original_corner is not corner:
            raise ValueError('Prepared circle piece changed its original Corner identity')
        ids=np.asarray(self.vertex_indices);points=np.asarray(self.reference_points)
        if (ids.dtype.kind not in 'iu' or np.any(ids<0) or np.any(ids>=len(curve.points))
                or not np.array_equal(curve.points[ids],points)
                or not np.array_equal(curve.stations[ids]/curve.length,self.prepared_parameters)):
            raise ValueError('Prepared circle piece changed its original source knots or vertices')
        c=self.original_corner;phase=np.asarray(self.circle_parameters)
        if c.center is None or c.radius<=0.:
            raise ValueError('Prepared circle piece requires its original actual circle')
        angle=c.start_angle+c.sweep*phase
        actual=c.center+c.radius*np.column_stack((np.cos(angle),np.sin(angle)))
        if np.max(np.linalg.norm(actual-points,axis=1))>1e-8:
            raise ValueError('Declared source-circle knots do not evaluate on the original analytical arc')
        for lead in self.auxiliary_leads:lead.validate(curve,c)
        return True


def declare(corner,curve,arc_phases,lead_tags,maximum_chord_error):
    """Only the producer calls this with tags retained during actual append."""
    ids=np.flatnonzero(np.isfinite(arc_phases))
    if len(ids)<2:return None
    if np.any(np.diff(ids)!=1):raise ValueError('One prepared circle piece has disconnected producer tags')
    leads=[]
    for role in (0,1):
        selected=np.array([i for i,tag in enumerate(lead_tags) if tag is not None and tag[0]==role],dtype=int)
        if len(selected)<2:continue
        if np.any(np.diff(selected)!=1):raise ValueError('Original auxiliary lead has disconnected producer tags')
        leads.append(PreparedAuxiliaryLead(corner.first if role==0 else corner.second,-1 if role==0 else 1,
            tuple(map(float,curve.stations[selected]/curve.length)),
            tuple(float(lead_tags[i][1]) for i in selected),tuple(map(int,selected)),
            tuple(tuple(map(float,point)) for point in curve.points[selected])))
    record=PreparedCircularInterval(corner,tuple(map(float,curve.stations[ids]/curve.length)),
        tuple(float(arc_phases[i]) for i in ids),tuple(map(int,ids)),
        tuple(tuple(map(float,point)) for point in curve.points[ids]),tuple(leads),float(maximum_chord_error))
    record.validate(curve,corner)
    return record
