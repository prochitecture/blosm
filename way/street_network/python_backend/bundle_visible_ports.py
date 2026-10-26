"""Original complete mouth samples and shared physical-owner cap resolution.

No generated vertex index becomes a source identity. Parameter maps are
ephemeral views over immutable cached quad patterns.
"""
from dataclasses import dataclass
import numpy as np
from .model import Key
from .ownership import PLAN, CompoundProduct, surface_key
from .separator_ports import PORTS
from .patterns import Grid
from .bundles import Cap


@dataclass(frozen=True)
class VisibleMouth:
    endpoint: object
    owner: object
    station: float
    parameters: np.ndarray
    physical_xyz: np.ndarray

    @staticmethod
    def orient(outward_u, physical_xyz, *, end, role, forward):
        """Order actual samples from outer to inner in the directed member."""
        u=np.asarray(outward_u,dtype=float)
        xyz=np.asarray(physical_xyz,dtype=float)
        if (end not in (0,1) or role not in (0,1) or u.ndim!=1
                or xyz.shape!=(len(u),3) or len(u)<2 or not np.all(np.isfinite(xyz))
                or not np.all(np.isfinite(u)) or u.min()<0. or u.max()>1.):
            raise ValueError('Invalid original compound mouth samples')
        original_q=(2.*u-1.)*(1 if end==0 else -1)
        v=.5*(1.+original_q*(1 if role==0 else -1)*(1 if forward else -1))
        order=np.argsort(v);v=v[order];xyz=xyz[order].copy()
        if v[0]!=0. or v[-1]!=1. or np.any(np.diff(v)<=0.):
            raise ValueError('A full visible profile needs unique outer/inner source samples')
        v.flags.writeable=xyz.flags.writeable=False
        return v,xyz

    @classmethod
    def from_compound(cls,context,endpoint,role,forward):
        port=context[PORTS].resolve(endpoint)
        if port is None:raise ValueError('No actual visible port for the original run endpoint')
        plan=context.products.get(PLAN)
        owner=(plan.owner(port.site.key) if plan is not None else None) or port.site.key
        product=context[surface_key(owner)]
        if not isinstance(product,CompoundProduct):return None
        identity=endpoint.section_id,endpoint.end
        if identity not in product.peaks:raise ValueError('Run endpoint is internal to a compound junction')
        station=product.peaks[identity][2];values=[]
        for (sid,end,u),(loop,node,at) in product.port_points.items():
            if (sid,end)!=identity:continue
            if abs(at-station)>1e-8:raise ValueError('A source profile has conflicting mouth stations')
            xyz=product.rims[loop][node].copy()
            xyz[2]=(product.field.datum+(xyz[2]-product.field.datum)/product.exaggeration
                    if product.exaggeration else float(product.field(xyz[None,:2])[0]))
            values.append((u,xyz))
        parameters,xyz=cls.orient([row[0] for row in values],[row[1] for row in values],
            end=endpoint.end,role=role,forward=forward)
        return cls(endpoint,owner,float(station),parameters,xyz)


class ProfileSamplingPatterns:
    """Ephemeral parameter maps; the shared immutable pattern cache is unchanged."""
    def __init__(self,patterns,parameters):
        values=np.asarray(parameters,dtype=float).copy()
        if values.ndim!=1 or len(values)<2 or values[0]!=0. or values[-1]!=1. or np.any(np.diff(values)<=0.):
            raise ValueError('Profile parameter map must retain an ordered complete cross-section')
        values.flags.writeable=False
        self.patterns,self.parameters=patterns,values
        self.uniform=np.linspace(0.,1.,len(values))
    def __getattr__(self,name):return getattr(self.patterns,name)
    def _map(self,pattern):
        if np.array_equal(self.parameters,self.uniform):return pattern
        uv=pattern.uv.copy();uv[:,1]=np.interp(uv[:,1],self.uniform,self.parameters)
        uv.flags.writeable=False
        return Grid(uv,pattern.quads,pattern.rows,pattern.columns)
    def grid(self,*args):return self._map(self.patterns.grid(*args))
    def quad(self,*args):return self._map(self.patterns.quad(*args))


class VisibleCaps:
    @staticmethod
    def common_compound(context,first,second,end,rim_coordinates=None):
        endpoints=first.ends[end],second.ends[end]
        ports=[context[PORTS].resolve(endpoint) for endpoint in endpoints]
        if not all(ports):return None
        plan=context.products.get(PLAN)
        owners=[(plan.owner(port.site.key) if plan is not None else None) or port.site.key for port in ports]
        if owners[0]!=owners[1]:return None
        product=context[surface_key(owners[0])]
        if not isinstance(product,CompoundProduct):return None
        q=(rim_coordinates if rim_coordinates is not None else
           tuple(context[Key('profile',endpoint.section_id)].peak for endpoint in endpoints))
        if any(value not in (-1.,1.) for value in q):
            raise ValueError('A separator cap needs explicit inner half-profile sides')
        outward=tuple(.5*(value*(1 if endpoint.end==0 else -1)+1.) for value,endpoint in zip(q,endpoints))
        return Cap.rounded(product.edge_cap(*((endpoint.section_id,endpoint.end) for endpoint in endpoints),*outward))
