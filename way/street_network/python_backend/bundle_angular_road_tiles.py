"""Local pavement charts retaining the exact existing grass/ridge interface."""
import math
import numpy as np
from .bundle_angular_curves import RayCurb
from .mesh import Assembler


class CutTooDense(ValueError):
    """The proposed computational cut would exceed the profile sample budget."""


class AngularRoadTiles:
    def __init__(self, builder, frame, field):
        self.builder,self.frame,self.field = builder,frame,field

    def _indices(self, ray, shared, stations):
        indices = [0]
        spacing = self.builder.spacing
        while indices[-1]<len(shared)-1:
            at = indices[-1]
            target = stations[at]+4.*spacing
            nearest = int(np.searchsorted(stations,target))
            if nearest>=len(shared)-2:
                indices.append(len(shared)-1)
                break
            choices = []
            for index in range(max(at+1,nearest-4),min(len(shared)-1,nearest+5)):
                station = stations[index]
                point = ray(station)
                lengths = np.linalg.norm(point-ray.points,axis=1)
                lengths = lengths[lengths>1e-8]
                clearance = min(float(lengths.min()),float(np.linalg.norm(shared[index]-shared[index-1])),
                                float(np.linalg.norm(shared[index+1]-shared[index])))
                choices.append((abs(station-target)/spacing+.05*spacing/max(clearance,1e-8),index))
            if not choices:
                raise ValueError('No bounded local road crown partition')
            indices.append(min(choices)[1])
            if len(indices)>512:
                raise ValueError('Local road body exceeds its tile budget')
        if len(indices)<3:
            raise ValueError('The road body has no useful local split')
        return indices

    def _fit_cuts(self, ray, shared, stations, indices):
        """Locate internal tile cuts on existing crown samples with enough clearance."""
        def fits(i):
            curb = ray(stations[i])
            # bounded_subchain reuses a source vertex within this station tolerance.
            separate = np.abs(ray.stations-stations[i]) > 1e-8
            distances = np.linalg.norm(ray.points[separate]-curb, axis=1)
            clearance = min(float(distances.min(initial=np.inf)),
                            float(np.linalg.norm(shared[i]-shared[i-1])),
                            float(np.linalg.norm(shared[i+1]-shared[i])))
            return np.linalg.norm(shared[i]-curb) <= 64*min(self.builder.spacing, 6.*clearance)

        result = [indices[0]]
        for index in indices[1:-1]:
            if index > result[-1] and fits(index):
                result.append(index)
                continue
            nearby = range(max(result[-1]+1, index-6), min(index+7, len(shared)-1))
            compatible = (i for i in nearby if fits(i))
            nearest = min(compatible, key=lambda i: abs(stations[i]-stations[index]), default=None)
            if nearest is not None:
                result.append(nearest)
        result.append(indices[-1])
        if len(result)<3:
            raise ValueError('The road body has no feasible local split')
        return result

    def build(self, outer, seam, name, *, fit_cuts=False):
        physical,kinds,keys = outer
        shared,shared_keys = seam
        curb_edges = [i for i,kind in enumerate(kinds) if kind[0]=='outer']
        first,last = min(curb_edges),max(curb_edges)+1
        ray = RayCurb(self.frame,physical[first:last+1])
        stations = self.frame.local(shared)[:,0]
        if np.any(np.diff(stations)<=1e-9):
            raise ValueError('Local road crown is not angularly monotone')
        indices = self._indices(ray,shared,stations)
        if fit_cuts:
            indices = self._fit_cuts(ray, shared, stations, indices)
        curbs = []
        for tile,(a,b) in enumerate(zip(indices,indices[1:])):
            start_key = keys[first] if tile==0 else (name,'local curb cut',tile)
            end_key = keys[last] if tile==len(indices)-2 else (name,'local curb cut',tile+1)
            curbs.append(self.builder.bounded_subchain(ray,keys[first:last+1],stations[a],stations[b],
                                                       start_key,end_key))
        cuts = [(physical[:first+1][::-1],list(keys[:first+1][::-1]))]
        for tile,index in enumerate(indices[1:-1],1):
            curb,crown = curbs[tile][0][0],shared[index]
            minimum = min(float(np.linalg.norm(curbs[tile-1][0][-1]-curbs[tile-1][0][-2])),
                          float(np.linalg.norm(curbs[tile][0][1]-curbs[tile][0][0])),
                          float(np.linalg.norm(shared[index]-shared[index-1])),
                          float(np.linalg.norm(shared[index+1]-shared[index])))
            count = max(2,int(math.ceil(np.linalg.norm(crown-curb)/min(self.builder.spacing,6.*minimum))))
            count += count%2
            if count>64:
                raise CutTooDense('Local road cut exceeds its bounded profile family')
            points = np.linspace(curb,crown,count+1)
            addresses = [(name,'local profile',tile,j) for j in range(count+1)]
            addresses[0],addresses[-1] = curbs[tile][1][0],shared_keys[index]
            cuts.append((points,addresses))
        cuts.append((physical[last:],list(keys[last:])))
        assembly,details = Assembler(),[]
        for tile,(a,b) in enumerate(zip(indices,indices[1:])):
            xy,labels = curbs[tile]
            start,end = cuts[tile],cuts[tile+1]
            boundary = np.vstack((start[0][::-1],xy[1:],end[0][1:]))
            addresses = start[1][::-1]+labels[1:]+end[1][1:]
            roles = ([('cut','local head')]*(len(start[0])-1)+[kinds[first]]*(len(xy)-1)
                     +[('cut','local tail')]*(len(end[0])-1))
            local_seam = shared[a:b+1],shared_keys[a:b+1]
            local_name = name+'/local '+str(tile)
            boundary = self.builder.balanced_physical((boundary,roles,addresses),local_seam,local_name)
            mesh,addresses,report = self.builder._mapped_single(self.frame,self.field,boundary,local_seam,local_name)
            if len(assembly.quads)+len(mesh.quads)>self.builder.max_cells:
                raise ValueError('Local road body exceeds its complete cell budget')
            assembly.add(mesh.xyz,mesh.quads,addresses)
            details.append(report)
        result = assembly.mesh()
        quality = result.validate()
        if quality['maximum_edge_ratio']>8.:
            raise ValueError('Local pavement world ratio exceeds 8')
        addresses = [None]*len(result.xyz)
        for key,index in assembly.lookup.items():
            addresses[index] = key
        return result,addresses,dict(construction='local pavement body tiles with existing crown cuts',
            world_quality=quality,patches=details,shared_crown_vertices_unchanged=True,local_tiles=len(details))
