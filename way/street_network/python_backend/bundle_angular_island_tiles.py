"""Local physical quad tiles for a finite angular separator island.

Original curb chains and their error budget remain authoritative. Crosscut
partitions are chosen together for adjacent tiles, so a narrow boundary cell
does not force a dense partition through the whole island. No triangles or
coordinate-based welding are used.
"""
import copy
import math
import numpy as np
from .mesh import Mesh, Assembler
from .bundle_axis import _simplify, _regular_rows
from .bundle_angular_curves import RayCurb
from .four_side_charts import FourSidedQuadChart, quality
from .quad_charts import FixedRimQuality
from .validation import validate_overlaps


class AngularIslandTiles:
    def __init__(self, source):
        self.source = source
        self.sampler = copy.copy(source)
        self.simplification = .1*source.curve_tolerance
        self.remaining = source.curve_tolerance-source.boundary_error-2.*self.simplification-1e-8
        if self.remaining<=0.:
            raise ValueError('Local island tiles leave no physical curve-error budget')
        self.sampler.body_tolerance = min(source.body_tolerance, self.remaining)

    def _trim(self, ray, a, b, smooth):
        return self.sampler._trimmed(ray, a, b, smooth, self.simplification)

    def _stations(self, rays, first, last, smooth):
        spacing = self.source.spacing
        knots = np.unique(np.concatenate([ray.stations for ray in rays]))
        middle = .5*(knots[:-1]+knots[1:])
        result = [first]
        while result[-1]<last-1e-10:
            width = np.linalg.norm(rays[0](result[-1])-rays[1](result[-1]))
            step = max(.25*spacing, min(2.*spacing, 2.5*width))
            target = min(last, result[-1]+step)
            if target<last-1e-10:
                candidates = middle[(middle>result[-1]+.5*step)&(middle<result[-1]+1.5*step)
                                    &(middle<last-.25*spacing)]
                if len(candidates):
                    index = np.clip(np.searchsorted(knots,candidates),1,len(knots)-1)
                    clearance = np.minimum(abs(candidates-knots[index-1]),abs(knots[index]-candidates))
                    score = abs(candidates-target)/step+.2*step/np.maximum(clearance,1e-8)
                    target = float(candidates[np.argmin(score)])
            result.append(target)
            if len(result)>4096:
                raise ValueError('Local island tile count exceeds its allocation budget')
        if len(result)>2 and result[-1]-result[-2]<.5*(result[-2]-result[-3]):
            result.pop(-2)
        for _ in range(2):
            for at in range(1,len(result)-1):
                choices = []
                base = result[at]
                for shift in (0.,-.05,.05,-.1,.1,-.2,.2,-.3,.3):
                    value = base+shift*spacing
                    if not result[at-1]+.1*spacing<value<result[at+1]-.1*spacing:
                        continue
                    lengths = []
                    for ray in rays:
                        previous = self._trim(ray,result[at-1],value,smooth)
                        following = self._trim(ray,value,result[at+1],smooth)
                        lengths.extend((np.linalg.norm(previous[-1]-previous[-2]),
                                        np.linalg.norm(following[1]-following[0])))
                    choices.append((min(lengths)/(1.+abs(shift)),value))
                if choices:
                    result[at] = max(choices)[1]
        return result

    def _pieces(self, rays, stations, smooth):
        pieces,total = [],0
        for a,b in zip(stations,stations[1:]):
            top,bottom = (self._trim(ray,a,b,smooth) for ray in rays)
            width = min(np.linalg.norm(rays[0](at)-rays[1](at)) for at in (a,b))
            spacing = min(self.source.spacing,2.5*width)
            if spacing<1e-6:
                raise ValueError('Local source strip approaches a zero-width contact')
            top,bottom = (_regular_rows(chain,spacing,multiple=1) for chain in (top,bottom))
            if len(top)==2:
                top = np.vstack((top[0],.5*(top[0]+top[1]),top[1]))
            if len(bottom)==2:
                bottom = np.vstack((bottom[0],.5*(bottom[0]+bottom[1]),bottom[1]))
            if (len(top)-len(bottom))%2:
                at = int(np.argmax(np.linalg.norm(np.diff(top,axis=0),axis=1)))
                top = np.insert(top,at+1,.5*(top[at]+top[at+1]),axis=0)
            total += len(top)+len(bottom)-2
            if total>2*self.source.max_cells+2:
                raise ValueError('Local island boundary exceeds its complete cell budget')
            pieces.append((top,bottom))
        return pieces

    def _cuts(self, rays, stations, pieces, distance, shear):
        cuts = []
        for at,station in enumerate(stations):
            high,low = (ray(station) for ray in rays)
            if at in (0,len(stations)-1):
                middle_station = station+distance*shear*(1 if at==0 else -1)
                middle = .5*(rays[0](middle_station)+rays[1](middle_station))
                cuts.append(np.vstack((low,middle,high)))
                continue
            minimum = min(np.linalg.norm(chain[-1]-chain[-2]) for chain in pieces[at-1])
            minimum = min(minimum,*(np.linalg.norm(chain[1]-chain[0]) for chain in pieces[at]))
            width = float(np.linalg.norm(high-low))
            count = max(2,int(math.ceil(width/min(self.source.spacing,6.*minimum))))
            count += count%2
            if count>64:
                raise ValueError('A local island crosscut exceeds its bounded profile family')
            cuts.append(low+np.linspace(0.,1.,count+1)[:,None]*(high-low))
        return cuts

    @staticmethod
    def _harmonic(mesh, rim_count):
        if len(mesh.xyz)>512:
            raise ValueError('Local harmonic initialization exceeds its vertex budget')
        edges = np.unique(np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1)
                                  .reshape(-1,2),axis=1),axis=0)
        matrix = np.zeros((len(mesh.xyz),len(mesh.xyz)))
        for a,b in edges:
            matrix[a,a] += 1.; matrix[b,b] += 1.
            matrix[a,b] -= 1.; matrix[b,a] -= 1.
        mesh.xyz[rim_count:,:2] = np.linalg.solve(matrix[rim_count:,rim_count:],
                                                -matrix[rim_count:,:rim_count]@mesh.xyz[:rim_count,:2])
        return mesh

    def _tile(self, rim, corners, budget):
        chart = FourSidedQuadChart(self.source.patterns,maximum_ratio=8.,max_cells=budget)
        origin = rim[0]
        try:
            mapped,counts = chart.map(rim-origin,corners)
            try:
                initial = quality(mapped)
            except ValueError:
                mapped = self._harmonic(mapped,len(rim))
                initial = quality(mapped)
            if not chart.accepted(initial):
                mapped,_ = FixedRimQuality(iterations=40,power=2).improve(mapped)
                if not chart.accepted(quality(mapped)):
                    raise ValueError('The local semantic atlas fails its shape limits')
        except (ValueError,np.linalg.LinAlgError):
            alternative = chart.build(rim-origin)
            mapped,counts = alternative.mesh,alternative.metadata['rim_intervals']
        xyz = np.c_[mapped.xyz[:,:2]+origin,np.full(len(mapped.xyz),self.source.height)]
        xyz[:len(rim),:2] = rim
        mesh = Mesh(xyz,mapped.quads,'grass')
        report = mesh.validate()
        if report['maximum_edge_ratio']>8.:
            raise ValueError('Local island tile exceeds the world quad ratio limit')
        return mesh,dict(quality=report,counts=counts)

    def _candidate(self, upper, lower, distance, shear, smooth):
        source = self.source
        rays = [RayCurb(source.frame,_simplify(p,self.simplification)) for p in (upper,lower)]
        start,finish = max(ray.stations[0] for ray in rays),min(ray.stations[-1] for ray in rays)
        first,last = start+distance,finish-distance
        if last-first<source.spacing:
            raise ValueError('Local island tips consume the body')
        error = max(source._cap_accuracy(ray,a,b) for ray in rays
                    for a,b in ((start,first),(last,finish)))
        if error>self.remaining:
            raise ValueError('Local island tip exceeds its physical curve-error budget')
        stations = self._stations(rays,first,last,smooth)
        pieces = self._pieces(rays,stations,smooth)
        cuts = self._cuts(rays,stations,pieces,distance,shear)
        offsets = [np.r_[0,np.cumsum([len(piece[role])-1 for piece in pieces])] for role in (0,1)]
        labels = lambda role,at:('island','upper' if role==0 else 'lower',int(at))
        cut_keys = []
        for i,cut in enumerate(cuts):
            keys = [('island','crosscut',i,j) for j in range(len(cut))]
            keys[0],keys[-1] = labels(1,offsets[1][i]),labels(0,offsets[0][i])
            if i in (0,len(cuts)-1):
                keys[1] = ('island','cut_middle',int(i!=0))
            cut_keys.append(keys)
        assembly,details = Assembler(),[]
        for i,(top,bottom) in enumerate(pieces):
            top_keys = [labels(0,offsets[0][i]+j) for j in range(len(top))]
            bottom_keys = [labels(1,offsets[1][i]+j) for j in range(len(bottom))]
            right = len(cuts[i+1])-1
            rim = np.vstack((bottom,cuts[i+1][1:-1],top[::-1],cuts[i][-2:0:-1]))
            keys = bottom_keys+cut_keys[i+1][1:-1]+top_keys[::-1]+cut_keys[i][-2:0:-1]
            corners = (0,len(bottom)-1,len(bottom)+right-1,len(bottom)+right+len(top)-2)
            budget = source.max_cells-len(assembly.quads)-2
            if budget<1:
                raise ValueError('Local island exceeds its complete cell budget')
            mesh,report = self._tile(rim,corners,budget)
            keys.extend(('island','tile',i,j) for j in range(len(rim),len(mesh.xyz)))
            assembly.add(mesh.xyz,mesh.quads,keys)
            details.append(report)
        top,bottom = (np.concatenate([p[role] if i==0 else p[role][1:] for i,p in enumerate(pieces)])
                      for role in (0,1))
        roots = np.asarray((upper[0],upper[-1]))
        for side in (0,1):
            if side==0:
                xy = np.vstack((roots[0],bottom[0],cuts[0][1],top[0]))
                keys = [('island_tip',0),labels(1,0),cut_keys[0][1],labels(0,0)]
            else:
                xy = np.vstack((roots[1],top[-1],cuts[-1][1],bottom[-1]))
                keys = [('island_tip',1),labels(0,offsets[0][-1]),cut_keys[-1][1],labels(1,offsets[1][-1])]
            assembly.add(np.c_[xy,np.full(4,source.height)],np.array(((0,1,2,3),)),keys)
        mesh = assembly.mesh('grass')
        report = mesh.validate()
        if len(mesh.quads)>source.max_cells or report['maximum_edge_ratio']>8.:
            raise ValueError('Local island exceeds its complete cell or world ratio limit')
        validate_overlaps({source.owner:mesh},cell_size=2.)
        rims = [(np.vstack((roots[:1],points,roots[1:])),
                 [('island_tip',0)]+[labels(role,i) for i in range(len(points))]+[('island_tip',1)])
                for role,points in enumerate((top,bottom))]
        return mesh,rims,dict(construction='local physical direct-quad island tiles',quality=report,
            quads=len(mesh.quads),cap_distance=distance,shear=shear,blocks=details,
            boundary_chord_error=error,source_boundary_error=source.boundary_error,
            source_simplification=self.simplification,body_tolerance=self.sampler.body_tolerance,
            complete_curve_error_bound=source.boundary_error+2.*self.simplification+max(error,self.sampler.body_tolerance))

    def build(self, upper, lower, *, smooth=False):
        failures = []
        for factor in (.5,1.,1.5,2.,2.5,3.):
            for shear in (1./6.,.2,.25):
                try:
                    mesh,rims,report = self._candidate(upper,lower,factor*self.source.spacing,shear,smooth)
                    report['failures'] = failures
                    return mesh,rims,report
                except ValueError as error:
                    failures.append(dict(factor=factor,shear=shear,error=str(error)))
        raise ValueError('No qualified local island tiles: '+str(failures))
