"""Open-ended separator tiles using the maintained finite-island pattern family."""
import numpy as np
from .bundle_angular_open_islands import OpenAngularIsland,clip_ray
from .bundle_angular_curves import RayCurb
from .bundle_angular_islands import AngularIslandBuilder
from .bundle_angular_island_tiles import AngularIslandTiles
from .mesh import Assembler
from .validation import validate_overlaps


class OpenRunIsland(OpenAngularIsland):
    def build(self,first,last,roots,index,*,distance=.5,shear=.2):
        try:
            return super().build(first,last,roots,index,distance=distance,shear=shear)
        except ValueError as error:
            old_failure=str(error)
        result=self.tiled(first,last,roots,index,distance=distance,shear=shear)
        result[2]['single_atlas_failure']=old_failure
        return result

    def tiled(self,first,last,roots,index,*,distance,shear):
        domain=self.domain;context=domain.context;frame=domain.frame;settings=context.network.settings
        closed=(first in roots,last in roots)
        start=first+(distance if closed[0] else 0.)
        finish=last-(distance if closed[1] else 0.)
        if finish<=start+settings.spacing:
            raise ValueError('Tiled open-island tips consume the positive body')
        source=AngularIslandBuilder(frame,context.patterns,domain.bundle.key,spacing=settings.spacing,
            curve_tolerance=settings.curve_tolerance,height=domain.field.ceiling,max_cells=5000,
            boundary_error=(0. if all(port.path.exact for port in domain.ports.values()) else .1*settings.curve_tolerance)+1e-7)
        tiles=AngularIslandTiles(source)
        rays=[RayCurb(frame,clip_ray(ray,first,last,roots)) for ray in domain.rays[1:3]]
        error=max((source._cap_accuracy(ray,*interval) for ray in rays
            for interval in ((first,start),(finish,last)) if interval[1]>interval[0]),default=0.)
        if error>tiles.remaining:
            raise ValueError('Tiled open island exceeds its physical tip-error budget')
        smooth=all(not port.path.exact for port in domain.ports.values())
        stations=tiles._stations(rays,start,finish,smooth)
        pieces=tiles._pieces(rays,stations,smooth)
        cuts=tiles._cuts(rays,stations,pieces,distance,shear)
        for end,at in ((0,0),(1,-1)):
            if closed[end]:continue
            high,low=(ray(stations[at]) for ray in rays)
            import math
            count=max(4,4*int(math.ceil(np.linalg.norm(high-low)/(4.*settings.spacing))))
            cuts[at]=np.linspace(low,high,count+1)
        offsets=[np.r_[0,np.cumsum([len(piece[role])-1 for piece in pieces])] for role in (0,1)]
        def label(role,at):return ('open island',index,'upper' if role==0 else 'lower',int(at))
        cut_keys=[]
        for i,cut in enumerate(cuts):
            keys=[('open island',index,'local cut',i,j) for j in range(len(cut))]
            keys[0],keys[-1]=label(1,offsets[1][i]),label(0,offsets[0][i])
            cut_keys.append(keys)
        assembly,details=Assembler(),[]
        for i,(top,bottom) in enumerate(pieces):
            top_keys=[label(0,offsets[0][i]+j) for j in range(len(top))]
            bottom_keys=[label(1,offsets[1][i]+j) for j in range(len(bottom))]
            rim=np.vstack((bottom,cuts[i+1][1:-1],top[::-1],cuts[i][-2:0:-1]))
            keys=bottom_keys+cut_keys[i+1][1:-1]+top_keys[::-1]+cut_keys[i][-2:0:-1]
            right=len(cuts[i+1])-1
            corners=(0,len(bottom)-1,len(bottom)+right-1,len(bottom)+right+len(top)-2)
            mesh,report=tiles._tile(rim,corners,source.max_cells-len(assembly.quads)-sum(closed))
            keys.extend(('open island',index,'tile',i,j) for j in range(len(rim),len(mesh.xyz)))
            assembly.add(mesh.xyz,mesh.quads,keys);details.append(report)
        top,bottom=(np.concatenate([p[role] if i==0 else p[role][1:] for i,p in enumerate(pieces)]) for role in (0,1))
        exported=[(top,[label(0,i) for i in range(len(top))]),(bottom,[label(1,i) for i in range(len(bottom))])]
        for end,at in ((0,0),(1,-1)):
            if not closed[end]:continue
            root=roots[first if end==0 else last]
            key=('open island',index,'tip',end)
            if end==0:
                xy=np.vstack((root,bottom[0],cuts[0][1],top[0]))
                keys=[key,label(1,0),cut_keys[0][1],label(0,0)]
            else:
                xy=np.vstack((root,top[-1],cuts[-1][1],bottom[-1]))
                keys=[key,label(0,offsets[0][-1]),cut_keys[-1][1],label(1,offsets[1][-1])]
            assembly.add(np.c_[xy,np.full(4,domain.field.ceiling)],np.array(((0,1,2,3),)),keys)
            for role,(points,names) in enumerate(exported):
                exported[role]=(np.vstack((root,points)) if end==0 else np.vstack((points,root)),
                                [key]+names if end==0 else names+[key])
        mesh=assembly.mesh('grass');quality=mesh.validate()
        if quality['maximum_edge_ratio']>8.:
            raise ValueError('Tiled open separator exceeds the world edge ratio limit')
        validate_overlaps({domain.bundle.key:mesh},cell_size=2.)
        self.keys=[None]*len(mesh.xyz)
        for key,at in assembly.lookup.items():self.keys[at]=key
        self.cuts=((cuts[0],cut_keys[0]),(cuts[-1],cut_keys[-1]))
        return mesh,exported,dict(construction='local open-island direct-quad tiles',quality=quality,
            closed_tips=closed,tip_error=error,tiles=details,all_direct_quads=True,
            flat_cut_intervals=[len(cuts[0 if end==0 else -1])-1 if not closed[end] else None for end in (0,1)])
