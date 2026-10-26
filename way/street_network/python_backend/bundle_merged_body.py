"""Direct quad pavement and separator bodies with a merged entrance."""
import numpy as np
from types import SimpleNamespace
from .bundle_angular_island_chain import MultiIslandSurfaceBuilder
from .bundle_angular_open_islands import OpenAngularIsland, gap_intervals, clip_ray
from .bundle_merged_islands import OpenRunIsland
from .mesh import Assembler, Mesh
from .model import Key
from .validation import validate_overlaps


class MergedRunBody:
    max_cells = 8000

    def build(self, context, bundle, approach, merged_cut, *, remote_inset=8., remote_station=None, remote_cuts=None):
        field, frame, settings = approach.field, approach.frame_chart, context.network.settings
        # R1 retains the compact default budget. A requested finer mesh has a
        # bounded larger allowance; physical/error/shape limits are unchanged.
        self.max_cells=8000*settings.resolution
        local = frame.local(merged_cut[:, :2])
        merged_station = float(local[:, 0].mean())
        if np.ptp(local[:, 0]) > 1e-8:
            raise ValueError('Merged run needs an actual common-fibre body cut')
        if remote_station is None:
            remote_station = (min(ray.stations[-1] for ray in field.curves)-remote_inset
                if approach.sign > 0. else max(ray.stations[0] for ray in field.curves)+remote_inset)
        first, last = sorted((merged_station, float(remote_station)))
        merged_end = 0 if approach.sign > 0. else 1
        if last <= first+2.*settings.spacing:
            raise ValueError('The declared source run has no positive body interval')
        domain = SimpleNamespace(context=context, frame=frame, field=field, rays=field.curves,
            bundle=bundle, ports={i:SimpleNamespace(path=path) for i,path in enumerate(approach.paths)})
        regions, roots = gap_intervals(domain, first, last)
        if regions[0 if merged_end == 0 else -1][2]:
            raise ValueError('The proposed merged cut contains a real separator')
        builder = MultiIslandSurfaceBuilder(context.patterns, spacing=settings.spacing,
            curve_tolerance=settings.curve_tolerance, max_cells=self.max_cells)
        cuts, cut_keys = {}, {}
        order = np.argsort(-local[:, 1])
        ordered = merged_cut[order].copy()
        peak = int(np.argmax(ordered[:, 2]))
        if peak in (0,len(ordered)-1):
            raise ValueError('A merged cut lost its prescribed common crest sample')
        complete_keys = [('merged cut', int(index)) for index in range(len(merged_cut))]
        groups = (np.arange(peak+1), np.arange(len(ordered)-1, peak-1, -1))
        for role, indices in enumerate(groups):
            cuts[role, merged_end] = ordered[indices].copy()
            cut_keys[role, merged_end] = [complete_keys[int(order[index])] for index in indices]
        remote_end = 1-merged_end
        remote = first if remote_end == 0 else last
        ao, ai, bi, bo = field.edges(remote)
        high, low = field.peaks(remote)
        if ai <= bi:
            raise ValueError('Remote body cut needs its merged-end policy too')
        for role, (outer, inner) in enumerate(((ao, high), (bo, low))):
            if remote_cuts is not None:
                points = np.asarray(remote_cuts[role]).copy()
                if np.max(abs(frame.local(points[:,:2])[:,0]-remote)) > 1e-8:
                    raise ValueError('An actual external connector has a different common-fibre cut')
                if np.max(abs(points[:,2]-field.values(points[:,:2]))) > 1e-8:
                    raise ValueError('External whole profile disagrees with the original Bundle field')
                cuts[role, remote_end] = points
                cut_keys[role, remote_end] = [('remote cut',role,i) for i in range(len(points))]
                continue
            count = len(cuts[role, merged_end])-1
            y = np.linspace(outer, inner, count+1)
            xy = frame.world(np.c_[np.full(len(y),remote),y])
            cuts[role, remote_end] = np.c_[xy,field.values(xy)]
            cut_keys[role, remote_end] = [('remote cut',role,i) for i in range(len(y))]

        grass, rims, islands, open_cuts = Assembler(), {}, [], {}
        for index,(a,b,positive) in enumerate(regions):
            if not positive:
                continue
            failures = []
            for distance in (.5,5./7.,1.,9./7.):
                try:
                    island = OpenAngularIsland(domain)
                    mesh, boundary, report = island.build(a,b,roots,index,distance=distance)
                    break
                except ValueError as error:
                    failures.append(str(error))
            else:
                # Preserve every established successful distance first. The
                # local-pattern family is a fallback for the whole old family,
                # not a reason to redistribute an already accepted island.
                for distance in (.5,5./7.,1.,9./7.):
                    try:
                        island=OpenRunIsland(domain)
                        mesh,boundary,report=island.tiled(a,b,roots,index,distance=distance,shear=.2)
                        report['selection']='local tiles after all established open-island choices'
                        break
                    except ValueError as error:
                        failures.append(str(error))
                else:
                    raise ValueError('Merged-run island construction failed: '+str(failures))
            report['prior_failures'] = failures
            grass.add(mesh.xyz,mesh.quads,island.keys)
            rims[index] = boundary
            islands.append(report)
            for end,at in ((0,first),(1,last)):
                if (a if end == 0 else b) != at:
                    continue
                open_cuts[end] = island.cuts[end]
                for role in (0,1):
                    point = boundary[role][0][0 if end == 0 else -1]
                    if np.linalg.norm(cuts[role,end][-1,:2]-point) > 1e-8:
                        raise ValueError('Separator and pavement disagree at the original body cut')
                    cuts[role,end][-1,:2] = point
                    cuts[role,end][-1,2] = field.ceiling
                    cut_keys[role,end][-1] = boundary[role][1][0 if end == 0 else -1]

        def endpoint(station, role, which):
            if station == first or station == last:
                end = 0 if station == first else 1
                return cuts[role,end][-1,:2], cut_keys[role,end][-1]
            adjacent = [index for index,(a,b,positive) in enumerate(regions)
                        if positive and (b == station if which == 0 else a == station)]
            if len(adjacent) != 1:
                raise ValueError('A merged-run crest needs its adjacent separator event')
            index = adjacent[0]
            point, keys = rims[index][role]
            at = -1 if which == 0 else 0
            return point[at], keys[at]

        chains = [[],[]]
        for index,(a,b,positive) in enumerate(regions):
            if positive:
                for role in (0,1):
                    chains[role].append(rims[index][role])
                continue
            start, sk = endpoint(a,0,0)
            finish, fk = endpoint(b,0,1)
            ridge = builder.ridge(frame,field,start,finish,('merged-run crest',index))
            ridge[1][0],ridge[1][-1] = sk,fk
            for role in (0,1):
                own_start, own_sk = endpoint(a,role,0)
                own_finish, own_fk = endpoint(b,role,1)
                if (np.linalg.norm(own_start-start) > 1e-8 or np.linalg.norm(own_finish-finish) > 1e-8
                        or own_sk != sk or own_fk != fk):
                    raise ValueError('Overlapping pavement halves lost their one shared crest')
                chains[role].append(ridge)
        seams = [builder.join_chains(parts) for parts in chains]
        pavement, reports = Assembler(), []
        prescribed = {key:point for token,points in cuts.items() for key,point in zip(cut_keys[token],points)}
        for role in (0,1):
            ray = field.curves[0 if role == 0 else 3]
            curve = clip_ray(ray,first,last,{})
            curve,kinds,names = builder.simplify_sample(curve,[('outer',role)]*(len(curve)-1),
                1e-9 if any(path.exact for path in approach.paths) else .6*settings.curve_tolerance,
                ('merged-run outer',role))
            head,tail = cuts[role,0],cuts[role,1]
            hk,tk = cut_keys[role,0],cut_keys[role,1]
            if max(np.linalg.norm(head[0,:2]-curve[0]),np.linalg.norm(tail[0,:2]-curve[-1])) > 1e-8:
                raise ValueError('Merged-run source curb lost a prescribed cross-section endpoint')
            curve[[0,-1]] = head[0,:2],tail[0,:2]
            names[0],names[-1] = hk[0],tk[0]
            boundary = (np.vstack((head[::-1,:2],curve[1:],tail[1:,:2])),
                [('cut','head')]*(len(head)-1)+kinds+[('cut','tail')]*(len(tail)-1),
                hk[::-1]+names[1:]+tk[1:])
            name = f'merged run {role}/body'
            boundary = builder.balanced_physical(boundary,seams[role],name)
            mesh,keys,report = builder._direct_single(frame,field,boundary,seams[role],name)
            if not set(hk+tk) <= set(keys):
                raise ValueError('Merged run dropped a mandatory body-cut sample')
            for i,key in enumerate(keys):
                if key in prescribed:
                    if np.linalg.norm(mesh.xyz[i]-prescribed[key]) > 1e-8:
                        raise ValueError('Merged run changed a prescribed complete profile')
                    mesh.xyz[i] = prescribed[key]
            pavement.add(mesh.xyz,mesh.quads,keys)
            reports.append(report)
        pave = pavement.mesh()
        grassy = grass.mesh('grass') if grass.quads else None
        if len(pave.quads)+(len(grassy.quads) if grassy is not None else 0) > self.max_cells:
            raise ValueError('Merged-run pavement and grass exceed the combined quad budget')
        pe = builder.boundary_edges(pave)
        available = {tuple(sorted((tuple(pave.xyz[a]),tuple(pave.xyz[b])))) for a,b in pe}
        reserved = {tuple(sorted((tuple(np.r_[a,field.ceiling]),tuple(np.r_[b,field.ceiling]))))
                    for points,keys in open_cuts.values() for a,b in zip(points,points[1:])}
        interfaces = 0
        if grassy is not None:
            for a,b in builder.boundary_edges(grassy):
                edge = tuple(sorted((tuple(grassy.xyz[a]),tuple(grassy.xyz[b]))))
                if edge in reserved:
                    continue
                if edge not in available:
                    raise ValueError('Merged-run grass has no identical whole pavement edge')
                interfaces += 1
        published = Mesh(pave.xyz.astype(np.float32),pave.quads)
        quality = published.validate()
        if quality['maximum_edge_ratio'] > 8.:
            raise ValueError('Merged-run publication exceeds edge ratio eight')
        slopes = published.slope_statistics()
        longitudinal = frame.longitudinal_bound(published)
        if slopes['maximum'] > .1+1e-6 or longitudinal > .03+1e-8:
            raise ValueError('Merged-run body exceeds its physical slope limits')
        payloads = {Key('body_probe',0):pave}
        if grassy is not None:
            payloads[Key('grass_probe',0)] = grassy
        overlap = validate_overlaps(payloads,cell_size=2.)
        self.pavement_keys = [None]*len(pave.xyz)
        for key,at in pavement.lookup.items():
            self.pavement_keys[at] = key
        self.cuts,self.cut_keys,self.open_grass_cuts = cuts,cut_keys,open_cuts
        self.grass_keys = [None]*len(grassy.xyz) if grassy is not None else []
        for key,at in grass.lookup.items():
            self.grass_keys[at] = key
        return pave,grassy,dict(pavement=pave.validate(),published_quality=quality,
            grass=grassy.validate() if grassy is not None else None, islands=islands,
            body=reports,regions=regions,shared_grass_edges=interfaces,
            physical_slopes=slopes,angular_longitudinal_bound=longitudinal,
            overlap_check=overlap,merged_end=merged_end,first=first,last=last,
            external_rounded_cap_included=False,normal_addon_integration=True,all_direct_quads=True)
