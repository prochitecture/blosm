"""Grass between a wholly source-owned member and its true free counterpart.

This is an explicit separator family, not a pavement union or a source
preparation. It never repairs a clipped owner: full consumed/prefix width must
already pass on its actual physical product. Private free closures end only
at existing, identity-proven original curb vertices.
"""
from dataclasses import dataclass
import math
import numpy as np
from .bundle_angular_curves import polyline_distance
from .bundle_merged_body import MergedRunBody
from .bundle_merged_plan import physical_context
from .bundles import _road_rim, CenterlineRange, certify_disjoint
from .coincident_source_curbs import _intervals, curve_parameter
from .compounds import rim_label
from .drainage import BoundaryDrainage
from .four_side_charts import quality
from .mesh import Mesh
from .model import End, Key
from .ownership import PLAN, CompoundProduct, surface_key
from .providers import CarriagewayProvider, Part
from .regions import Region, RegionUnion, signed_area
from .routes import RouteTopology
from .separator_surfaces import FixedRimSeparator, SeparatorSurface
from .fixed_rim_bands import BandedSeparatorChart
from .profiles import plateau_extension
from .wholly_owned_bundle_preparation import WhollyOwnedBundleUnavailable
from .updates import BuildContext
from .validation import validate_overlaps


def _covered(metadata, length):
    rows = metadata.get('compound_spans', ())
    owners = {row['owner'] for row in rows}
    if len(owners) != 1 or metadata['end'] > metadata['start']+1e-7:
        return None
    last = 0.
    for row in sorted(rows, key=lambda item: item['start']):
        if row['start'] > last+1e-7: return None
        last = max(last, row['end'])
    return Key.parse(next(iter(owners))) if abs(last-length) <= 1e-7 else None


def _runs(mask):
    """Existing complete edge runs, with no numerical bridge across a gap."""
    if np.all(mask): return ()
    result = []; n = len(mask)
    for start in range(n):
        if not mask[start] or mask[(start-1) % n]: continue
        ids = [start]; at = start
        while mask[at]:
            at = (at+1) % n; ids.append(at)
            if len(ids) > n: raise ValueError('An original owned bank has no finite end')
        result.append(ids)
    return tuple(result)


def _original_edge(context, source, span, sources, path, offset, reference):
    label = rim_label(span, sources)
    if label.kind != 'curb': return False
    layout = sources[span.owner].layout
    curve = layout.boundaries[label.index]
    for row in _intervals(curve, label.first, label.last, .5*(label.first+label.last)):
        if row['path'] is not path or abs(row['offset']-offset) > 1e-8: continue
        bounds = sorted(curve_parameter(row, station) for station in (row['low'], row['high']))
        if bounds[0] <= min(label.first,label.last)+1e-9 and bounds[1] >= max(label.first,label.last)-1e-9:
            return True
    # A full source-strip operand records both original root incidences.
    # Its Path/width was already certified continuously; distance below only
    # chooses the declared physical side of that SAME original source.
    caps = tuple(getattr(layout, 'original_root_caps', {}).values())
    if ({end.section_id for end in caps} != {source.id} or {end.end for end in caps} != {0,1}):
        return False
    points = curve(np.array((label.first, .5*(label.first+label.last), label.last)))
    limit = .07*context.network.settings.curve_tolerance+1e-7
    return float(polyline_distance(points, reference).max()) <= limit


def _plane_anchor(points, origin, tangent, spacing, end):
    signed = (points[:, :2]-origin)@tangent
    hits = []
    for i, (a,b) in enumerate(zip(signed,signed[1:])):
        if abs(a)<=1e-10 and abs(b)<=1e-10:
            raise WhollyOwnedBundleUnavailable('A free original plane coincides with a positive owned-bank interval')
        if a*b > 0. or abs(b-a) < 1e-12: continue
        fraction = float(np.clip(a/(a-b),0.,1.))
        xy = (1.-fraction)*points[i,:2]+fraction*points[i+1,:2]
        if hits and np.linalg.norm(hits[-1][1]-xy) < 1e-9: continue
        if abs(a) <= 1e-10:index=i
        elif abs(b) <= 1e-10:index=i+1
        else:index=i+1 if end==0 else i
        edge_length = float(np.linalg.norm(points[i+1,:2]-points[i,:2]))
        if edge_length > spacing+1e-7:
            raise WhollyOwnedBundleUnavailable('A free closure needs a more closely sampled actual owned curb')
        hits.append((index,xy,dict(edge=i,fraction=fraction,
            exact_plane_intersection=xy.tolist(),existing_vertex=index,
            selection='first existing vertex inside the two true free planes',
            existing_vertex_plane_error=abs(float(signed[index])),whole_owned_edge_length=edge_length,
            added_owner_vertices=0),b>a))
    if len(hits) != 1 or not hits[0][3]:
        raise WhollyOwnedBundleUnavailable('A free original profile has no unique finite owned-bank crossing')
    return hits[0][:3]


def _physical_mask(rim, loop, eligible):
    """Map only roundoff-equal subdivisions of this same owned physical rim."""
    a, b = loop.xy, np.roll(loop.xy,-1,axis=0)
    edge = b-a; squared=np.sum(edge*edge,axis=1)
    rounding=512.*np.finfo(float).eps*max(1.,float(np.max(abs(rim[:,:2]))))
    low,high=np.minimum(a,b)-rounding,np.maximum(a,b)+rounding
    result=[]
    for first,last in zip(rim[:,:2],np.roll(rim[:,:2],-1,axis=0)):
        ids=np.flatnonzero(np.all(first>=low,axis=1)&np.all(first<=high,axis=1)
            &np.all(last>=low,axis=1)&np.all(last<=high,axis=1))
        valid=[]
        for i in ids:
            offsets=np.array((first,last))-a[i]
            fraction=offsets@edge[i]/squared[i]
            if (np.all(fraction>=-1e-10) and np.all(fraction<=1.+1e-10)
                    and np.max(np.linalg.norm(offsets-fraction[:,None]*edge[i],axis=1))<=rounding):
                valid.append(int(i))
        if len(valid)!=1:
            raise ValueError('Owned bank lost its exact whole physical-rim subdivision provenance')
        result.append(bool(eligible[valid[0]]))
    return np.asarray(result)


@dataclass
class WhollyOwnedBundleDomain:
    owner: object
    owned_role: int
    owned_section: int
    free_section: int
    owner_bank: np.ndarray
    free_bank: np.ndarray
    closures: tuple
    physical: object
    free_road: object
    proof: dict

    @classmethod
    def build(cls, context, bundle):
        if (len(bundle.members) != 2 or not all(m.complete and len(m.path)==1 for m in bundle.members)):
            raise WhollyOwnedBundleUnavailable('Owned/free separators need two complete single original Sections')
        if PLAN not in context.products:
            raise WhollyOwnedBundleUnavailable('An owned/free separator needs actual source ownership')
        topology = RouteTopology(context.network); attachments = topology.attachments(context.network)
        visits = tuple(member.path[0] for member in bundle.members)
        metadata = [CarriagewayProvider.span_metadata(context,context.network.section_index[v.section_id],
                    attachments[v.section_id]) for v in visits]
        owners = [_covered(row,context[Key('path',v.section_id)].length) for row,v in zip(metadata,visits)]
        if sum(owner is not None for owner in owners) != 1:
            raise WhollyOwnedBundleUnavailable('An owned/free separator needs exactly one wholly owned original member')
        role = next(i for i,owner in enumerate(owners) if owner is not None); free = 1-role; owner = owners[role]
        source, counterpart = (context.network.section_index[visits[i].section_id] for i in (role,free))
        if (counterpart.scope_ends or counterpart.selection_ends
                or any(end.section_id==counterpart.id for site in (*context.network.intersections,*context.network.connections)
                       for end in site.ends)
                or metadata[free]['start'] != 0.
                or metadata[free]['end'] != context[Key('path',counterpart.id)].length):
            raise WhollyOwnedBundleUnavailable('The counterpart must have two actual unscoped free original ends')
        incident = {end:site for site in context.network.intersections for end in site.ends}
        original_ends = (End(source.id,0),End(source.id,1))
        if any(end not in incident or (context[PLAN].owner(incident[end].key) or incident[end].key)!=owner
               for end in original_ends):
            raise WhollyOwnedBundleUnavailable('Whole member ownership needs both real original incidences')
        physical = physical_context(context,required_owners={incident[end].key for end in original_ends})
        product = physical[surface_key(owner)]
        if not isinstance(product,CompoundProduct) or source.id not in {c.section_id for c in product.group.claims}:
            raise WhollyOwnedBundleUnavailable('The complete original member needs its actual compound source claim')
        consumed = {claim.section_id for claim in product.group.claims}
        # Fatal if the accepted owner clips ANY consumed or incoming source.
        from .wholly_owned_bundle_preparation import certify as certify_free_counterpart
        certificate=certify_free_counterpart(physical,product.group,store=False)
        if not isinstance(certificate,dict) or not any(row['bundle']==bundle.id for row in certificate['candidates']):
            raise ValueError('The owned/free source has no complete independent counterpart certificate')
        width_proof=certificate['complete_consumed_and_prefix_width']
        path = physical[Key('path',source.id)]
        sign = (1. if role==0 else -1.)*(1. if visits[role].forward else -1.)
        offset = source.offset+sign*.5*source.width
        stations = path.samples(0.,path.length,min(.25,physical.network.settings.spacing),
            tolerance=.02*physical.network.settings.curve_tolerance,max_samples=20000)
        reference = path.point(stations,offset)
        original_profile = physical[Key('profile',source.id)].oriented(int(not visits[role].forward))
        free_profile = physical[Key('profile',counterpart.id)].oriented(int(not visits[free].forward))
        if (original_profile.half != (1 if role==0 else -1)
                or free_profile.half != (1 if free==0 else -1)
                or abs(original_profile.datum-free_profile.datum)+abs(original_profile.rise-free_profile.rise)>1e-9):
            raise WhollyOwnedBundleUnavailable('Owned/free banks need their original matching inward half crowns')
        free_path = physical[Key('path',counterpart.id)]
        separation = certify_disjoint((CenterlineRange.build(path,0.,path.length),),
            (CenterlineRange.build(free_path,0.,free_path.length),))
        road = CarriagewayProvider.build(physical,counterpart,attachments[counterpart.id])
        if road.mesh is None: raise ValueError('The true free counterpart has no normal physical carriageway')
        products = dict(physical.products);products[Key('road',counterpart.id)] = road
        sampled = BuildContext(physical.network,products,physical.patterns)
        side = (1. if free==0 else -1.)*(1. if visits[free].forward else -1.)
        free_bank = _road_rim(sampled,counterpart.id,side)
        if not visits[free].forward: free_bank = free_bank[::-1]
        endpoints = []
        for end in (0,1):
            at = (0. if end==0 else free_path.length) if visits[free].forward else (free_path.length if end==0 else 0.)
            point,tangent,_ = free_path.frame(at)
            endpoints.append((point,tangent*(1. if visits[free].forward else -1.)))
        sources = {f.region.owner:f for f in product.group.footprints};candidates=[]
        for li,loop in enumerate(product.group.loops):
            if loop.hole: continue
            rim = product.physical_rims[li]
            eligible = np.asarray([_original_edge(physical,source,span,sources,path,offset,reference) for span in loop.spans])
            mask = _physical_mask(rim,loop,eligible)
            for ids in _runs(mask):
                points = rim[ids].copy()
                if role==1: points=points[::-1]
                try:
                    first = _plane_anchor(points,*endpoints[0],physical.network.settings.spacing,0)
                    last = _plane_anchor(points,*endpoints[1],physical.network.settings.spacing,1)
                except WhollyOwnedBundleUnavailable: continue
                if first[0] >= last[0]: continue
                bank = points[first[0]:last[0]+1].copy()
                candidates.append((bank,(first[2],last[2]),li))
        if len(candidates) != 1:
            raise WhollyOwnedBundleUnavailable('The original owned member has no unique continuous finite inner bank')
        bank,closures,loop = candidates[0]
        return cls(owner,role,source.id,counterpart.id,bank,free_bank.copy(),closures,physical,road,
            dict(source_width=width_proof,source_identity='original Path/offset leads or full original root-cap strip',
                complete_free_ribbon=certificate,
                complete_original_centerline_separation=separation,
                original_member=source.id,owner_loop=loop,bank_vertices=len(bank),
                free_original_ends=True,scope_or_selection_connections=0,source_aliases_created=0))


class WhollyOwnedBundleSeparator:
    @classmethod
    def build(cls, context, bundle):
        domain = WhollyOwnedBundleDomain.build(context,bundle)
        physical = domain.physical;settings = physical.network.settings
        profile = physical[Key('profile',domain.free_section)]
        first,last = domain.owner_bank,domain.free_bank
        counts = [max(2,int(math.ceil(np.linalg.norm(first[i,:2]-last[i,:2])/settings.spacing))) for i in (0,-1)]
        total = len(first)-1+len(last)-1+sum(counts)
        if total%2: counts[0]+=1
        head = np.linspace(last[0],first[0],counts[0]+1)
        tail = np.linspace(first[-1],last[-1],counts[1]+1)
        boundary = np.vstack((first,tail[1:],last[-2::-1],head[1:-1]))
        if signed_area(boundary[:,:2]) < 0.:boundary=boundary[::-1].copy()
        loops = RegionUnion(max_segments=8192).build((Region(bundle.key,boundary[:,:2]),))
        if len(loops)!=1 or loops[0].hole: raise ValueError('Owned/free separator needs one actual simple finite domain')
        try:
            surface = FixedRimSeparator(context.patterns,spacing=settings.spacing,max_cells=50000).build(boundary,
                datum=profile.datum,height=profile.rise,exaggeration=1.)
        except ValueError as previous:
            failures=[str(previous)]
            for bands in (2,4,8,16):
                try:
                    chart=BandedSeparatorChart(context.patterns,spacing=settings.spacing,max_cells=50000,bands=bands).build(boundary[:,:2])
                    candidate=chart.mesh
                    candidate.xyz[:,2],height=plateau_extension(candidate.xyz[:,:2],boundary,profile.datum,profile.rise)
                    candidate.xyz[:len(boundary)]=boundary;candidate.material='grass'
                    stored_candidate=Mesh(candidate.xyz.astype(np.float32).astype(float),candidate.quads,'grass')
                    if max(candidate.validate()['maximum_edge_ratio'],stored_candidate.validate()['maximum_edge_ratio'])>8.:
                        raise ValueError('An owned/free band needs better complete stored grading')
                    surface=SeparatorSurface(candidate,dict(chart=chart.metadata,height=height,prior_failures=failures))
                    break
                except ValueError as error:failures.append(str(error))
            else:raise ValueError('No accepted owned/free fixed-rim grass chart: '+str(failures)) from previous
        mesh = surface.mesh
        stored = Mesh(mesh.xyz.astype(np.float32).astype(float),mesh.quads,'grass')
        checks = [quality(mesh),quality(stored)]
        if max(row['maximum_edge_ratio'] for row in checks)>8.:raise ValueError('Owned/free grass exceeds stored ratio eight')
        stored.validate()
        slope = max(mesh.slope_statistics()['maximum'],stored.slope_statistics()['maximum'])
        height_recovery=None
        if slope>.1+1e-6:
            from .fixed_grass_heights import reconcile
            mesh,_,_,height_recovery=reconcile(mesh,datum=profile.datum,ceiling=profile.datum+profile.rise,
                max_cells=50000,policy='owned_free_actual_whole_banks')
            stored=Mesh(mesh.xyz.astype(np.float32).astype(float),mesh.quads,'grass')
            checks=[quality(mesh),quality(stored)]
            slope=max(mesh.slope_statistics()['maximum'],stored.slope_statistics()['maximum'])
        if slope>.1+1e-6:raise ValueError('Owned/free grass exceeds its unchanged physical slope limit')
        fixed=np.zeros(len(mesh.xyz),bool);fixed[np.unique(MergedRunBody._boundary_vertices(mesh))]=True
        drainage=BoundaryDrainage(max_cells=50000).fill(stored,fixed)[1]
        if drainage['maximum_fill']>1e-6:raise ValueError('Owned/free grass creates an interior depression')
        owner=physical[surface_key(domain.owner)]
        payloads={Key('owned_free_grass',bundle.id):mesh,Key('owned_free_road',domain.free_section):domain.free_road.mesh,
                  Key('owned_free_owner',domain.owner.id):owner.pavement.mesh}
        if owner.grass is not None:payloads[Key('owned_free_owner_grass',domain.owner.id)]=owner.grass.mesh
        overlap=validate_overlaps(payloads,cell_size=2.)
        token=lambda a,b:tuple(sorted((tuple(a),tuple(b))))
        available={token(pave.xyz[a],pave.xyz[b]) for pave in (owner.pavement.mesh,domain.free_road.mesh)
                   for a,b in MergedRunBody._boundary_vertices(pave)}
        private={token(a,b) for curve in (head,tail) for a,b in zip(curve,curve[1:])}
        for a,b in MergedRunBody._boundary_vertices(mesh):
            if token(mesh.xyz[a],mesh.xyz[b]) not in available|private:
                raise ValueError('Owned/free grass changed an actual whole pavement boundary edge')
        xyz=mesh.xyz.copy();xyz[:,2]=profile.datum+context.network.settings.exaggeration*(xyz[:,2]-profile.datum)
        return Part(bundle.key,Mesh(xyz,mesh.quads.copy(),'grass'),dict(surface.metadata,
            construction='whole original owned member to true free counterpart separator',
            source_sections=list(bundle.sections),physical_owner=domain.owner.text(),source_proof=domain.proof,
            free_closures=domain.closures,private_closure_counts=counts,physical_and_stored_quality=checks,
            maximum_physical_slope=slope,drainage=drainage,overlap=overlap,
            fixed_boundary_height_recovery=height_recovery,
            original_bank_whole_edges=True,source_aliases_created=0,pavement_republished=False,
            runtime_discovery_enabled=True))
