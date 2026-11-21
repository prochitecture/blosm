"""Failed-only fixed-rim quad recovery for a finite Bundle pavement half.

An uneven protected separator bank must not dictate a tiny row across the
whole road width. Reuse the semantic direct-quad topology, then the maintained
fixed-rim energy, rotations and interior refinement. Physical bank/profile
vertices, their keys and the dynamic source height field stay authoritative.
"""
import numpy as np
from .bundle_terminal_caps import seed_terminal_chart
from .drainage import BoundaryDrainage
from .four_side_charts import FourSidedQuadChart, quality
from .mesh import Mesh, publication_coordinates
from .model import Key
from .partition import CompoundQuadChart
from .quad_charts import FixedRimQuality
from .quality_correction import QuadShapeCorrection
from .regions import signed_area
from .validation import validate_overlaps


class FixedRimBodyRecovery:
    """Keep the existing finite allocation and independently check each result."""

    def __init__(self, builder, frame, field):
        self.builder, self.frame, self.field = builder, frame, field

    @staticmethod
    def boundary(outer, seam, field):
        physical, kinds, keys = outer
        shared, shared_keys = seam
        physical, shared = np.asarray(physical), np.asarray(shared)
        if (physical.ndim != 2 or physical.shape[1] != 2 or shared.ndim != 2
                or shared.shape[1] != 2 or len(physical) != len(keys)
                or len(kinds) != len(physical)-1 or len(shared) != len(shared_keys)
                or min(len(physical),len(shared)) < 2
                or max(np.linalg.norm(physical[0]-shared[0]),
                       np.linalg.norm(physical[-1]-shared[-1])) > 1e-8):
            raise ValueError('Fixed-rim body needs complete matching physical chains')
        curb = [i for i, kind in enumerate(kinds) if kind[0] == 'outer']
        if not curb:
            raise ValueError('Fixed-rim body has no original outer curb')
        first,last = min(curb),max(curb)+1
        if any(kind[0] != 'outer' for kind in kinds[first:last]):
            raise ValueError('Fixed-rim body has a disconnected outer curb')
        # Shared profiles and banks cannot receive unilateral parity samples.
        if ((last-first-(len(shared)-1)) % 2
                or (first-(len(physical)-1-last)) % 2):
            raise ValueError('Fixed-rim body needs coordinated semantic parity')
        keys = list(keys)
        keys[0],keys[-1] = shared_keys[0],shared_keys[-1]
        rim = np.vstack((physical,shared[-2:0:-1]))
        keys += list(shared_keys[-2:0:-1])
        values = field.values(rim)
        values[first:last+1] = field.datum
        values[len(physical):] = field.ceiling
        values[0] = values[len(physical)-1] = field.ceiling
        corners = (0,first,last,len(physical)-1)
        if signed_area(rim) < 0.:
            rim,keys,values = rim[::-1].copy(),keys[::-1],values[::-1].copy()
            corners = tuple(sorted(len(rim)-1-i for i in corners))
        if len(rim) % 2 or len(set(keys)) != len(keys):
            raise ValueError('Fixed-rim body needs one even, uniquely keyed boundary')
        return rim,keys,values,corners

    def _surface(self, planar, rim, keys, values, name):
        count = len(rim)
        if (not np.array_equal(planar.xyz[:count,:2],rim)
                or not np.array_equal(planar.seams[0],np.arange(count))):
            raise ValueError('Fixed-rim body correction changed a prescribed boundary')
        xyz = planar.xyz.copy()
        xyz[:,2] = self.field.values(xyz[:,:2])
        xyz[:count,2] = values
        mesh = planar.with_positions(xyz)
        edges = np.sort(np.stack((mesh.quads,np.roll(mesh.quads,-1,axis=1)),axis=-1).reshape(-1,2),axis=1)
        unique,incidence = np.unique(edges,axis=0,return_counts=True)
        expected = {tuple(sorted((i,(i+1)%count))) for i in range(count)}
        if (not set(incidence).issubset({1,2})
                or {tuple(edge) for edge in unique[incidence == 1]} != expected
                or len(xyz)-len(unique)+len(mesh.quads) != 1
                or len(np.unique(mesh.quads)) != len(xyz)):
            raise ValueError('Fixed-rim body is not one conforming physical quad disk')
        stored = mesh.with_positions(publication_coordinates(xyz,1.,physical=True))
        measured,published = quality(mesh),quality(stored)
        if (max(measured['maximum_edge_ratio'],published['maximum_edge_ratio']) > 8.
                or min(measured['minimum_corner_sine'],published['minimum_corner_sine']) < .03):
            raise ValueError('Fixed-rim body fails physical or stored pavement grading')
        if (np.any(xyz[:,2] < self.field.datum-1e-10)
                or np.any(xyz[:,2] > self.field.ceiling+1e-10)):
            raise ValueError('Fixed-rim body leaves its original profile height range')
        slopes,stored_slopes = mesh.slope_statistics(),stored.slope_statistics()
        along = max(self.frame.longitudinal_bound(mesh),self.frame.longitudinal_bound(stored))
        if (max(slopes['maximum'],stored_slopes['maximum']) > self.builder.maximum_slope+1e-6
                or along > self.builder.maximum_longitudinal_slope+1e-8):
            raise ValueError('Fixed-rim body exceeds its unchanged physical slope limits')
        fixed = np.zeros(len(xyz),bool)
        fixed[:count] = True
        _,drainage = BoundaryDrainage(max_cells=self.builder.max_cells).fill(mesh,fixed)
        if drainage['maximum_fill'] > 1e-6:
            raise ValueError('Fixed-rim body has an interior depression')
        overlap = validate_overlaps({Key('body_chart',0):mesh},cell_size=max(1.,2.*self.builder.spacing))
        result_keys = list(keys)+[(name,'recovered interior',i) for i in range(count,len(xyz))]
        return mesh,result_keys,dict(world_quality=measured,stored_quality=published,
            physical_slopes=slopes,stored_slopes=stored_slopes,angular_longitudinal_bound=along,
            drainage=drainage,overlap_check=overlap,boundary_vertices=count,
            exact_prescribed_xyz=True,all_direct_quads=True,source_field_unchanged=True,
            exact_whole_boundary_edges=True,conforming_quad_disk=True)

    def build(self, outer, seam, name):
        rim,keys,values,corners = self.boundary(outer,seam,self.field)
        builder = self.builder
        atlas = FourSidedQuadChart(builder.patterns,max_cells=builder.max_cells)
        mapped,counts = atlas.map(self.frame.local(rim),corners)
        mapped,seed = seed_terminal_chart(mapped,len(rim))
        # Reproduce the maintained semantic seed before trying new search
        # order; its previous shape correction is a useful positive seed.
        mapped,previous = QuadShapeCorrection(max_cells=builder.max_cells,max_iterations=250).improve(mapped)
        xy = self.frame.world(mapped.xyz[:,:2])
        xy[:len(rim)] = rim
        initial = Mesh(np.c_[xy,np.zeros(len(xy))],mapped.quads,
            seams={0:np.arange(len(rim),dtype=np.int32)})
        try:
            initial.validate()
        except ValueError:
            from .quad_seed_untangling import positive_seeds
            seeds = positive_seeds(initial,max_cells=builder.max_cells)
        else:
            seeds = ((initial,None),)
        chart = CompoundQuadChart(builder.patterns,spacing=builder.spacing,max_cells=builder.max_cells)
        failures = []
        for seed_index,(initial,untangling) in enumerate(seeds):
            for power in (0,1,2):
                try:
                    if power:
                        candidate,relaxed = FixedRimQuality(iterations=500,power=power).improve(initial)
                    else:
                        candidate,relaxed = initial,None
                    candidate,corrected = chart._correct(candidate,retry_rotations=True)
                    mesh,result_keys,report = self._surface(candidate,rim,keys,values,name)
                    return mesh,result_keys,dict(report,
                        construction='failed-only semantic fixed-rim pavement recovery',
                        rim_intervals=counts,terminal_seed=seed,previous_semantic_correction=previous,
                        positive_seed=untangling,energy_power=power,energy_seed=relaxed,
                        correction=corrected,prior_failures=failures)
                except ValueError as error:
                    failures.append(dict(seed=seed_index,energy_power=power,reason=str(error)))
        raise ValueError('No qualified fixed-rim pavement body recovery: '+str(failures))
