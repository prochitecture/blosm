"""Prescribed hole samples preserve source chords and a conforming quad mesh."""
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.conforming_boundary_refinement import refine,_template
from tested_python_streets.compounds import SampledFootprint,RimSegment,CompoundGroup
from tested_python_streets.mesh import Mesh
from tested_python_streets.model import Key
from tested_python_streets.partition import BoundaryRef
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.regions import Region,BoundaryLoop,SourceSpan


def fixture():
    xyz=np.array([(float(x),float(y),.1*x+.2*y) for y in range(4) for x in range(4)])
    quads=np.array([(4*y+x,4*y+x+1,4*(y+1)+x+1,4*(y+1)+x)
        for y in range(3) for x in range(3) if (x,y)!=(1,1)])
    outer=np.array((0,1,2,3,7,11,15,14,13,12,8,4));hole=np.array((5,9,10,6))
    seams={0:outer,1:hole};loops=[];footprints=[];refs={}
    for li,ids in seams.items():
        owner=Key('fixture',li);xy=xyz[ids,:2]
        polygon=xy if li==0 else xy[::-1]
        labels=tuple(RimSegment('curb',0,0.,1.) for _ in ids)
        footprints.append(SampledFootprint(Region(owner,polygon),labels,SimpleNamespace()))
        loops.append(BoundaryLoop(xy,tuple(SourceSpan(owner,i,0.,1.) for i in range(len(ids)))))
        refs.update({int(v):BoundaryRef(li,i,0.) for i,v in enumerate(ids)})
    mesh=Mesh(xyz,quads,seams=seams).with_uv_layer('UVMap',xyz[quads,:2])
    return mesh,refs,CompoundGroup((Key('intersection',1),),(),tuple(footprints),tuple(loops))


def requested(mesh,parts):
    result=[];ids=mesh.seams[1]
    for i,(a,b) in enumerate(zip(ids,np.roll(ids,-1))):
        result.append(mesh.xyz[a,:2])
        for t in parts.get(i,()):result.append((1.-t)*mesh.xyz[a,:2]+t*mesh.xyz[b,:2])
    return {1:np.asarray(result)[::-1]}


class ConformingBoundaryRefinementTests(unittest.TestCase):
    def test_skew_six_vertex_pavement_rim_accepts_three_direct_quads(self):
        from tested_python_streets.conforming_boundary_refinement import _graded_face
        corners=np.array(((1.25,-.51),(0.,0.),(.39,-.36),(.80,-.46)))
        rim=np.array((corners[0],.5*(corners[0]+corners[1]),corners[1],
                      corners[2],corners[3],.5*(corners[3]+corners[0])))+[1027.,-1128.]
        mesh,details=_graded_face(rim,PatternLibrary(),3,corners.mean(axis=0)+[1027.,-1128.])
        self.assertEqual(len(mesh.quads),3)
        np.testing.assert_array_equal(mesh.xyz[mesh.seams[0],:2],rim)
        self.assertLessEqual(details['stored_quality']['maximum_edge_ratio'],8.)
        self.assertGreaterEqual(details['stored_quality']['minimum_corner_sine'],.03)

    def test_all_small_even_rims_are_positive_cached_direct_quads(self):
        patterns=PatternLibrary()
        for counts in np.ndindex(4,4,4,4):
            counts=tuple(n+1 for n in counts)
            if sum(counts)%2:continue
            pattern=_template(patterns,counts)
            patterns._positive(pattern)
            self.assertEqual(pattern.quads.shape[1],4)
            for side,count in enumerate(counts):self.assertEqual(len(patterns.boundary(pattern,side)),count+1)
            self.assertIs(_template(patterns,counts),pattern)

    def test_adjacent_odd_requests_propagate_inside_without_changing_outer_profiles(self):
        mesh,refs,group=fixture();raw=requested(mesh,{0:(.3,),3:(.7,)})
        output,newrefs,mapping,report=refine(mesh,refs,group,PatternLibrary(),raw)
        np.testing.assert_array_equal(output.xyz[mapping],mesh.xyz)
        np.testing.assert_array_equal(output.seams[0],mapping[mesh.seams[0]])
        actual=output.xyz[output.seams[1],:2][::-1]
        start=int(np.flatnonzero(np.all(actual==raw[1][0],axis=1))[0])
        np.testing.assert_array_equal(np.roll(actual,-start,axis=0),raw[1])
        self.assertGreater(report['propagated_interior_edges'],0)
        self.assertTrue(report['all_direct_quads'])
        edges=np.sort(np.stack((output.quads,np.roll(output.quads,-1,axis=1)),axis=-1),axis=-1)
        unique,counts=np.unique(edges.reshape(-1,2),axis=0,return_counts=True)
        self.assertTrue(set(counts)<={1,2})
        boundary={tuple(e) for e,n in zip(unique,counts) if n==1}
        declared={tuple(sorted((int(a),int(b)))) for ids in output.seams.values() for a,b in zip(ids,np.roll(ids,-1))}
        self.assertEqual(boundary,declared)
        self.assertEqual(set(newrefs),{int(v) for ids in output.seams.values() for v in ids})
        np.testing.assert_allclose(output.uv_layers['UVMap'],output.xyz[output.quads,:2],atol=1e-14)
        output.validate()

    def test_nonuniform_even_subdivisions_keep_refs_and_texture_seams(self):
        mesh,refs,group=fixture();raw=requested(mesh,{0:(.1,.8),2:(.2,.6)})
        output,newrefs,mapping,report=refine(mesh,refs,group,PatternLibrary(),raw)
        self.assertEqual(report['requested_boundary_vertices'],4)
        self.assertEqual(report['propagated_interior_edges'],0)
        for node,ref in newrefs.items():
            loop=group.loops[ref.loop];a,b=loop.xy[ref.edge],loop.xy[(ref.edge+1)%len(loop.xy)]
            np.testing.assert_allclose(output.xyz[node,:2],(1.-ref.fraction)*a+ref.fraction*b,atol=1e-14)
        np.testing.assert_allclose(output.uv_layers['UVMap'],output.xyz[output.quads,:2],atol=1e-14)
        # Face-corner UV discontinuities remain independent of geometric welds.
        shifted=mesh.with_uv_layer('UVMap',mesh.uv_layers['UVMap']+np.arange(len(mesh.quads))[:,None,None]*10.)
        textured,_,_,_=refine(shifted,refs,group,PatternLibrary(),raw)
        self.assertGreater(float(np.ptp(textured.uv_layers['UVMap'][:,:,0])),50.)

    def test_missing_points_off_chord_odd_total_and_budget_fail_transactionally(self):
        mesh,refs,group=fixture();before=mesh.xyz.copy();faces=mesh.quads.copy()
        raw=requested(mesh,{0:(.5,),3:(.5,)})
        variants=[{1:raw[1][1:]},requested(mesh,{0:(.5,)})]
        moved=raw[1].copy();moved[1]+=(.01,.02);variants.append({1:moved})
        for request in variants:
            with self.assertRaises(ValueError):refine(mesh,refs,group,PatternLibrary(),request)
        with self.assertRaisesRegex(ValueError,'cell budget'):refine(mesh,refs,group,PatternLibrary(),raw,max_cells=8)
        with self.assertRaisesRegex(ValueError,'insertion budget'):refine(mesh,refs,group,PatternLibrary(),raw,max_insertions=1)
        np.testing.assert_array_equal(mesh.xyz,before);np.testing.assert_array_equal(mesh.quads,faces)

    def test_no_request_preserves_every_source_point_and_uses_no_odd_transition(self):
        mesh,refs,group=fixture()
        output,newrefs,mapping,report=refine(mesh,refs,group,PatternLibrary(),{})
        self.assertEqual(len(output.quads),len(mesh.quads))
        np.testing.assert_array_equal(output.xyz[mapping],mesh.xyz)
        self.assertEqual(report['relaxation_old_vertex_ids'],[])

    def test_short_prescribed_corner_edges_receive_a_graded_local_collar(self):
        for shift in ((0.,0.),(1200.,-900.)):
            mesh,refs,group=fixture()
            xyz=mesh.xyz.copy();xyz[:,:2]+=shift
            mesh=mesh.with_positions(xyz)
            group=replace(group,loops=tuple(replace(loop,xy=loop.xy+shift) for loop in group.loops))
            raw=requested(mesh,{0:(.02,),3:(.98,)})
            output,_,mapping,report=refine(mesh,refs,group,PatternLibrary(),raw)
            self.assertTrue(report['local_grading'])
            np.testing.assert_array_equal(output.xyz[mapping],mesh.xyz)
            self.assertLessEqual(output.validate()['maximum_edge_ratio'],8.)
            stored=Mesh(output.xyz.astype(np.float32).astype(float),output.quads)
            self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)
            # The source's planar face-corner texture remains exact despite
            # the local atlas having a different number of interior vertices.
            np.testing.assert_allclose(output.uv_layers['UVMap'],output.xyz[output.quads,:2]-shift,atol=1e-10)


if __name__=='__main__':unittest.main()
