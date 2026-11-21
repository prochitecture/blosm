"""Original body references retain finite ownership and physical-only outlets."""
from dataclasses import replace
from types import SimpleNamespace
import math
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section,Network,Settings,Bundle,BundleMember,Traversal,Key
from tested_python_streets.paths import Path
from tested_python_streets.profiles import ProfileProvider,memberships,Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.mesh import Mesh
from tested_python_streets.compound_cycles import _footprint
from tested_python_streets.compounds import CompoundGroup
from tested_python_streets.regions import RegionUnion
from tested_python_streets.source_body_references import (
    BodyGuideLimits,BodyWork,OriginalBodyReference,PairedBodyReference,SourceBodyUnavailable,declared_frame)
from tested_python_streets.source_body_heights import prepare,_influence,_represented


def context_for(sections,bundles=()):
    network=Network('finite_original_body',tuple(sections),(),bundles=tuple(bundles),
                    settings=replace(Settings(),cell_size=1.))
    ownership=memberships(network)
    products={}
    for section in network.sections:
        products[Key('path',section.id)]=Path(section,network.settings.curve_tolerance)
        products[Key('profile',section.id)]=ProfileProvider.build(
            section,ownership.get(section.id),network.section_index)
    return BuildContext(network,products,PatternLibrary())


def strip(context,identity):
    path=context[Key('path',identity)];at=np.array((0.,path.length));s=path.section
    points=np.vstack((path.point(at,s.offset+.5*s.width),path.point(at[::-1],s.offset-.5*s.width)))
    return _footprint(Key('cycle_strip',identity),points,internal_edges=(1,3))


def group_for(context,*,foreign=(),ranges=None):
    sources=tuple(strip(context,section.id) for section in context.network.sections)+tuple(foreign)
    if ranges is None:
        ranges={s.id:((0.,context[Key('path',s.id)].length),) for s in context.network.sections}
    layout=SimpleNamespace(source_cycle_full_width_proof=dict(sections=[
        dict(section=sid,interval=(first,last)) for sid,owned in ranges.items() for first,last in owned]))
    return CompoundGroup((),(),sources,RegionUnion().build(s.region for s in sources),layout)


def straight_fixture(*,angle=0.,translation=(0.,0.),reverse=False,foreign=False):
    rotation=np.array(((math.cos(angle),-math.sin(angle)),(math.sin(angle),math.cos(angle))))
    controls=np.array(((0.,0.),(10.,0.)))
    if reverse:controls=controls[::-1]
    controls=controls@rotation.T+translation
    section=Section(1,1,tuple((float(x),float(y),0.) for x,y in controls),width=4.,interpolation=0)
    context=context_for((section,));path=context[Key('path',1)];profile=context[Key('profile',1)]
    at=np.linspace(0.,10.,21);q=np.linspace(1.,-1.,9)
    xy=path.point(at[:,None],2.*q[None,:]).reshape(-1,2)
    z=np.full((len(at),len(q)),.02);z[[0,-1]]=profile.datum+profile.value(q);z[:,[0,-1]]=0.
    faces=np.array([(i*len(q)+j,(i+1)*len(q)+j,(i+1)*len(q)+j+1,i*len(q)+j+1)
                    for i in range(len(at)-1) for j in range(len(q)-1)])
    mesh=Mesh(np.c_[xy,z.ravel()],faces,seams={'actual_end':np.arange(len(q))},
              uv_layers={'StreetUV':xy[faces]})
    fixed=np.zeros((len(at),len(q)),bool);fixed[[0,-1]]=True;fixed[:,[0,-1]]=True
    areas=()
    if foreign:
        points=np.array(((4.2,-.8),(5.8,-.8),(5.8,.8),(4.2,.8)))@rotation.T+translation
        areas=(_footprint(Key('actual_corner_area',8),points),)
    return context,group_for(context,foreign=areas),mesh,np.flatnonzero(fixed.ravel()),at,q


class SourceBodyHeightTests(unittest.TestCase):
    def test_transformed_reversed_full_source_guides_preserve_mesh_contract_and_outlets(self):
        for angle,translation,reverse in ((0.,(0.,0.),False),(.71,(300.,-700.),True)):
            context,group,mesh,boundary,at,q=straight_fixture(angle=angle,translation=translation,reverse=reverse)
            guide=prepare(context,group,mesh,boundary)
            self.assertEqual(len(guide.core_nodes),(len(at)-2)*(len(q)-2))
            self.assertEqual(set(guide.dependencies),{Key('path',1),Key('profile',1)})
            self.assertFalse(guide.ledger['runtime_registered'])
            self.assertEqual(guide.ledger['physical_boundary_nodes'],len(boundary))
            self.assertTrue(guide.ledger['physical_only_outlets'])
            self.assertFalse(guide.physical_fixed[guide.core_nodes].any())
            np.testing.assert_array_equal(np.flatnonzero(guide.physical_fixed),boundary)
            source=context[Key('profile',1)]
            expected=np.tile(source.value(q),(len(at),1)).ravel()[guide.core_nodes]
            np.testing.assert_allclose(guide.core_targets,expected,atol=1e-14,rtol=0.)
            guide.assert_unchanged(mesh)
            self.assertEqual(guide.height_slabs.violation(guide.seed.xyz[:,2]),0.)
            field=guide.physical_field(guide.seed,datum=0.,ceiling=.08)
            self.assertTrue(np.isfinite(field.evaluate(np.array((mesh.xyz[4*len(q)+4,:2],)))).all())
            wrong=guide.seed.with_positions(guide.seed.xyz.copy());wrong.xyz[guide.core_nodes[0],2]+=.001
            with self.assertRaisesRegex(ValueError,'hard body guide'):
                guide.physical_field(wrong,datum=0.,ceiling=.08)
            wrong=guide.seed.with_positions(guide.seed.xyz.copy());wrong.xyz[boundary[0],2]+=.001
            with self.assertRaisesRegex(ValueError,'physical boundary'):
                guide.assert_unchanged(wrong)

    def test_disjoint_owned_ranges_do_not_claim_or_seed_the_unowned_middle(self):
        context,group,mesh,boundary,at,q=straight_fixture()
        group=group_for(context,ranges={1:((0.,4.),(6.,10.))})
        guide=prepare(context,group,mesh,boundary)
        middle=np.flatnonzero((mesh.xyz[:,0]>4.)&(mesh.xyz[:,0]<6.))
        self.assertFalse(set(map(int,middle))&set(map(int,guide.core_nodes)))
        np.testing.assert_array_equal(guide.seed.xyz[middle],mesh.xyz[middle])
        self.assertEqual([r.first for r in guide.references],[0.,6.])
        self.assertEqual([r.last for r in guide.references],[4.,10.])
        overlap=group_for(context,ranges={1:((0.,6.),(4.,10.))})
        with self.assertRaisesRegex(ValueError,'disjoint ranges'):prepare(context,overlap,mesh,boundary)

    def test_foreign_junction_support_blocks_the_whole_row_and_both_seed_and_hard_nodes(self):
        context,group,mesh,boundary,at,q=straight_fixture(foreign=True)
        guide=prepare(context,group,mesh,boundary)
        junction_rows=(mesh.xyz[:,0]>4.)&(mesh.xyz[:,0]<6.)
        self.assertFalse(set(map(int,np.flatnonzero(junction_rows)))&set(map(int,guide.core_nodes)))
        inside=junction_rows&(abs(mesh.xyz[:,1])<=.8)
        self.assertGreater(int(inside.sum()),0)
        np.testing.assert_array_equal(guide.seed.xyz[inside],mesh.xyz[inside])
        influence=guide.ledger['source_influence'][0]
        self.assertEqual(influence['core_ranges'],[(0.,4.),(6.,10.)])
        self.assertEqual(influence['paired_body_sources_excluded'],[])
        self.assertTrue(all(row['owner']=='actual_corner_area:8' for row in influence['influence_witnesses']))

    def test_only_the_declared_opposite_role_is_excluded_from_body_influence(self):
        sections=(Section(1,1,((0.,1.5,0.),(5.,1.5,0.)),width=4.,interpolation=0),
                  Section(2,1,((5.,1.5,0.),(10.,1.5,0.)),width=4.,interpolation=0),
                  Section(3,2,((0.,-1.5,0.),(10.,-1.5,0.)),width=4.,interpolation=0))
        bundle=Bundle(7,(BundleMember((Traversal(1),Traversal(2))),BundleMember((Traversal(3),))))
        context=context_for(sections,(bundle,));group=group_for(context);work=BodyWork(BodyGuideLimits())
        originals=tuple(OriginalBodyReference.build(context,s.id,0.,context[Key('path',s.id)].length,work) for s in sections)
        represented={s.region.owner:_represented(context,s) for s in group.area_footprints}
        cores,reports=_influence(context,group,originals,represented,work)
        first=next(row for row in reports if row['section']==1)
        self.assertEqual(first['paired_body_sources_excluded'],[3])
        self.assertTrue(any(row['owner']=='cycle_strip:2' for row in first['influence_witnesses']))
        self.assertLess(cores[originals[0]][-1][1],5.)

    def test_declared_pairs_have_finite_crown_domains_and_retain_roles_and_remote_dependencies(self):
        for distance in (5.,3.):
            sections=(Section(1,1,((0.,.5*distance,0.),(10.,.5*distance,0.)),width=4.,interpolation=0),
                      Section(2,2,((10.,-.5*distance,0.),(0.,-.5*distance,0.)),width=4.,interpolation=0))
            bundle=Bundle(7,(BundleMember((Traversal(1),)),BundleMember((Traversal(2,False),))))
            context=context_for(sections,(bundle,));work=BodyWork(BodyGuideLimits())
            originals=tuple(OriginalBodyReference.build(context,s.id,0.,10.,work) for s in sections)
            frame,proof=declared_frame(context,bundle)
            pair=PairedBodyReference.build(context,bundle,originals,frame,proof)
            self.assertIs(pair.bundle,bundle)
            self.assertEqual(pair.forward,(True,False))
            a,b=originals
            ranges={a:((1.,3.),(5.,9.)), b:((1.,4.),(5.,8.))}
            # In the shared direction the second road's ranges are [2,5] and [6,9].
            # The contact at station 5 alone is not a positive-length body.
            intervals=pair.whole_body_intervals(ranges)
            expected=frame.local(a.path.point(np.array([2.,3.,6.,9.]),0.))[:,0]
            np.testing.assert_allclose(intervals,expected.reshape(2,2),atol=1e-12)
            self.assertEqual(pair.whole_body_intervals({a:((0.,1.),),b:((0.,1.),)}),())
            self.assertEqual(set(pair.dependencies),{bundle.key,Key('path',1),Key('profile',1),Key('path',2),Key('profile',2)})
            points=np.mean(np.stack([source.path.point(np.array((5.,)),source.path.section.offset)
                                     for source in originals]),axis=0)
            self.assertAlmostEqual(float(pair.values(points)[0]),.08)
            with self.assertRaisesRegex(SourceBodyUnavailable,'extrapolate'):
                pair.values(frame.world(np.array(((pair.last+.001,0.),))))
            wrong=replace(originals[1],profile=replace(originals[1].profile,half=-originals[1].profile.half))
            with self.assertRaisesRegex(SourceBodyUnavailable,'facing half-profile'):
                PairedBodyReference.build(context,bundle,(originals[0],wrong),frame,proof)
            for profile in (replace(originals[1].profile,rise=.09),replace(originals[1].profile,datum=.01)):
                with self.assertRaisesRegex(SourceBodyUnavailable,'common datum and rise'):
                    PairedBodyReference.build(context,bundle,(originals[0],replace(originals[1],profile=profile)),frame,proof)

    def test_wrong_original_profiles_false_outlets_and_unbounded_work_are_rejected(self):
        context,group,mesh,boundary,at,q=straight_fixture()
        original=context.products[Key('profile',1)]
        context.products[Key('profile',1)]=Crown(.09,0.)
        with self.assertRaisesRegex(SourceBodyUnavailable,'original Path or Profile'):
            prepare(context,group,mesh,boundary)
        context.products[Key('profile',1)]=original
        for false in (boundary[:-1],np.r_[boundary,10],boundary.astype(float)):
            with self.assertRaisesRegex(ValueError,'physical.*boundary'):prepare(context,group,mesh,false)
        for kwargs in (dict(maximum_samples=True),dict(maximum_cells=8193),dict(source_tolerance=.001),dict(maximum_source_error=.003)):
            with self.assertRaises(ValueError):BodyGuideLimits(**kwargs)
        with self.assertRaisesRegex(SourceBodyUnavailable,'samples budget'):
            prepare(context,group,mesh,boundary,limits=BodyGuideLimits(maximum_samples=2))
        with self.assertRaisesRegex(SourceBodyUnavailable,'classification_work budget'):
            prepare(context,group,mesh,boundary,limits=BodyGuideLimits(maximum_classification_work=1))


if __name__=='__main__':unittest.main()
