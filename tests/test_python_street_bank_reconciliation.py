"""Physical source banks reconcile without changing original graph ownership."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from _python_street_bank_reconciliation_fixture import fixture
from tested_python_streets.model import Key
from tested_python_streets.paths import Path
from tested_python_streets.source_bank_reconciliation import (
    SourceBankDiscovery,fit_switch,reconcile,signature,validate_prepared_neighbors,
    SourceBankReconciliationError,PreparedBankFillet)
from tested_python_streets.generated_junction_curbs import is_generated_junction_area
from tested_python_streets.compound_ports import CompoundPortPlanner
from tested_python_streets.cycle_dimensional_retry import certify_layout


class BankReconciliationTests(unittest.TestCase):
    def test_unequal_incident_banks_can_cross_beyond_the_narrow_half_width(self):
        from tested_python_streets.model import End, Section
        from tested_python_streets.source_bank_reconciliation import BankOrigin, NativeBankWitness
        planner, group, _, _ = fixture()
        scan = SourceBankDiscovery(planner.context, group)
        site = group.members[0]
        for end in (0,1):
            banks, origins = [], []
            for sid, tip, width, side, distance in (
                    (100,(20.,0.,0.),5.,-1,5.4), (101,(0.,20.,0.),10.8,1,2.5)):
                points = ((0.,0.,0.),tip)
                path = Path(Section(sid,sid,points[::-1] if end else points,width=width,scope_ends=(1-end,)))
                origin = BankOrigin(site,End(sid,end),side)
                limits = (path.length-10.,path.length) if end else (0.,10.)
                banks.append(NativeBankWitness(path,side*(1-2*end),path.length-distance if end else distance,
                    (limits,),(origin,),Key('cycle_port',sid),0,'EXPLICIT_PREPARED_NATIVE_VERTEX_PARAMETERS'))
                origins.append(origin)
            first, last = banks
            self.assertTrue(scan._incident_bank_crossing(first,last,*origins))
            self.assertFalse(scan._incident_bank_crossing(replace(first,station=first.station+.02),last,*origins))
            reaches_cut = replace(last,intervals=((0.,last.path.length),))
            self.assertFalse(scan._incident_bank_crossing(first,reaches_cut,*origins))
            local_cut = Path(replace(last.path.section,scope_ends=(),selection_ends=(end,)))
            self.assertFalse(scan._incident_bank_crossing(first,replace(last,path=local_cut),*origins))

    def test_one_street_can_keep_a_submillimetre_sampled_seam(self):
        from tested_python_streets.model import End, Section
        from tested_python_streets.profiles import Crown
        from tested_python_streets.source_bank_reconciliation import BankOrigin, BankSwitch, NativeBankWitness

        planner, group, _, _ = fixture()
        scan = SourceBankDiscovery(planner.context, group)
        # Two real adjoining Section directions. Their 2 m circular fillet
        # departs from its chord by less than 7 nm.
        points = ((-915.5676879882812,-1026.4765625,0.), (-923.6573486328125,-1014.9426879882812,0.),
                  (-934.7105102539062,-999.1781616210938,0.))
        paths = [Path(Section(100+i,100,points[i:i+2],width=7.)) for i in range(2)]
        site = group.members[0]
        pivot = paths[0].point(paths[0].length,3.5)
        gap = float(np.linalg.norm(pivot-paths[1].point(0.,3.5)))
        ends = (End(100,1),End(101,0))
        banks = []
        for path,end,station in zip(paths,ends,(paths[0].length,gap)):
            sid = path.section.id
            planner.context.products[Key('path',sid)] = path
            planner.context.products[Key('profile',sid)] = Crown(.08,0.,1)
            banks.append(NativeBankWitness(path,1,station,((0.,path.length),),
                (BankOrigin(site,end,1 if not end.end else -1),),Key('cycle_strip',sid),0,
                'EXPLICIT_PREPARED_NATIVE_VERTEX_PARAMETERS'))
        seam = BankSwitch(0,0,'DIRECT_ORIGINAL_BANK_SWITCH',*banks,site,2.,'STYLE',
                          ((ends[0],2.,'STYLE'),),tuple(pivot))
        self.assertTrue(scan.sampled_street_continuation(seam))
        self.assertFalse(scan.sampled_street_continuation(replace(seam,mode='STRICT')))
        self.assertFalse(scan.sampled_street_continuation(replace(seam,last=replace(banks[1],station=1.))))
        self.assertFalse(scan.sampled_street_continuation(replace(seam,radius=20.)))
        planner.context.products[Key('profile',101)] = Crown(.08,0.,-1)
        self.assertFalse(scan.sampled_street_continuation(seam))

    def test_can_cancel_during_bank_fillet_union(self):
        from tested_python_streets.regions import RegionUnion
        from tested_python_streets.source_bank_reconciliation import reconcile_steps
        planner, group, grids, _ = fixture()
        marker = object()
        closed = []
        original = RegionUnion.build_steps

        def observe_union(union, regions):
            regions = tuple(regions)
            if any(r.owner.kind == 'source_bank_fillet' for r in regions):
                try:
                    yield marker
                finally:
                    closed.append(True)
            return (yield from original(union, regions))

        with patch.object(RegionUnion, 'build_steps', observe_union):
            steps = reconcile_steps(planner, group, grids)
            try:
                self.assertTrue(any(event is marker for event in steps))
                previously_completed = len(closed)
            finally:
                steps.close()
        self.assertEqual(len(closed), previously_completed + 1)
        self.assertIsNone(getattr(group.physical_layout, 'source_bank_reconciliation', None))

    def test_exposed_zero_banks_require_a_live_side_entrance_and_style_radii(self):
        from tested_python_streets.model import End
        from tested_python_streets.profiles import Crown
        from tested_python_streets.source_bank_reconciliation import BankOrigin, NativeBankWitness
        planner, group, _, _ = fixture()
        scanner = SourceBankDiscovery(planner.context, group)
        arms = group.physical_layout.approaches
        first, via, last = [End(a.source.id, a.end) for a in arms[:3]]
        banks, origins = [], []
        site = scanner.sites[group.members[0]]
        for end, half in ((first, -1), (last, 0)):
            path = scanner.context[Key('path', end.section_id)]
            station = .1 if end.end == 0 else path.length-.1
            origin = BankOrigin(site.key, end, 1 if end.end == 0 else -1)
            banks.append(NativeBankWitness(path,1,station,(),(origin,),None,0,
                                           'EXPLICIT_PREPARED_NATIVE_VERTEX_PARAMETERS'))
            origins.append(origin)
            scanner.context.products[Key('profile', end.section_id)] = Crown(.08,0.,half)
        walk = [first,via]
        def accepted():
            return scanner._exposed_profile_join(*banks,*origins,walk,site)
        with patch.object(scanner,'tangent',side_effect=lambda bank: np.array((1.,0.) if bank is banks[0] else (.99,-.1))):
            self.assertTrue(accepted())
            banks[0] = replace(banks[0],proof='inferred')
            self.assertFalse(accepted())
            banks[0] = replace(banks[0],proof='EXPLICIT_PREPARED_NATIVE_VERTEX_PARAMETERS')
            with patch.object(type(site),'corner_radius',return_value=(2.,'STRICT')):
                self.assertFalse(accepted())
            scanner.internal.add(via)
            self.assertFalse(accepted())

    def test_only_unfittable_style_tips_on_a_verified_bundle_plateau_are_retained(self):
        from tested_python_streets.corner_radii import CornerContactError
        from tested_python_streets.profiles import Crown
        from tested_python_streets.source_bank_reconciliation import BankSwitch, NativeBankWitness
        planner,group,_,_=fixture()
        scanner=SourceBankDiscovery(planner.context,group)
        paths=[scanner.context[Key('path',sid)] for sid in sorted(scanner.owned)[:2]]
        banks=[NativeBankWitness(path,side,0.,(),(),None,0,'fixture') for path,side in zip(paths,(1,-1))]
        switches=[BankSwitch(1,i,'DIRECT_ORIGINAL_BANK_SWITCH',a,b,group.members[0],2.,'STYLE',(),p)
                  for i,a,b,p in ((0,*banks,(0.,0.)),(2,*reversed(banks),(10.,0.)))]
        loop=SimpleNamespace(hole=True,xy=np.array(((0.,0.),(5.,.2),(10.,0.),(5.,-.2))),spans=(0,0,1,1))
        bundle=SimpleNamespace(id=7,members=(SimpleNamespace(complete=True),SimpleNamespace(complete=True)))
        scanner.group.physical_layout.source_cycle_proofs={7:True}
        roles={bank.path.section.id:(bundle,role,True) for role,bank in enumerate(banks)}
        for bank in banks:
            scanner.context.products[Key('profile',bank.path.section.id)]=Crown(.08,0.,bank.side)
        prefix='tested_python_streets.'
        with patch.object(scanner,'describe',side_effect=lambda span,tail:banks[span]), \
             patch(prefix+'profiles.memberships',return_value=roles), \
             patch(prefix+'source_bank_reconciliation.fit_exposed_switch',
                   side_effect=CornerContactError('no permitted contact')) as fitting:
            self.assertEqual(scanner.natural_separator_tips(loop,switches),[0,2])
            fitting.side_effect=None
            self.assertFalse(scanner.natural_separator_tips(loop,switches))
            fitting.reset_mock()
            self.assertFalse(scanner.natural_separator_tips(loop,[replace(switches[0],mode='STRICT'),switches[1]]))
            loop.spans=(0,1,0,1)
            self.assertFalse(scanner.natural_separator_tips(loop,switches))
            loop.spans=(0,0,1,1)
            loop.hole=False
            self.assertFalse(scanner.natural_separator_tips(loop,switches))
            loop.hole=True
            scanner.group.physical_layout.source_cycle_proofs={}
            self.assertFalse(scanner.natural_separator_tips(loop,switches))
            fitting.assert_not_called()

    def test_actual_finite_banks_keep_originals_and_fresh_width_profile_guards(self):
        for reverse,angle,shift in ((False,0.,(0.,0.)),(True,.7,(12.,-8.))):
            planner,group,grids,others=fixture(reverse=reverse,angle=angle,shift=shift)
            proposals=SourceBankDiscovery(planner.context,group).discover()
            self.assertTrue(proposals)
            result=reconcile(planner,group,grids)
            self.assertIs(result.physical_layout.corners,group.physical_layout.corners)
            self.assertIs(result.physical_layout.approaches,group.physical_layout.approaches)
            self.assertEqual(result.physical_layout.mouths,group.physical_layout.mouths)
            self.assertEqual(result.members,group.members);self.assertEqual(result.claims,group.claims)
            self.assertFalse(CompoundPortPlanner.incomplete_ports(result))
            self.assertLessEqual(certify_layout(planner.context,result,())['maximum_error'],.002)
            self.assertFalse(SourceBankDiscovery(planner.context,result).discover())
            generated=[source for source in result.area_footprints if is_generated_junction_area(source)]
            self.assertTrue(generated)
            self.assertTrue(all(type(source.layout.source_bank_fillet) is PreparedBankFillet for source in generated))
            result.physical_layout.source_parent_chain_reconciliation=dict(
                records=('retained source chain',),prepared_owner_area_guard='PENDING_FINAL_PLAN_INVENTORY')
            qualified,=validate_prepared_neighbors(planner.context,(result,),planner.layouts)
            self.assertEqual(signature(result),signature(qualified))
            self.assertEqual(qualified.physical_layout.generated_curb_prepared_neighbor_guard['tested_pairs'],0)
            chain=qualified.physical_layout.source_parent_chain_reconciliation
            self.assertEqual(chain['records'],('retained source chain',))
            self.assertIs(chain['prepared_owner_area_guard'],qualified.physical_layout.generated_curb_prepared_neighbor_guard)
            self.assertEqual(result.physical_layout.source_parent_chain_reconciliation['prepared_owner_area_guard'],
                             'PENDING_FINAL_PLAN_INVENTORY')
            self.assertIs(reconcile(planner,qualified,grids),qualified)

    def test_original_path_identity_and_disjoint_support_remain_mandatory(self):
        planner,group,grids,_=fixture();context=planner.context
        proposal=SourceBankDiscovery(context,group).discover()[0]
        with self.assertRaisesRegex(ValueError,'original Path identity'):
            fit_switch(context,replace(proposal,first=replace(proposal.first,path=Path(proposal.first.path.section,.005))))
        fitted=fit_switch(context,proposal);arm=fitted.corner.first
        a,b=[value if arm.end==0 else arm.length-value for value in (fitted.corner.stations[0],fitted.intersection[0])]
        middle=.5*(a+b);gap=.1*abs(a-b);ranges=[]
        for low,high in proposal.first.intervals:
            if low<middle<high:ranges.extend(((low,middle-gap),(middle+gap,high)))
            else:ranges.append((low,high))
        with self.assertRaises(ValueError):
            fit_switch(context,replace(proposal,mode='STRICT',first=replace(proposal.first,intervals=tuple(ranges))))

    def test_prepared_neighbor_area_and_opposite_mouth_cannot_be_stolen(self):
        planner,group,grids,others=fixture(neighbor=(4.,4.))
        result=reconcile(planner,group,grids)
        # A prepared owner must supersede its obsolete/raw map entry.
        foreign,=others;layouts=dict(planner.layouts);layouts[foreign.members[0]]=SimpleNamespace(mouths={})
        with self.assertRaisesRegex(SourceBankReconciliationError,'prepared physical owner'):
            validate_prepared_neighbors(planner.context,(result,foreign),layouts)
        far_planner,far,far_grids,far_others=fixture(neighbor=(100.,100.))
        far_result=reconcile(far_planner,far,far_grids)
        self.assertEqual(len(validate_prepared_neighbors(far_planner.context,(far_result,*far_others),far_planner.layouts)),2)

    def test_unrelated_ordinary_groups_and_budgets_are_unchanged(self):
        from tested_python_streets.source_bank_preparation import SourceBankReconciler
        planner,group,grids,_=fixture()
        # No source ancestry means no new inferred bank assignment.
        stripped=[]
        for source in group.footprints:
            layout=copy(source.layout);curves=[]
            for old in layout.boundaries:
                curve=copy(old)
                if hasattr(curve,'source_coordinates'):del curve.source_coordinates
                curves.append(curve)
            layout.boundaries=tuple(curves);stripped.append(replace(source,layout=layout))
        unlabelled=replace(group,footprints=tuple(stripped))
        self.assertIs(reconcile(planner,unlabelled,grids),unlabelled)
        for budget in (0,65,True):
            with self.assertRaises(ValueError):SourceBankReconciler(planner,group,grids,maximum_fillets=budget)

    def test_explicit_strict_input_is_not_replaced_by_a_fitted_style_radius(self):
        from tested_python_streets.updates import BuildContext
        from tested_python_streets.corner_radii import CornerContactError
        planner,group,grids,_=fixture();site,=planner.context.network.intersections
        first=SourceBankDiscovery(planner.context,group).discover()[0]
        endpoint,_,_=first.sectors[0]
        site=site.with_corner_radius(endpoint,2000.,mode='STRICT')
        context=BuildContext(replace(planner.context.network,intersections=(site,)),
            planner.context.products,planner.context.patterns)
        discovery=SourceBankDiscovery(context,group)
        selected=next(proposal for proposal in discovery.discover() if proposal.site==site.key and proposal.mode=='STRICT')
        self.assertEqual(selected.radius,2000.)
        with self.assertRaises(CornerContactError):fit_switch(context,selected)

    def test_generated_geometry_never_adds_original_body_contact_bounds(self):
        from tested_python_streets.source_body_heights import _contact_bounds,_represented
        from tested_python_streets.source_body_references import OriginalBodyReference
        planner,group,grids,_=fixture();context=planner.context
        result=reconcile(planner,group,grids)
        originals=tuple(OriginalBodyReference(context[Key('path',section.id)],
            context[Key('profile',section.id)],0.,context[Key('path',section.id)].length,
            np.array((0.,context[Key('path',section.id)].length))) for section in context.network.sections)
        bounds=_contact_bounds(context,result.area_footprints,originals)
        self.assertTrue(all(bounds[source]==[source.first,source.last] for source in originals))
        generated=[source for source in result.area_footprints if is_generated_junction_area(source)]
        self.assertTrue(generated)
        self.assertTrue(all(not _represented(context,source) for source in generated))
        # An arbitrary Boolean marker cannot hide real original contacts.
        untyped=[]
        for source in generated:
            layout=copy(source.layout);layout.generated_junction_area=True
            untyped.append(replace(source,layout=layout))
        altered=_contact_bounds(context,untyped,originals)
        self.assertTrue(any(altered[source]!=[source.first,source.last] for source in originals))

    def test_final_neighbor_guard_uses_prepared_original_mouth_aliases(self):
        from tested_python_streets.model import End
        from tested_python_streets.routes import Route,Piece,RouteTopology
        from tested_python_streets.intersections import Approach,Mouth
        from tested_python_streets.updates import BuildContext
        from tested_python_streets.compound_cycles import _port_footprint
        from tested_python_streets.regions import RegionUnion
        planner,group,grids,others=fixture(neighbor=(20.,0.));foreign,=others
        context=planner.context;path=context[Key('path',1)]
        arm=Approach(Route((Piece(path,False),)),1,context[Key('profile',1)].oriented(1))
        sources=[]
        for source in foreign.footprints:
            prior=source.layout.approaches[0]
            if prior.source.id==5:
                source=_port_footprint(Key('cycle_port',3),arm,15.,8,.7,.00025,
                    {1:np.linspace(0.,path.length,81)})
            sources.append(source)
        layout=copy(foreign.physical_layout)
        layout.approaches=tuple(arm if prior.source.id==5 else prior for prior in layout.approaches)
        layout.mouths={key:value for key,value in layout.mouths.items() if key[0]!=5}
        layout.mouths[1,1]=Mouth(5.,16)
        foreign=replace(foreign,footprints=tuple(sources),physical_layout=layout,
            loops=RegionUnion().build(source.region for source in sources))
        network=context.network
        sites=tuple(replace(site,ends=tuple(End(1,1) if end.section_id==5 else end for end in site.ends))
                    if site.key in foreign.members else site for site in network.intersections)
        network=replace(network,sections=tuple(section for section in network.sections if section.id!=5),intersections=sites)
        context=BuildContext(network,context.products,context.patterns)
        raw_neighbor=SimpleNamespace(mouths={(1,1):Mouth(19.75,16)})
        planner=SimpleNamespace(context=context,sites={site.key:site for site in sites},
            layouts={group.members[0]:group.physical_layout,foreign.members[0]:raw_neighbor},
            topology=RouteTopology(network))
        result=reconcile(planner,group,grids)
        # The old 19.75 m cut is clear; the eventual owner's real 5 m cut
        # overlaps this owner's 15 m prefix. It must not steal that interval.
        with self.assertRaisesRegex(SourceBankReconciliationError,'prepared neighbor mouth'):
            validate_prepared_neighbors(context,(result,foreign),planner.layouts)


if __name__=='__main__':unittest.main()
