"""Original lead/phase pieces keep source profiles and whole-curve authority."""
from copy import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_run_inventory import prepared,side_junction
from tested_python_streets.ownership import footprint_key,surface_key
from tested_python_streets.model import Key
from tested_python_streets.source_curb_conditions import SourceCurbConditions
from tested_python_streets.source_offset_rules import SourceOffsetRule,SourceCurvePieces,whole_original_curve
from tested_python_streets.surface_fields import BoundaryJetField
from tested_python_streets.compound_surfaces import CompoundSurfaceBuilder
from tested_python_streets.updates import BuildContext


class SourceOffsetRuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.context,_=prepared(side_junction(False))
        cls.product=cls.context[surface_key(Key('intersection',200))]
        cls.group=cls.product.group

    def test_actual_gapped_original_core_keeps_whole_source_and_all_final_guards(self):
        context,group=self.context,self.group;sources={s.region.owner:s for s in group.footprints}
        source=next(iter(sources.values()));before=context.network.to_dict()
        for curve in source.layout.boundaries:
            self.assertTrue(whole_original_curve(context,group,sources,curve))
        conditions,_,limit=CompoundSurfaceBuilder(context.patterns,source_heights=True).boundary_conditions(group,context)
        self.assertLessEqual(BoundaryJetField(conditions).metadata['maximum_sampled_fit_error'],1e-8)
        self.assertTrue(self.product.pavement.metadata['source_height_conditions'])
        self.assertLessEqual(float(self.product.physical_pavement.quad_slopes().max()),limit+1e-8)
        self.assertEqual(context.network.to_dict(),before)

    def test_clipped_union_wrong_scope_or_stale_layout_never_gain_whole_curve_authority(self):
        context,group=self.context,self.group;sources={s.region.owner:s for s in group.footprints}
        owner=group.members[0];source=sources[owner];curve=source.layout.boundaries[0]
        for bad in (replace(group,members=(Key('intersection',999),)),
                    replace(group,claims=(object(),)),replace(group,physical_layout=source.layout),
                    replace(group,footprints=(source,source))):
            self.assertFalse(whole_original_curve(context,bad,sources,curve))
        loops=list(group.loops);spans=list(loops[0].spans)
        target=next(i for i,span in enumerate(spans) if span.label.kind=='curb' and span.label.index==0)
        spans.pop(target);loops[0]=replace(loops[0],spans=tuple(spans))
        self.assertFalse(whole_original_curve(context,replace(group,loops=tuple(loops)),sources,curve))
        products=dict(context.products);products[footprint_key(owner)]=copy(source.layout)
        stale=BuildContext(context.network,products,context.patterns)
        self.assertFalse(whole_original_curve(stale,group,sources,curve))
        original=next(node for node in context.network.intersections if node.key==owner)
        changed=replace(context.network,intersections=tuple(replace(node,radius=node.radius+.1) if node.key==owner else node for node in context.network.intersections))
        self.assertFalse(whole_original_curve(BuildContext(changed,context.products,context.patterns),group,sources,curve))

    def test_cross_contact_chord_uses_declared_original_parts_without_resampling(self):
        context,group=self.context,self.group;sources={s.region.owner:s for s in group.footprints}
        resolver=SourceCurbConditions(context,group);checked=False
        for source in sources.values():
            for curve in source.layout.boundaries:
                if curve.corner.center is None:continue
                contact=float(curve.lengths[1]/curve.length)
                if not 0.<contact<.5:continue
                first,last=0.,min(1.,contact+.01)
                before=curve(np.linspace(0.,1.,101)).copy()
                rule=resolver.derived_rule(curve,first,last)
                self.assertIsInstance(rule,SourceCurvePieces)
                self.assertIn((contact-first)/(last-first),rule.polynomial_breaks())
                for a,b in zip((0.,*rule.polynomial_breaks()),(*rule.polynomial_breaks(),1.)):
                    height,gradient=rule.sample_piece(np.linspace(0.,1.,9),a,b)
                    self.assertTrue(np.isfinite(height).all());self.assertTrue(np.isfinite(gradient).all())
                np.testing.assert_array_equal(curve(np.linspace(0.,1.,101)),before)
                checked=True;break
        self.assertTrue(checked)

    def test_original_offset_jacobian_and_map_knots_are_retained(self):
        context,group=self.context,self.group;source=group.footprints[0]
        checked=0
        for curve in source.layout.boundaries:
            if curve.corner.center is None:continue
            for index,(approach,parameters,cumulative,offset) in enumerate(curve.leads):
                if cumulative[-1]<=1e-8:continue
                lead=SimpleNamespace(approach=approach,parameters=parameters,cumulative=cumulative,
                    offset=offset,length=float(cumulative[-1]))
                rule=SourceOffsetRule(context,lead,0.,1.)
                for value in cumulative[1:-1]/lead.length:self.assertIn(float(value),rule.polynomial_breaks())
                bounds=(0.,*rule.polynomial_breaks(),1.)
                for a,b in zip(bounds,bounds[1:]):
                    station=float(rule._stations(.5*(a+b)));height,gradient=rule.sample_piece(.5,a,b)
                    normal=approach.frame(station)[2]
                    self.assertAlmostEqual(float(gradient@(.5*approach.width*normal)),float(approach.profile.derivative(rule.q)),places=10)
                    self.assertEqual(float(height),approach.profile.datum+float(approach.profile.value(rule.q)))
                checked+=1
        self.assertGreater(checked,0)


if __name__=='__main__':unittest.main()
