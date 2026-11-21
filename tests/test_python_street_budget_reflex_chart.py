"""An acute half-profile fork needs bounded strip charts, not dense radial layers."""
from dataclasses import replace
import unittest
from unittest.mock import patch
import numpy as np
import _python_street_package
from tested_python_streets.model import Network,Section,Intersection,End,Key,Settings
from tested_python_streets.paths import Path
from tested_python_streets.routes import Route,Piece
from tested_python_streets.profiles import Crown
from tested_python_streets.patterns import PatternLibrary
from tested_python_streets.updates import BuildContext
from tested_python_streets.intersections import Approach,IntersectionProvider
from tested_python_streets.ownership import PLAN,OwnershipPlan,JunctionOwnership,footprint_key
from tested_python_streets.compounds import CompoundGroup,sample_footprint
from tested_python_streets.regions import RegionUnion
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.reflex_cuts import ReflexCurbCuts
from tested_python_streets.quad_charts import QuadChartCellBudgetError,QuadChartLayerBudgetError


def fixture():
    # Resolved local source geometry from a real acute fork. Source identities
    # are deliberately local; no runtime construction knows an extract or ID.
    bend=((1413.0421142578125,19.429636001586914),(1408.093017578125,27.3099365234375),
          (1406.272216796875,29.51365089416504),(1404.031494140625,31.19407081604004),
          (1365.7080078125,54.161895751953125),(1341.519775390625,72.51322937011719),
          (1334.936767578125,82.45264434814453),(1331.40283203125,93.85100555419922),
          (1330.1324462890625,105.01607513427734),(1332.47900390625,117.17266845703125),
          (1337.929443359375,128.75106811523438),(1344.8275146484375,136.42245483398438),
          (1352.993896484375,142.01246643066406),(1364.1856689453125,149.88519287109375))
    shapes=(bend,(bend[0],(1121.023193359375,199.20689392089844)),
            ((1433.4566650390625,6.643689155578613),bend[0]))
    rise=float(np.float32(.08))
    sections=tuple(Section(i+1,i+1,tuple((x,y,0.) for x,y in points),
        width=3.5 if i==0 else 7.,rise=rise) for i,points in enumerate(shapes))
    site=Intersection(1,(End(1,0),End(2,0),End(3,1)),radius=10.,radius_modes=('STYLE',)*3)
    network=Network('acute_resolved_half_profile_fork',sections,(site,),settings=Settings(exaggeration=1.))
    products={PLAN:OwnershipPlan((),{})}
    for section in sections:
        products[Key('path',section.id)]=Path(section)
        products[Key('profile',section.id)]=Crown(rise,0.,0 if section.id==1 else -1)
    context=BuildContext(network,products,PatternLibrary())
    approaches=[Approach(Route((Piece(products[Key('path',end.section_id)],end.end==0),)),
                         end.end,products[Key('profile',end.section_id)])
                for end in site.ends]
    physical=IntersectionProvider.curved_boundary(context,site,approaches)
    products[footprint_key(site.key)]=physical
    try:layout=IntersectionProvider.layout(context,site,physical)
    except ValueError:layout=replace(physical,mesh_strategy='boundary_chart')
    products[site.layout_key]=layout
    return context,site


class BudgetReflexChartTests(unittest.TestCase):
    def test_layer_budget_enters_the_existing_real_reflex_surface_family(self):
        context,site=fixture();before=context.network.to_dict()
        original=CompoundQuadChart._build_steps;attempts=[]
        def typed_layer_failure(chart,group,ctx,*,require_quality):
            if not chart.four_sided:
                attempts.append('radial_layer_budget')
                raise QuadChartLayerBudgetError('Direct quad chart exceeds its layer budget')
            attempts.append('existing_reflex_chart')
            return original(chart,group,ctx,require_quality=require_quality)
        with patch.object(CompoundQuadChart,'_build_steps',typed_layer_failure):
            product=JunctionOwnership.surface(context,site.key)
        self.assertIn('existing_reflex_chart',attempts)
        self.assertIn('budget_reflex_fallback',product.pavement.metadata['planar_chart'])
        stored=product.pavement.mesh.with_positions(product.pavement.mesh.xyz.astype(np.float32).astype(float))
        self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)
        self.assertEqual(before,context.network.to_dict())

    def test_radial_budget_failure_uses_same_fixed_rim_and_keeps_surface_contracts(self):
        context,site=fixture();original=context.network.to_dict()
        footprint=sample_footprint(site.key,context[footprint_key(site.key)])
        group=CompoundGroup((site.key,),(),(footprint,),RegionUnion().build((footprint.region,)))
        self.assertFalse(ReflexCurbCuts.eligible(group))
        chart=CompoundQuadChart(context.patterns,spacing=context.network.settings.spacing)
        with self.assertRaises(QuadChartCellBudgetError):chart._build(group,context,require_quality=True)
        product=JunctionOwnership.surface(context,site.key)
        metadata=product.pavement.metadata
        self.assertIn('budget_reflex_fallback',metadata['planar_chart'])
        self.assertLess(len(product.pavement.mesh.quads),8000)
        stored=product.pavement.mesh.with_positions(product.pavement.mesh.xyz.astype(np.float32).astype(float))
        self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)
        self.assertEqual(metadata['source_cycle_reconciliation']['joins']['isolated_interior_minima'],0)
        self.assertLessEqual(metadata['source_cycle_reconciliation']['joins']['maximum_entrance_longitudinal_slope'],.03)
        for rim in product.physical_rims.values():
            self.assertGreaterEqual(float(rim[:,2].min()),-1e-10)
            self.assertLessEqual(float(rim[:,2].max()),float(np.float32(.08))+1e-10)
        self.assertEqual(context.network.to_dict(),original)


if __name__=='__main__':unittest.main()
