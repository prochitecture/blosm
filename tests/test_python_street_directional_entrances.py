"""An asymmetric four-arm junction must keep profiles and gentle entrances."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.model import Section,End,Intersection,Network,Settings,Key
from tested_python_streets.mesh import publication_coordinates
from tested_python_streets.ownership import surface_key
from tested_python_streets.four_side_charts import quality


SHAPES=(
 ((475.5946960449219,490.6888732910156),(611.8163452148438,496.1198425292969)),
 ((475.5946960449219,490.6888732910156),(478.2325439453125,501.87664794921875),(476.8374938964844,552.8497924804688)),
 ((470.6297302246094,552.6489868164062),(472.2808532714844,502.3549499511719),(475.5946960449219,490.6888732910156)),
 ((475.5946960449219,490.6888732910156),(449.7989501953125,482.9283447265625)))


class DirectionalEntranceTests(unittest.TestCase):
    def test_actual_source_controls_pass_the_normal_backend_without_profile_changes(self):
        rise=float(np.float32(.08))
        sections=tuple(Section(i,i,tuple((x,y,0.) for x,y in xy),width=5. if i==4 else 2.,
            rise=rise,interpolation=3,category='residential' if i==4 else 'service')
            for i,xy in enumerate(SHAPES,1))
        site=Intersection(1,(End(1,0),End(2,0),End(3,1),End(4,0)),radius=1.,
            radii=(1.,1.,1.,2.),radius_modes=('STYLE',)*4)
        network=Network('asymmetric_service_entrances',sections,(site,),settings=Settings(exaggeration=1.))
        before=network.to_dict()
        result=PythonMeshBackend(network.settings).prepare(network)
        product=result.products[surface_key(site.key)]
        metadata=product.pavement.metadata
        recovery=metadata['deferred_projection_reconciliation']['directional_entrance_recovery']
        self.assertLessEqual(recovery['joins']['maximum_entrance_longitudinal_slope'],.03)
        self.assertEqual(recovery['joins']['isolated_interior_minima'],0)
        self.assertEqual(len(metadata['deferred_projection_reconciliation']['discarded_refinement_attempts']),3)
        self.assertEqual(network.to_dict(),before)
        mesh=product.pavement.mesh
        stored=mesh.with_positions(publication_coordinates(mesh.xyz,1.,physical=True))
        for candidate in (mesh,stored):
            shape=quality(candidate)
            self.assertLessEqual(shape['maximum_edge_ratio'],8.)
            self.assertGreaterEqual(shape['minimum_corner_sine'],.03)
            self.assertLessEqual(float(candidate.quad_slopes().max()),.176+1e-7)
        for sid,end,u in product.port_points:
            loop,node,station=product.port_points[sid,end,u]
            expected=rise*(1.-(2.*u-1.)**2)
            self.assertAlmostEqual(product.rims[loop][node,2],expected,places=12)
        self.assertEqual(set(sid for sid,end in product.peaks),{1,2,3,4})
        published=result.products[Key('junction_mesh',site.id)].mesh
        self.assertTrue(np.isfinite(published.uv_layers['UVMap']).all())


if __name__=='__main__':unittest.main()
