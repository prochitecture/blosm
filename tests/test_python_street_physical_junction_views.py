"""Inspection exaggeration must not erase authoritative junction query meshes."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_fixtures import compound_cross_network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_merged_plan import physical_context
from tested_python_streets.ownership import CompoundProduct
from tested_python_streets.updates import BuildContext


class PhysicalJunctionViewTests(unittest.TestCase):
    def test_zero_and_exaggerated_views_keep_physical_pavement_grass_and_release(self):
        source=compound_cross_network(model)
        reference=None
        for exaggeration in (1.,0.,10.):
            network=replace(source,settings=replace(source.settings,exaggeration=exaggeration))
            engine=PythonMeshBackend(network.settings)
            result=engine.prepare(network)
            context=BuildContext(network,result.products,engine.patterns)
            physical=physical_context(context)
            products={key:part for key,part in result.products.items() if isinstance(part,CompoundProduct)}
            self.assertTrue(products)
            current={}
            for key,product in products.items():
                view=physical.products[key]
                self.assertEqual(view.exaggeration,1.)
                self.assertIs(view.pavement.mesh,product.physical_pavement)
                self.assertIsNotNone(product.grass)
                self.assertIs(view.grass.mesh,product.physical_grass)
                current[key]=(view.pavement.mesh.xyz.copy(),view.grass.mesh.xyz.copy())
                if exaggeration==1.:
                    self.assertIs(product.physical_pavement,product.pavement.mesh)
                    self.assertIs(product.physical_grass,product.grass.mesh)
                elif exaggeration==0.:
                    self.assertEqual(float(product.pavement.mesh.xyz[:,2].max()),0.)
                    self.assertGreater(float(view.pavement.mesh.xyz[:,2].max()),.07)
                    self.assertGreater(float(view.grass.mesh.xyz[:,2].max()),.07)
                for li,rim in product.physical_rims.items():
                    np.testing.assert_array_equal(view.rims[li],rim)
                height=product.height(np.array(((0.,0.),)))
                product.release_geometry()
                self.assertIsNone(product.physical_pavement)
                self.assertIsNone(product.physical_grass)
                self.assertIsNone(product.pavement.mesh)
                self.assertIsNone(product.grass.mesh)
                np.testing.assert_array_equal(product.height(np.array(((0.,0.),))),height)
            if reference is None:reference=current
            else:
                self.assertEqual(set(reference),set(current))
                for key,(pavement,grass) in current.items():
                    np.testing.assert_array_equal(pavement,reference[key][0])
                    np.testing.assert_array_equal(grass,reference[key][1])


if __name__=='__main__':unittest.main()
