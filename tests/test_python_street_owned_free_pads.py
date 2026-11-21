"""A tiny true dead end can belong wholly to a smooth through-junction pool."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.model import Section, End, Intersection, Network, Key
from tested_python_streets.backend import PythonMeshBackend


def original_pad():
    # Original residential junction controls, translated into a local metre plane.
    return Network('finite_free_pad', (
        Section(1, 1, ((-12.413330078125, 25.61102294921875, 0.), (0., 0., 0.)), width=5.),
        Section(2, 2, ((0., 0., 0.), (-2.2708740234375, -.92462158203125, 0.)), width=5.),
        Section(3, 3, ((0., 0., 0.), (41.0765380859375, -82.07501220703125, 0.)), width=5.)),
        (Intersection(1, (End(1, 1), End(2, 0), End(3, 0)), radius=2.,
            radius_modes=('STYLE',)*3),))


class OwnedFreePadTests(unittest.TestCase):
    def test_narrow_dead_end_is_contained_in_a_wider_bent_through_pool(self):
        # Same finite geometry as the source-backed service-stub regression,
        # translated locally; the original rectangle is not shortened.
        network=Network('narrow_bent_pad',(
            Section(1,1,((7.1761474609375,-9.38232421875,0.),(0.,0.,0.)),width=7.),
            Section(2,2,((0.,0.,0.),(.2325439453125,3.172607421875,0.)),width=2.),
            Section(3,3,((0.,0.,0.),(-261.468994140625,68.2178955078125,0.)),width=7.)),
            (Intersection(1,(End(1,1),End(2,0),End(3,0)),radius=2.,
                radii=(2.,1.,2.),radius_modes=('STYLE',)*3),))
        before=network.to_dict()
        update=PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(),before)
        layout=update.products[Key('mouths',1)]
        claim,=layout.finite_free_pads
        self.assertEqual(claim.section,2)
        self.assertGreater(claim.length,3.)
        self.assertGreater(claim.proof_margin,.02)
        self.assertIsNone(getattr(update.products[Key('road',2)],'mesh',None))
        part=update.products[Key('junction_mesh',1)]
        stored=part.mesh.with_positions(part.mesh.xyz.astype(np.float32).astype(float))
        self.assertLessEqual(stored.validate()['maximum_edge_ratio'],8.)
        self.assertTrue(part.metadata['finite_port_reconciliation']['exact_boundary'])
        joins=part.metadata['finite_port_reconciliation']['joins']
        self.assertEqual(joins['isolated_interior_minima'],0)
        self.assertLessEqual(joins['maximum_entrance_longitudinal_slope'],.03)

    def test_section_reversal_and_exact_sources_keep_the_same_whole_pad_owner(self):
        source=original_pad()
        for reversed_ids,mode in (({2},3),({1,2,3},0)):
            with self.subTest(reverse=reversed_ids,interpolation=mode):
                sections=tuple(replace(s,points=s.points[::-1] if s.id in reversed_ids else s.points,
                    interpolation=mode) for s in source.sections)
                swap=lambda e:End(e.section_id,1-e.end) if e.section_id in reversed_ids else e
                network=replace(source,sections=sections,intersections=tuple(replace(j,
                    ends=tuple(map(swap,j.ends))) for j in source.intersections))
                original=network.to_dict()
                update=PythonMeshBackend(network.settings).prepare(network)
                self.assertEqual(network.to_dict(),original)
                layout=update.products[Key('mouths',1)]
                claim,=layout.finite_free_pads
                self.assertEqual(claim.section,2)
                self.assertEqual(claim.end,1)
                self.assertEqual(layout.mouths[2,1].station,0.)
                self.assertIsNone(getattr(update.products[Key('road',2)],'mesh',None))
                part=update.products[Key('junction_mesh',1)]
                rounded=part.mesh.with_positions(part.mesh.xyz.astype(np.float32).astype(float))
                self.assertLessEqual(rounded.validate()['maximum_edge_ratio'],8.)
                self.assertLessEqual(float(rounded.quad_slopes(network.settings.exaggeration).max()),.1+1e-7)
                check=part.metadata['finite_port_reconciliation']
                self.assertTrue(check['exact_boundary'])
                self.assertEqual(check['joins']['isolated_interior_minima'],0)
                self.assertLessEqual(check['joins']['maximum_entrance_longitudinal_slope'],.03)

    def test_selection_cut_is_not_reclassified_as_a_free_pad(self):
        source=original_pad()
        selected=replace(source,sections=(source.sections[0],
            replace(source.sections[1],selection_ends=(1,)),source.sections[2]))
        with self.assertRaises(ValueError):PythonMeshBackend(selected.settings).prepare(selected)

    def test_original_short_pad_has_one_complete_physical_owner(self):
        network = original_pad()
        before = network.to_dict()
        prepared = PythonMeshBackend(network.settings).prepare(network)
        self.assertEqual(network.to_dict(), before)
        layout = prepared.products[Key('mouths', 1)]
        claim, = layout.finite_free_pads
        self.assertEqual(claim.section, 2)
        self.assertGreater(claim.proof_margin, .02)
        self.assertGreater(layout.mouths[2, 0].station, 2.4)
        self.assertIsNone(getattr(prepared.products[Key('road', 2)], 'mesh', None))
        part = prepared.products[Key('junction_mesh', 1)]
        self.assertLessEqual(part.mesh.validate()['maximum_edge_ratio'], 8.)
        self.assertGreaterEqual(float(part.mesh.xyz[:, 2].min()), -1e-8)
        self.assertLessEqual(float(part.mesh.xyz[:, 2].max()), .08*network.settings.exaggeration+1e-8)

    def test_known_scope_end_is_not_reclassified_as_a_free_pad(self):
        source = original_pad()
        scoped = replace(source, sections=(source.sections[0],
            replace(source.sections[1], scope_ends=(1,)), source.sections[2]))
        with self.assertRaises(ValueError):
            PythonMeshBackend(scoped.settings).prepare(scoped)

    def test_complete_pad_and_profiles_survive_rotated_kilometre_coordinates(self):
        source = original_pad()
        angle = .713
        rotation = np.array(((np.cos(angle), -np.sin(angle)), (np.sin(angle), np.cos(angle))))
        sections = []
        for s in source.sections:
            xy = np.asarray(s.points)[:, :2] @ rotation.T + (1795.819824, 679.240967)
            sections.append(replace(s, points=tuple((x, y, 0.) for x, y in xy)))
        network = replace(source, sections=tuple(sections))
        update = PythonMeshBackend(network.settings).prepare(network)
        part = update.products[Key('junction_mesh', 1)]
        published = part.mesh.with_positions(part.mesh.xyz.astype(np.float32).astype(float))
        self.assertLessEqual(published.validate()['maximum_edge_ratio'], 8.)
        layout = update.products[Key('mouths', 1)]
        self.assertEqual(len(layout.finite_free_pads), 1)
        for approach, station in zip(layout.approaches, layout.outward_stations):
            q = np.linspace(-1., 1., 2*layout.radial+1)
            expected = np.column_stack((approach.point(station, .5*approach.width*q),
                approach.profile.datum + approach.profile.value(q)*network.settings.exaggeration))
            distances = np.linalg.norm(part.mesh.xyz[:, None]-expected[None], axis=-1)
            self.assertLess(float(distances.min(axis=0).max()), 1e-8)
        self.assertLessEqual(part.metadata['physical_slopes']['maximum'], .1)

    def test_curved_envelope_preserves_end_jets_and_certifies_curvature(self):
        from tested_python_streets.continuations import WidthContinuation
        from tested_python_streets.curb_envelope import CurbEnvelope
        base = WidthContinuation((-5., 0.), (5., 0.), (1., 0.), (1., 0.), 2.)
        curve = CurbEnvelope(base, (0., 1.), .5, 2.)
        np.testing.assert_array_equal(curve([0., 1.]), base([0., 1.]))
        np.testing.assert_allclose(curve.derivative([0., 1.]), base.derivative([0., 1.]), atol=1e-13)
        np.testing.assert_allclose(curve.second_derivative[[0, -1]], base.second_derivative[[0, -1]], atol=1e-12)
        self.assertAlmostEqual(curve(.5)[1], .5)
        t = np.linspace(0., 1., 10001)
        velocity = curve.derivative(t)
        acceleration = np.gradient(velocity, t, axis=0, edge_order=2)
        curvature = abs(velocity[:, 0]*acceleration[:, 1]-velocity[:, 1]*acceleration[:, 0]) / np.linalg.norm(velocity, axis=1)**3
        self.assertLessEqual(float(curvature.max()), curve.curvature_bound+1e-7)
        self.assertLessEqual(curve.curvature_bound, .5)
        for radius in (0., -1., float('nan'), float('inf')):
            with self.subTest(radius=radius), self.assertRaises(ValueError):
                CurbEnvelope(base, (0., 1.), .5, radius)
        for outward,distance in (((0.,0.),.5),((0.,2.),.5),((0.,float('nan')),.5),
                                  ((0.,1.),-.1),((0.,1.),float('inf'))):
            with self.subTest(outward=outward,distance=distance),self.assertRaises(ValueError):
                CurbEnvelope(base,outward,distance,2.)
        # Degree elevation is defined for the original quintic only. Reusing
        # an already elevated envelope would silently drop a control point.
        with self.assertRaises(ValueError):CurbEnvelope(curve,(0.,1.),.1,2.)

    def test_free_pool_never_bypasses_an_absorbed_dimensional_seam_reserve(self):
        from tested_python_streets.free_pads import _candidate
        from tested_python_streets.intersections import Approach, CornerConstructionError
        from tested_python_streets.paths import Path
        from tested_python_streets.profiles import Crown
        from tested_python_streets.routes import Piece, Route
        from types import SimpleNamespace
        n = original_pad()
        approaches = [Approach(Route((Piece(Path(s), e.end == 0),), terminal_free=True),
            e.end, Crown(s.rise, s.datum)) for s, e in zip(n.sections, n.intersections[0].ends)]
        context = SimpleNamespace(network=n)
        failure = CornerConstructionError('test', n.intersections[0], approaches)
        self.assertIsNotNone(_candidate(context, approaches, failure))
        approaches[0].path.minimum_mouth = 12.
        self.assertIsNone(_candidate(context, approaches, failure))

    def test_local_update_and_furniture_follow_the_one_owned_surface(self):
        from tested_python_streets.model import Placement
        from tested_python_streets.ownership import surface_key
        source = replace(original_pad(), furniture=(Placement(1, 2, 1., offset=-2.5),))
        engine = PythonMeshBackend(source.settings)
        initial = engine.prepare(source)
        engine.commit(initial)
        pad = source.sections[1]
        edited_pad = replace(pad, points=(pad.points[0], (-2.17, -.884, 0.)))
        edited = replace(source, sections=(source.sections[0], edited_pad, source.sections[2]))
        local = engine.prepare(edited)
        fresh = PythonMeshBackend(edited.settings).prepare(edited)
        for key, part in local.products.items():
            if getattr(part, 'mesh', None) is not None:
                np.testing.assert_array_equal(part.mesh.xyz, fresh.products[key].mesh.xyz)
                np.testing.assert_array_equal(part.mesh.quads, fresh.products[key].mesh.quads)
        anchor = local.products[Key('placement', 1)]
        query = local.products[surface_key(Key('intersection', 1))]
        self.assertAlmostEqual(anchor.position[2], float(query.height(np.asarray([anchor.position[:2]]))[0]))
        engine.commit(local)
        restored = engine.prepare(source)
        np.testing.assert_array_equal(restored.products[Key('junction_mesh', 1)].mesh.xyz,
                                      initial.products[Key('junction_mesh', 1)].mesh.xyz)


if __name__ == '__main__':
    unittest.main()
