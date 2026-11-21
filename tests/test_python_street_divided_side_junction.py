"""A full crown joining a divided road with a clipped far end in Beijing."""
import unittest

import numpy as np

import _python_street_package
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.drainage import BoundaryDrainage
from tested_python_streets.publication_slopes import PublicationSlopeGuard
from tested_python_streets.updates import ProductReady
from tested_python_streets.uv import PlanarMapping, apply_mapping


def side_junction():
    junction = (2049.818603515625, 309.9551696777344)
    slip_end = (2102.93505859375, 299.02691650390625)
    side_end = (2050.284912109375, 290.2851257324219)
    west = (1931.8204345703125, 325.16400146484375)
    east = (2047.525634765625, 326.1849365234375)
    arterial = 10.800000190734863

    def section(sid, street, category, width, points, **ends):
        return model.Section(sid, street, tuple((*p, 0.) for p in points),
            category=category, width=width, rise=.07999999821186066, interpolation=3, **ends)

    sections = (
        section(456, 456, 'primary_link', 2.700000047683716, (junction, slip_end)),
        section(460, 460, 'tertiary', 5.400000095367432,
                (slip_end, (2232.87255859375, 300.2667541503906)), selection_ends=(1,)),
        section(457, 2327, 'residential', 5., (junction, side_end)),
        section(462, 2327, 'residential', 5., (side_end, (2052.184814453125, 210.90371704101562),
                                                (2052.955810546875, 106.50849914550781))),
        section(106, 2445, 'primary', arterial, ((2126.61474609375, 327.94317626953125), east), scope_ends=(0,)),
        section(105, 2445, 'primary', arterial, (east, west)),
        section(107, 2445, 'primary', arterial, (west, (1794.2041015625, 322.347412109375))),
        section(285, 2446, 'primary', arterial, ((1778.01806640625, 304.4100341796875), junction)),
        section(455, 2446, 'primary', arterial, (junction, (2127.046142578125, 312.4142150878906)), scope_ends=(1,)))
    ends = tuple(model.End(*pair) for pair in ((285, 1), (455, 0), (457, 0), (456, 0)))
    site = model.Intersection(912, ends, radius=4., radii=(4., 4., 2., 4.), radius_modes=('STYLE',)*4)
    bundle = model.Bundle(49, (
        model.BundleMember(tuple(model.Traversal(sid, False) for sid in (107, 105, 106))),
        model.BundleMember(tuple(model.Traversal(sid, True) for sid in (285, 455)))))
    connections = (model.Connection(0, (model.End(460, 0), model.End(456, 1))),
                   model.Connection(78, (model.End(457, 1), model.End(462, 0))))
    return model.Network('beijing_side_junction', sections, (site,), (bundle,), connections=connections,
                         settings=model.Settings(cell_size=.699999988079071, exaggeration=1.))


class DividedSideJunctionTests(unittest.TestCase):
    def test_source_profiles_rounded_banks_and_stored_pavement_drain(self):
        network = side_junction()
        # The crop includes the original Bundle paths. Junctions 215 and
        # 219 on its opposite carriageway are outside this surface fixture.
        steps = PythonMeshBackend(network.settings).prepare_steps(network, emit_products=True)
        try:
            surface = next(event.product for event in steps if isinstance(event, ProductReady)
                           and event.key == model.Key('junction_surface', 912))
        finally:
            steps.close()
        part = surface.pavement
        mesh = apply_mapping(part.mesh, PlanarMapping())
        stored = PublicationSlopeGuard.physical_buffer(mesh, part.metadata['datum'], 1.)
        self.assertLessEqual(stored.validate()['maximum_edge_ratio'], 8.)
        self.assertLessEqual(stored.quad_slopes().max(), part.metadata['measured_slope_limit'] + 1e-7)
        fixed = np.zeros(len(stored.xyz), dtype=bool)
        for loop, vertices in stored.seams.items():
            np.testing.assert_allclose(part.mesh.xyz[vertices], surface.physical_rims[loop], rtol=0., atol=1e-10)
            fixed[vertices] = True
        self.assertLessEqual(BoundaryDrainage().fill(stored, fixed)[1]['maximum_fill'], 1e-6)
        self.assertTrue(np.isfinite(mesh.uv_layers['UVMap']).all())
        fillets = surface.group.physical_layout.source_bank_fillets
        self.assertEqual(sorted(fillet.fitted.corner.radius for fillet in fillets), [2., 4.])


if __name__ == '__main__':
    unittest.main()
