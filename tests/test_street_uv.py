"""Numerical contracts for UV scale, seams, payload ownership and strategies."""
import unittest
import numpy as np
import _python_street_package
from tested_python_streets.mesh import Mesh
from tested_python_streets.uv import (DEFAULT_UV_MAP, PlanarMapping, RoadMapping,
                                    apply_mapping, publication_uv_coordinates)


def strip():
    return Mesh(np.array(((-4., -2., .1), (0., -2., .2), (4., -2., .1),
                          (-4., 2., .1), (0., 2., .2), (4., 2., .1))),
                np.array(((0, 1, 4, 3), (1, 2, 5, 4))))


class StreetUVTests(unittest.TestCase):
    def test_planar_metres_negative_coordinates_and_shared_phase(self):
        mesh = apply_mapping(strip(), PlanarMapping())
        uv = mesh.uv_layers[DEFAULT_UV_MAP]
        np.testing.assert_array_equal(uv[0, 0], (-4., -2.))
        np.testing.assert_array_equal(uv[0, 1], uv[1, 0])
        np.testing.assert_array_equal(uv[0, 2], uv[1, 3])
        np.testing.assert_array_equal(np.linalg.norm(np.roll(uv, -1, axis=1) - uv, axis=-1), 4.)

    def test_planar_mapping_uses_explicit_frame_not_part_bounds(self):
        mesh = strip()
        mapping = PlanarMapping(origin=(-4., -2.), rotation=np.pi / 2.)
        uv = mapping.coordinates(mesh)
        np.testing.assert_allclose(uv[0, 0], (0., 0.), atol=1e-14)
        np.testing.assert_allclose(uv[0, 1], (0., -4.), atol=1e-14)
        part = Mesh(mesh.xyz[[1, 2, 4, 5]], np.array(((0, 1, 3, 2),)))
        np.testing.assert_array_equal(mapping.coordinates(part)[0], uv[1])

    def test_road_coordinates_preserve_station_transverse_offset_and_width(self):
        mesh = strip()
        stations = np.array((10., 14., 18., 10., 14., 18.))
        lateral = np.array((-1., -1., -1., 3., 3., 3.))
        uv = RoadMapping(stations, lateral).coordinates(mesh)
        np.testing.assert_array_equal(uv[0], ((10., -1.), (14., -1.), (14., 3.), (10., 3.)))
        # A wider profile changes only transverse extent, never station scale.
        wider = RoadMapping(stations, lateral * 2.).coordinates(mesh)
        np.testing.assert_array_equal(wider[..., 0], uv[..., 0])
        np.testing.assert_array_equal(wider[..., 1], uv[..., 1] * 2.)
        phase = RoadMapping(stations, lateral, station_origin=10., transverse_origin=1.).coordinates(mesh)
        np.testing.assert_array_equal(phase[0, 0], (0., -2.))

    def test_explicit_corner_seam_does_not_split_geometry(self):
        mesh = strip()
        uv = PlanarMapping().coordinates(mesh)
        uv[1, :, 0] += 20.
        mapped = mesh.with_uv_layer('Seamed', uv)
        self.assertEqual(len(mapped.xyz), len(mesh.xyz))
        self.assertEqual(mapped.quads[0, 1], mapped.quads[1, 0])
        self.assertNotEqual(mapped.uv_layers['Seamed'][0, 1, 0], mapped.uv_layers['Seamed'][1, 0, 0])

    def test_owned_immutable_maps_clone_and_memory_accounting(self):
        mesh = strip()
        uv = PlanarMapping().coordinates(mesh)
        mapped = mesh.with_uv_layer('UVMap', uv)
        uv[:] = 99.
        self.assertNotEqual(mapped.uv_layers['UVMap'][0, 0, 0], 99.)
        with self.assertRaises(ValueError):
            mapped.uv_layers['UVMap'][0, 0, 0] = 4.
        self.assertEqual(mapped.nbytes, mesh.nbytes + 2 * 4 * 2 * 8)
        second = mapped.with_uv_layer('Detail', mapped.uv_layers['UVMap'] * .25)
        self.assertEqual(set(second.uv_layers), {'UVMap', 'Detail'})
        np.testing.assert_array_equal(second.uv_layers['UVMap'], mapped.uv_layers['UVMap'])

    def test_vertical_placement_preserves_uvs_and_topology(self):
        mapped = apply_mapping(strip(), PlanarMapping())
        xyz = mapped.xyz.copy()
        xyz[:, 2] += (10., 9., 12., 10., 9., 12.)
        placed = mapped.with_positions(xyz)
        np.testing.assert_array_equal(placed.uv_layers['UVMap'], mapped.uv_layers['UVMap'])
        np.testing.assert_array_equal(placed.quads, mapped.quads)
        np.testing.assert_array_equal(placed.xyz, xyz)
        with self.assertRaisesRegex(ValueError, 'vertex count'):
            mapped.with_positions(xyz[:-1])

    def test_placement_shares_owned_maps_but_defends_against_external_aliases(self):
        mapped = apply_mapping(strip(), PlanarMapping())
        placed = mapped.with_positions(mapped.xyz + (0., 0., 1.))
        self.assertIs(placed.uv_layers['UVMap'], mapped.uv_layers['UVMap'])
        self.assertIsNot(placed.uv_layers, mapped.uv_layers)
        with self.assertRaises(ValueError):
            placed.uv_layers['UVMap'].setflags(write=True)
        external = mapped.uv_layers['UVMap'].copy()
        readonly_view = external.view()
        readonly_view.setflags(write=False)
        mapped.uv_layers['External'] = readonly_view
        copied = mapped.with_positions(mapped.xyz)
        self.assertIsNot(copied.uv_layers['External'], readonly_view)
        external[:] = 99.
        self.assertFalse(np.any(copied.uv_layers['External'] == 99.))
        self.assertIs(copied.uv_layers['UVMap'], mapped.uv_layers['UVMap'])
        detail = copied.with_uv_layer('Detail', copied.uv_layers['UVMap'] * .25)
        self.assertIs(detail.uv_layers['UVMap'], copied.uv_layers['UVMap'])
        self.assertEqual(set(detail.uv_layers), {'UVMap', 'External', 'Detail'})

    def test_publication_rejects_float32_overflow_and_uses_saved_precision(self):
        uv = PlanarMapping().coordinates(strip())
        uv[0, 0, 0] = 1e40
        mesh = strip().with_uv_layer('UVMap', uv)
        mesh.validate()  # finite numerical payload; storage range is separate.
        with self.assertRaisesRegex(ValueError, 'float32'):
            publication_uv_coordinates(mesh.uv_layers['UVMap'])
        valid = np.array((((.1, -2.3), (4.25, 8.75)),), dtype=np.float64)
        published = publication_uv_coordinates(valid[:, ::-1])
        self.assertEqual(published.dtype, np.float32)
        self.assertTrue(published.flags.c_contiguous)
        np.testing.assert_array_equal(publication_uv_coordinates(published.astype(float)), published)
        for invalid in (np.array((0.,)), np.zeros((4, 3)), np.array(((np.nan, 0.),))):
            with self.assertRaisesRegex(ValueError, 'float32'):
                publication_uv_coordinates(invalid)

    def test_fallback_preserves_maps_and_explicit_replacement_is_supported(self):
        mapped = apply_mapping(strip(), PlanarMapping())
        self.assertIs(apply_mapping(mapped, PlanarMapping(origin=(10., 10.))), mapped)
        revised = apply_mapping(mapped, PlanarMapping(origin=(10., 10.)), overwrite=True)
        np.testing.assert_array_equal(revised.uv_layers['UVMap'], mapped.uv_layers['UVMap'] - 10.)

    def test_invalid_uv_payloads_and_strategy_inputs_are_diagnostics(self):
        mesh = strip()
        for shape in ((6, 2), (2, 3, 2), (1, 4, 2)):
            with self.assertRaisesRegex(ValueError, 'face-corner'):
                mesh.with_uv_layer('UVMap', np.zeros(shape))
        bad = np.zeros((2, 4, 2)); bad[0, 0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, 'finite'):
            mesh.with_uv_layer('UVMap', bad)
        for name in ('', 'a\0b', 'a' * 64):
            with self.assertRaisesRegex(ValueError, 'names'):
                mesh.with_uv_layer(name, np.zeros((2, 4, 2)))
        with self.assertRaisesRegex(ValueError, 'origin'):
            PlanarMapping(origin=(0., np.nan)).coordinates(mesh)
        with self.assertRaisesRegex(ValueError, 'stations'):
            RoadMapping([0.], [0.]).coordinates(mesh)
        # Providers may edit the layer dictionary; publication validation catches it.
        mesh.uv_layers['Broken'] = np.zeros((1, 4, 2))
        with self.assertRaisesRegex(ValueError, 'UV'):
            mesh.validate()


if __name__ == '__main__':
    unittest.main()
