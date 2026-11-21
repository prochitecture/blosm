"""Long grass separators with exact pavement-adjacent boundary samples."""
import numpy as np
from .four_side_charts import FourSidedQuadChart, quality
from .hole_banded_chart import _Bands, audit
from .mesh import Mesh, publication_coordinates
from .positive_quad_patches import positive_patches
from .quad_charts import ChartResult, FixedRimQuality, LayeredQuadChart, QuadChartQualityError
from .quality_correction import QuadShapeCorrection
from .regions import signed_area
from .separator_surfaces import PartitionedSeparatorChart


def _qualified(mesh):
    try:
        physical = quality(mesh)
        stored = quality(mesh.with_positions(publication_coordinates(mesh.xyz, 1., physical=True)))
    except ValueError:
        return False
    return (max(physical['maximum_edge_ratio'], stored['maximum_edge_ratio']) <= 8.
            and min(physical['minimum_corner_sine'], stored['minimum_corner_sine']) > 0.)


def _correct(mesh, max_cells):
    return QuadShapeCorrection(max_iterations=1000, max_cells=max_cells,
        minimum_sine=1e-6, target_sine=1e-5).improve(mesh)


class _GrassLayered(LayeredQuadChart):
    def build(self, boundary):
        try:
            return super().build(boundary)
        except QuadChartQualityError as previous:
            candidate, correction = _correct(previous.result.mesh, self.max_cells)
            if not _qualified(candidate):
                raise ValueError('Fixed-rim grass band still needs another positive chart') from previous
            return ChartResult(candidate, dict(previous.result.metadata,
                positive_grass_correction=correction))


class _GrassFour(FourSidedQuadChart):
    def build(self, boundary):
        seeds = []
        options = self.corners(boundary, include_smooth_corners=True,
            minimum_polygon_sine=.001, minimum_corner_angle=.001)
        for score, corners in options[:24]:
            initial = None
            try:
                initial, counts = self.map(boundary, corners)
                report = quality(initial)
            except ValueError:
                if initial is None:
                    continue
                try:
                    variants = positive_patches(initial, max_cells=self.max_cells)
                except ValueError:
                    continue
                seeds.extend((quality(mesh)['maximum_edge_ratio'], mesh, counts,
                              dict(atlas_corners=corners, positive_patch=proof))
                             for mesh, proof in variants)
            else:
                seeds.append((report['maximum_edge_ratio'], initial, counts, dict(atlas_corners=corners)))
        for ratio, initial, counts, construction in sorted(seeds, key=lambda item: item[0])[:4]:
            for power in (1, 2):
                candidate, energy = FixedRimQuality(iterations=500, power=power).improve(initial)
                candidate, correction = _correct(candidate, self.max_cells)
                if _qualified(candidate):
                    return ChartResult(candidate, dict(construction,
                        construction='fixed-rim positive grass band', quality=quality(candidate),
                        rim_intervals=counts, energy=energy, correction=correction))
        raise ValueError('No qualified bounded positive grass band atlas')


class BandedSeparatorChart(PartitionedSeparatorChart):
    """A failed-only family: at most 64 bands, without adding any rim point."""
    def __init__(self, *args, bands=64, **kwargs):
        if type(bands) is not int or not 2 <= bands <= 64:
            raise ValueError('Fixed-rim separator needs two to sixty-four computational bands')
        super().__init__(*args, **kwargs)
        self.bands = bands

    def cut_policy(self):
        return _Bands(self.bands)

    def block_builders(self, remaining):
        return (*super().block_builders(remaining),
            ('positive-grass-layered', _GrassLayered(self.patterns, spacing=self.interior_spacing,
                max_cells=remaining, quality_iterations=100)),
            ('positive-grass-four', _GrassFour(self.patterns, max_cells=remaining)))

    def build(self, boundary):
        boundary = np.asarray(boundary, dtype=float)
        if (boundary.ndim != 2 or boundary.shape[1] != 2 or not 4 <= len(boundary) <= 8192
                or len(boundary) % 2 or not np.isfinite(boundary).all()
                or signed_area(boundary) <= 0. or not 4 <= self.max_cells <= 50000):
            raise ValueError('Invalid bounded fixed-rim separator inputs')
        result = super().build(boundary)
        result.mesh.seams = {0: np.arange(len(boundary), dtype=np.int32)}
        certificate = audit(result.mesh, boundary)
        if certificate['added_boundary_samples'] or not np.array_equal(result.mesh.xyz[:len(boundary), :2], boundary):
            raise ValueError('Standalone grass bands must retain every whole pavement-adjacent edge')
        if len(result.mesh.quads) > self.max_cells:
            raise ValueError('Standalone grass bands exceed their unchanged allocation budget')
        result.metadata = dict(result.metadata, construction='exact-rim bounded direct-quad separator bands',
            bands=self.bands, exact_whole_boundary_edges=True, rim_certificate=certificate)
        return result


class GradedRimSeparatorChart(BandedSeparatorChart):
    """One interior row with extra quads where adjacent fixed edge lengths jump."""

    def build(self, boundary):
        from .model import Key
        from .validation import validate_overlaps
        boundary = np.asarray(boundary, dtype=float)
        if (boundary.ndim != 2 or boundary.shape[1] != 2 or not 4 <= len(boundary) <= 8192
                or len(boundary) % 2 or not np.isfinite(boundary).all() or signed_area(boundary) <= 0.
                or not 4 <= self.max_cells <= 50000):
            raise ValueError('Graded separator needs an even counterclockwise physical rim')
        count = len(boundary)
        edges = np.roll(boundary, -1, axis=0) - boundary
        lengths = np.linalg.norm(edges, axis=1)
        if lengths.min() <= 1e-8:
            raise ValueError('Graded separator has a degenerate physical edge')
        previous = np.roll(lengths, 1)
        jumps = np.maximum(previous, lengths) > 8. * np.minimum(previous, lengths)
        if not jumps.any():
            raise ValueError('Graded separator has no abrupt boundary spacing change')
        remaining = self.max_cells - count - int(jumps.sum())
        if remaining < 4:
            raise ValueError('Graded separator row exceeds its cell budget')
        directions = edges / lengths[:, None]
        tangent = directions + np.roll(directions, 1, axis=0)
        norms = np.linalg.norm(tangent, axis=1)
        if norms.min() <= 1e-8:
            raise ValueError('Graded separator has a reversing boundary')
        tangent /= norms[:, None]
        inward = np.column_stack((-tangent[:, 1], tangent[:, 0]))
        points = list(boundary)
        faces, first, last = [], [], []
        for i, point in enumerate(boundary):
            scale = np.sqrt(previous[i] * lengths[i])
            start = len(points)
            first.append(start)
            if jumps[i]:
                # A single inward edge shared by the two boundary faces cannot
                # grade an arbitrarily large jump. A parallelogram gives each
                # face its own inward edge. Their lengths grade more slowly
                # than the prescribed curb lengths.
                sizes = np.array((previous[i], lengths[i]))
                radii = .6 * scale * (sizes / scale) ** .35
                along = np.minimum(.65 * radii, .25 * sizes)
                across = np.sqrt(radii**2 - along**2)
                left = across[0] * inward[i] - along[0] * tangent[i]
                right = across[1] * inward[i] + along[1] * tangent[i]
                points.extend((point + left, point + left + right, point + right))
                faces.append((i, start + 2, start + 1, start))
            else:
                points.append(point + .6 * scale * inward[i])
            last.append(len(points) - 1)
        for i in range(count):
            following = (i + 1) % count
            faces.append((i, following, first[following], last[i]))
        points = np.asarray(points)
        collar = Mesh(np.column_stack((points, np.zeros(len(points)))), np.asarray(faces))
        if not _qualified(collar):
            raise ValueError('Graded separator row failed physical or stored quad quality')
        inner = BandedSeparatorChart(self.patterns, spacing=self.spacing, max_cells=remaining,
                                     bands=self.bands).build(points[count:])
        # The core stores its boundary first, in the order of the row's inner
        # vertices. A constant index offset welds every shared edge.
        xyz = np.vstack((collar.xyz, inner.mesh.xyz[len(points)-count:]))
        mesh = Mesh(xyz, np.vstack((collar.quads, inner.mesh.quads + count)), 'grass')
        mesh.seams = {0: np.arange(count, dtype=np.int32)}
        certificate = audit(mesh, boundary)
        validate_overlaps({Key('grass_chart', 0): mesh})
        return ChartResult(mesh, dict(construction='graded fixed-rim separator', bands=self.bands,
            spacing_jumps=int(jumps.sum()), exact_whole_boundary_edges=True,
            rim_certificate=certificate, core=inner.metadata))
