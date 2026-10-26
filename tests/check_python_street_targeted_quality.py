"""Target actual quad-shape limits on an assembled real junction chart."""
import argparse
import json
from pathlib import Path
import time
from unittest.mock import patch
import numpy as np
import _python_street_package
from _python_targeted_quality import TargetedQuality
from check_python_street_chart_failure import context_for
from tested_python_streets import model
from tested_python_streets.compounds import CompoundGroup, sample_footprint
from tested_python_streets.ownership import footprint_key
from tested_python_streets.partition import CompoundQuadChart
from tested_python_streets.regions import RegionUnion
from tested_python_streets.quad_rotations import QuadEdgeRotations
from tested_python_streets.quad_charts import FixedRimQuality
from tested_python_streets.validation import validate_overlaps
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.providers import Part


class TargetedChart(CompoundQuadChart):
    def build(self, group, context, *, require_quality=True):
        try:
            return super().build(group, context, require_quality=require_quality)
        except ValueError as exc:
            if not require_quality or 'needs better grading' not in str(exc):
                raise
            original_error = str(exc)
        chart = super().build(group, context, require_quality=False)
        seed, rotations = QuadEdgeRotations(max_passes=8, max_cells=self.max_cells).improve(
            chart.mesh, FixedRimQuality(iterations=300, power=2))
        mesh, report = TargetedQuality().improve(seed)
        quality = mesh.validate()
        if quality['maximum_edge_ratio'] > 8. or report['final']['minimum_sine'] < .03:
            raise ValueError('Targeted planar correction did not meet the unchanged quality limits')
        chart.mesh = mesh
        chart.metadata = dict(chart.metadata, quality=quality, targeted_quality=report,
                              seed_rotations=rotations, initial_chart_failure=original_error)
        return chart


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--junction', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--try-network', action='store_true')
    args = parser.parse_args();args.output.mkdir(parents=True, exist_ok=True)
    network, = [n for n in map(model.Network.from_dict, json.loads(args.input.read_text()))
                if any(j.id == args.junction for j in n.intersections)]
    site = next(j for j in network.intersections if j.id == args.junction)
    context = context_for(network)
    footprint = sample_footprint(site.key, context[footprint_key(site.key)])
    group = CompoundGroup((site.key,), (), (footprint,), RegionUnion().build((footprint.region,)))
    chart = CompoundQuadChart(context.patterns, spacing=network.settings.spacing).build(group, context, require_quality=False)
    mesh = chart.mesh
    initial = mesh.validate()
    mesh, rotations = QuadEdgeRotations(max_passes=8).improve(mesh, FixedRimQuality(iterations=300, power=2))
    start = time.perf_counter()
    candidate, details = TargetedQuality().improve(mesh, lambda r: print(json.dumps(r), flush=True))
    report = dict(initial=initial, seed=mesh.validate(), rotations=rotations, targeted=details,
        final=candidate.validate(), seconds=time.perf_counter() - start,
        planar_only=True, full_network_qualified=False, overlaps=validate_overlaps({site.key: candidate}, cell_size=1.4))
    np.savez_compressed(args.output / 'candidate.npz', xyz=candidate.xyz, quads=candidate.quads,
                        boundary=np.concatenate([mesh.seams[i] for i in mesh.seams]))
    (args.output / 'targeted_checks.json').write_text(json.dumps(report, indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ('rotations', 'targeted')}), flush=True)
    if args.try_network:
        start = time.perf_counter()
        try:
            with patch('tested_python_streets.compound_surfaces.CompoundQuadChart', TargetedChart):
                update = PythonMeshBackend(network.settings).prepare(network)
            parts = {key: p for key, p in update.products.items() if isinstance(p, Part) and p.mesh is not None}
            result = dict(passed=True, numerical_network_only=True, blender_publication_qualified=False,
                quads=sum(len(p.mesh.quads) for p in parts.values()), seconds=time.perf_counter() - start,
                overlap=update.timings['overlap_check'],
                parts={key.text(): dict(quality=p.mesh.validate(), metadata=p.metadata) for key, p in parts.items()})
            for key, p in parts.items():
                np.savez_compressed(args.output / (key.text().replace(':', '_') + '.npz'), xyz=p.mesh.xyz, quads=p.mesh.quads)
        except ValueError as exc:
            result = dict(passed=False, reason=str(exc), seconds=time.perf_counter() - start)
        (args.output / 'network_checks.json').write_text(json.dumps(result, indent=2))
        print(json.dumps({key: value for key, value in result.items() if key != 'parts'}), flush=True)


if __name__ == '__main__':
    main()
