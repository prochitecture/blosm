"""Optional chart-attempt ledger and bounded kernel replay inputs for benchmarks.

Nested call times are cumulative, not additive. This instrumentation never
changes construction decisions and must not be included in timing baselines.
"""
from contextlib import ExitStack
import importlib
import json
import sys
import time
from unittest.mock import patch

import numpy as np


class ConstructionDiagnostics:
    def __init__(self, output, replay_bytes=16 * 1024 * 1024):
        self.output, self.budget = output, replay_bytes
        self.records, self.polygons = [], []
        self.separators, self.quality_meshes, self.active = [], [], []
        self.bytes = 0
        self.stack = ExitStack()

    def wrapper(self, original, name):
        def call(*args, **kwargs):
            record = dict(operation=name, parent=self.active[-1] if self.active else None)
            if name == 'PatternLibrary.quad':
                record['counts'] = list(args[1:5])
            if name == 'visibility_center':
                polygon = np.asarray(args[0])
                record['vertices'] = len(polygon)
                record['epsilon'] = kwargs.get('epsilon', args[1] if len(args) > 1 else 1e-10)
                if self.bytes + polygon.nbytes <= self.budget:
                    record['replay'] = len(self.polygons)
                    self.polygons.append(polygon.copy())
                    self.bytes += polygon.nbytes
            else:
                for value in args[1:]:
                    if isinstance(value, np.ndarray):
                        record['input_shape'] = list(value.shape)
                        break
                    if hasattr(value, 'quads'):
                        record['input_quads'] = len(value.quads)
                        break
            if name == 'FixedRimSeparator.build':
                boundary = np.asarray(args[1])
                if self.bytes + boundary.nbytes <= self.budget:
                    record['separator_replay'] = len(self.separators)
                    record['settings'] = dict(spacing=args[0].spacing, max_cells=args[0].max_cells,
                                              **kwargs)
                    self.separators.append(boundary.copy())
                    self.bytes += boundary.nbytes
            if name == 'FixedRimQuality.improve':
                mesh = args[1]
                size = mesh.xyz.nbytes + mesh.quads.nbytes
                if len(mesh.quads) >= 128 and self.bytes + size <= self.budget:
                    record['quality_replay'] = len(self.quality_meshes)
                    record['quality_settings'] = dict(iterations=args[0].iterations,
                        target_ratio=args[0].target_ratio, power=args[0].power)
                    self.quality_meshes.append((mesh.xyz.copy(), mesh.quads.copy()))
                    self.bytes += size
            # Append before calling so each event has a stable preorder index.
            self.active.append(len(self.records))
            self.records.append(record)
            start = time.perf_counter()
            try:
                result = original(*args, **kwargs)
                record['status'] = 'accepted'
                metadata = getattr(result, 'metadata', None)
                if isinstance(result, tuple) and len(result) == 2 and isinstance(result[1], dict):
                    metadata = result[1]
                if isinstance(metadata, dict):
                    record['iterations'] = {k: metadata[k] for k in ('steps', 'iterations', 'stop') if k in metadata}
                    if name == 'FixedRimSeparator.build':
                        record['result'] = metadata
                return result
            except Exception as error:
                record.update(status='rejected', error=str(error))
                raise
            finally:
                record['seconds'] = time.perf_counter() - start
                self.active.pop()
        return call

    def install(self, package):
        prefix = package + '.way.street_network.python_backend.'
        targets = [('patterns', 'PatternLibrary', 'quad'),
                   ('quad_charts', 'LayeredQuadChart', 'build'),
                   ('quad_charts', 'FixedRimQuality', 'improve'),
                   ('decomposition', 'VisibilitySplitter', 'split'),
                   ('four_side_charts', 'FourSidedQuadChart', 'build'),
                   ('separator_surfaces', 'FixedRimSeparator', 'build'),
                   ('partition', 'CompoundQuadChart', 'build'),
                   ('slope_projection', 'SlopeProjection', 'solve')]
        for module, cls, method in targets:
            owner = getattr(importlib.import_module(prefix + module), cls)
            self.stack.enter_context(patch.object(owner, method,
                self.wrapper(getattr(owner, method), cls + '.' + method)))
        module = importlib.import_module(prefix + 'quad_charts')
        original = module.visibility_center
        wrapped = self.wrapper(original, 'visibility_center')
        for module in tuple(sys.modules.values()):
            if (getattr(module, '__name__', '').startswith(prefix)
                    and getattr(module, 'visibility_center', None) is original):
                self.stack.enter_context(patch.object(module, 'visibility_center', wrapped))

    def close(self):
        self.stack.close()
        self.output.mkdir(parents=True, exist_ok=True)
        (self.output / 'construction.json').write_text(json.dumps(dict(records=self.records,
            replay_bytes=self.bytes, replay_budget=self.budget,
            qualification='Instrumented nested call times; not a timing baseline'), indent=2), encoding='utf8')
        np.savez_compressed(self.output / 'kernel_inputs.npz',
                            **{'p' + str(i): value for i, value in enumerate(self.polygons)})
        np.savez_compressed(self.output / 'separator_inputs.npz',
                            **{'p' + str(i): value for i, value in enumerate(self.separators)})
        np.savez_compressed(self.output / 'quality_inputs.npz',
            **{f'q{i}_{name}': values[j] for i, values in enumerate(self.quality_meshes)
               for j, name in enumerate(('xyz', 'quads'))})
