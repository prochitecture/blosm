"""Compare candidate UV output geometry with a specified pre-UV Git revision.

Host Python + NumPy; copied reference sources and results are local artifacts.
No Blender or saved scene is required. Never changes the selected Git checkout.
"""
import argparse
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import time
import numpy as np
import _python_street_package
import tested_python_streets as candidate
from tested_python_streets.backend import PythonMeshBackend
from _python_street_fixtures import fixture_networks, transition_fixture_networks, short_junction_pair


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', default='e67035c')
    parser.add_argument('--output', type=Path,
        default=root.parent / 'artifacts/streets/uv_mapping/geometry_identity')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    reference = args.output / 'reference_source'
    reference.mkdir(exist_ok=True)
    archive = subprocess.check_output(['git', 'archive', args.reference,
        'way/street_network/python_backend'], cwd=root)
    with tarfile.open(fileobj=io.BytesIO(archive)) as source:
        for entry in source.getmembers():
            target = (reference / entry.name).resolve()
            if not target.is_relative_to(reference.resolve()) or entry.issym() or entry.islnk():
                raise ValueError('Unexpected reference archive entry')
        source.extractall(reference)
    package = reference / 'way/street_network/python_backend'
    spec = importlib.util.spec_from_file_location('_street_uv_reference', package / '__init__.py',
        submodule_search_locations=[str(package)])
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    from _street_uv_reference import model as reference_model
    from _street_uv_reference.backend import PythonMeshBackend as ReferenceBackend
    from tested_python_streets import model

    def networks(m):
        return (*fixture_networks(m), *transition_fixture_networks(m)[:2], short_junction_pair(m))

    report = []
    for old, new in zip(networks(reference_model), networks(model)):
        start = time.perf_counter()
        before = ReferenceBackend(old.settings).prepare(old)
        after = PythonMeshBackend(new.settings).prepare(new)
        meshes = lambda update: {k.text(): p.mesh for k, p in update.products.items()
                                 if getattr(p, 'mesh', None) is not None}
        a, b = meshes(before), meshes(after)
        assert a.keys() == b.keys()
        for key in a:
            np.testing.assert_array_equal(a[key].xyz, b[key].xyz, err_msg=key)
            np.testing.assert_array_equal(a[key].quads, b[key].quads, err_msg=key)
            assert a[key].material == b[key].material
            assert 'UVMap' in b[key].uv_layers
        result = dict(fixture=new.uuid, meshes=len(a), quads=sum(len(m.quads) for m in a.values()),
                      exact_positions=True, exact_topology=True, seconds=time.perf_counter()-start)
        report.append(result)
        print('UV_GEOMETRY_IDENTITY', json.dumps(result), flush=True)
    (args.output / 'checks.json').write_text(json.dumps(dict(reference=args.reference,
        candidate=candidate.VERSION, passed=True, fixtures=report), indent=2), encoding='utf8')


if __name__ == '__main__':
    main()
