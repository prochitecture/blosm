"""Inspect independent construction failures in captured numerical inputs.

No Blender data or mesh is published. Failed products remain absent, and their
dependent tasks are skipped. This is an inventory of independent failures, never
a way to bypass the normal transactional import or qualify partial geometry.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import time
import traceback

import _python_street_package
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.model import Network, fingerprint
from tested_python_streets.updates import BuildContext, DependencyIndex


def inspect(network):
    backend = PythonMeshBackend(network.settings)
    graph = DependencyIndex(tuple(t for p in backend.providers for t in p.tasks(network)))
    # Include source preparation only; a failed independent path must not hide
    # every other path/contact failure behind the first global ownership task.
    wanted = {k for k in graph.tasks if k.kind in ('path', 'profile', 'physical_intersection')}
    pending = list(wanted)
    while pending:
        key = pending.pop()
        for dependency in graph.tasks[key].dependencies:
            if dependency in graph.tasks and dependency not in wanted:
                wanted.add(dependency)
                pending.append(dependency)
    sources = {k: fingerprint(v) for k, v in network.records().items()}
    products, digests, failed, skipped, events = {}, {}, set(), set(), []
    context = BuildContext(network, products, backend.patterns, digests, sources)
    for key in graph.order:
        if key not in wanted:
            continue
        task = graph.tasks[key]
        unavailable = set(task.dependencies) & (failed | skipped)
        if unavailable:
            skipped.add(key)
            continue
        started = time.perf_counter()
        try:
            dependencies = [(d.text(), sources.get(d, digests.get(d))) for d in sorted(task.dependencies)]
            if task.signature is not None:
                dependencies = ([d.text() for d in sorted(task.dependencies)], fingerprint(task.signature(context)))
            value = task.build(context)
            products[key] = value
            digests[key] = fingerprint((task.version, dependencies))
        except Exception as exc:
            failed.add(key)
            events.append(dict(task=key.text(), error=str(exc), traceback=traceback.format_exc()))
            print(network.uuid, key.text(), str(exc), flush=True)
        finally:
            if key in failed:
                events[-1]['seconds'] = time.perf_counter() - started
    return dict(network=network.uuid, sections=len(network.sections),
        built=len(products), failed=len(failed), skipped=len(skipped), failures=events,
        failure_kinds=dict(Counter(k.kind for k in failed)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding='utf-8'))
    networks = [Network.from_dict(row) for row in (data if isinstance(data, list) else [data])]
    reports = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for network in networks:
        reports.append(inspect(network))
        args.output.write_text(json.dumps(dict(qualified=False, published=False, reports=reports), indent=2),
                               encoding='utf-8')


if __name__ == '__main__':
    main()
