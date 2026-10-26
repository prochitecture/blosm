"""Explicit, transactional update planning with persistent reverse dependencies.

No depsgraph handler or polling timer is installed by the Python backend.
Call prepare after reading the authored inputs, publish replacements, then
commit. Failed preparation/publication cannot advance the last-built manifest.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import time
import heapq
from .model import Key, fingerprint
from . import VERSION


@dataclass(frozen=True)
class Task:
    key: Key
    dependencies: tuple
    build: object
    version: str = '1'
    signature: object = None


class FeatureProvider(ABC):
    source_kinds = frozenset()

    @abstractmethod
    def tasks(self, network):
        """Declare outputs and every source/derived input they consume."""


class DependencyIndex:
    def __init__(self, tasks):
        self.tasks = {task.key: task for task in tasks}
        if len(self.tasks) != len(tasks):
            raise ValueError('Two feature providers claim the same output')
        self.reverse = {}
        for task in tasks:
            if len(set(task.dependencies)) != len(task.dependencies):
                raise ValueError('Repeated task dependency')
            for key in task.dependencies:
                self.reverse.setdefault(key, set()).add(task.key)
        self.order = self._topological_order()

    def _topological_order(self):
        pending = {key: len(set(t.dependencies) & self.tasks.keys()) for key, t in self.tasks.items()}
        ready = sorted(key for key, deps in pending.items() if not deps)
        heapq.heapify(ready)
        order = []
        while ready:
            key = heapq.heappop(ready)
            order.append(key)
            for child in self.reverse.get(key, ()):
                pending[child] -= 1
                if pending[child] == 0:
                    heapq.heappush(ready, child)
        if len(order) != len(self.tasks):
            raise ValueError('Feature dependency cycle')
        return tuple(order)

    def affected(self, changed):
        result, pending = set(), list(changed)
        while pending:
            key = pending.pop()
            for child in self.reverse.get(key, ()):
                if child not in result:
                    result.add(child)
                    pending.append(child)
        return result


@dataclass
class BuildContext:
    network: object
    products: dict
    patterns: object
    task_digests: dict = None
    source_digests: dict = None

    def __getitem__(self, key):
        return self.products[key]

    def digest(self, key):
        """Digest of an already evaluated, declared upstream dependency."""
        if self.source_digests is not None and key in self.source_digests:
            return self.source_digests[key]
        if self.task_digests is not None and key in self.task_digests:
            return self.task_digests[key]
        raise ValueError(f'Dependency digest is not available: {key.text()}')


@dataclass
class PreparedUpdate:
    network_id: str
    source_fingerprints: dict
    task_fingerprints: dict
    products: dict
    rebuilt: tuple
    removed: tuple
    changed: tuple
    dependencies: dict
    timings: dict = field(default_factory=dict)
    base_revision: int = 0


class UpdateCoordinator:
    def __init__(self, providers, patterns):
        self.providers = tuple(providers)
        self.patterns = patterns
        self.source_fingerprints = {}
        self.task_fingerprints = {}
        self.products = {}
        self.network_id = None
        self.revision = 0
        self.dependencies = {}

    def prepare(self, network, *, force=False):
        start = time.perf_counter()
        if self.network_id is not None and self.network_id != network.uuid:
            raise ValueError('An updater cannot reuse outputs from another network')
        source = {key: fingerprint(item) for key, item in network.records().items()}
        settings = fingerprint(network.settings)
        source[Key('settings', 0)] = settings
        changed = {key for key in source.keys() | self.source_fingerprints.keys()
                   if source.get(key) != self.source_fingerprints.get(key)}
        tasks = tuple(task for provider in self.providers for task in provider.tasks(network))
        graph = DependencyIndex(tasks)
        for task in tasks:
            missing = set(task.dependencies) - graph.tasks.keys() - source.keys()
            if missing:
                raise ValueError(f'Missing feature dependencies: {sorted(missing)}')
        products = {key: value for key, value in self.products.items() if key in graph.tasks}
        digests, rebuilt, durations = {}, [], {}
        context = BuildContext(network, products, self.patterns, digests, source)
        for key in graph.order:
            task = graph.tasks[key]
            dependencies = [(dep.text(), source[dep] if dep in source else digests[dep])
                            for dep in sorted(task.dependencies)]
            if task.signature is not None:
                # A provider can project broad source records onto exactly the
                # fields it consumes. Keep dependency identities/order in the
                # digest and manifest; a vertex move must not invalidate an
                # unchanged shared profile merely because rise and points live
                # in the same immutable Section record.
                dependencies = ([dep.text() for dep in sorted(task.dependencies)],
                                fingerprint(task.signature(context)))
            digest = fingerprint((VERSION, task.version, dependencies, settings))
            digests[key] = digest
            if force or digest != self.task_fingerprints.get(key) or key not in products:
                before = time.perf_counter()
                try:
                    products[key] = task.build(context)
                except ValueError as exc:
                    raise ValueError(f'{key.text()}: {exc}') from exc
                durations[key.text()] = time.perf_counter() - before
                rebuilt.append(key)
        return PreparedUpdate(network.uuid, source, digests, products, tuple(rebuilt),
            tuple(sorted(self.task_fingerprints.keys() - digests.keys())), tuple(sorted(changed)),
            {key: task.dependencies for key, task in graph.tasks.items()},
            dict(total=time.perf_counter() - start, tasks=durations,
                 pattern_cache=self.patterns.stats()), self.revision)

    def commit(self, update):
        if update.base_revision != self.revision:
            raise ValueError('An update was superseded; prepare it again')
        if self.network_id is not None and update.network_id != self.network_id:
            raise ValueError('Mismatched update ownership')
        self.network_id = update.network_id
        self.source_fingerprints = dict(update.source_fingerprints)
        self.task_fingerprints = dict(update.task_fingerprints)
        self.dependencies = dict(update.dependencies)
        self.products = dict(update.products)
        self.revision += 1

    def manifest(self):
        return dict(version=1, implementation=VERSION, network_uuid=self.network_id, revision=self.revision,
            source_fingerprints={k.text(): v for k, v in self.source_fingerprints.items()},
            task_fingerprints={k.text(): v for k, v in self.task_fingerprints.items()},
            dependencies={k.text(): [d.text() for d in deps] for k, deps in self.dependencies.items()})

    def restore_manifest(self, manifest):
        """Restore saved identities; derived numerical data is rebuilt on demand.

        The scene adapter hydrates saved output metadata/identity separately.
        Missing numerical paths/layouts are reconstructed without forcing
        unchanged published meshes to be replaced.
        """
        if manifest.get('version') != 1:
            raise ValueError('Unsupported Python generation manifest')
        self.network_id = manifest['network_uuid']
        self.revision = manifest['revision']
        self.source_fingerprints = {Key.parse(k): v for k, v in manifest['source_fingerprints'].items()}
        self.task_fingerprints = {Key.parse(k): v for k, v in manifest['task_fingerprints'].items()}
        self.dependencies = {Key.parse(k): tuple(Key.parse(d) for d in deps)
                             for k, deps in manifest['dependencies'].items()}
        self.products = {}
