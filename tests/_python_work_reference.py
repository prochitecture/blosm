"""Previous work paths for controlled tests; never imported by the addon."""
import numpy as np
from types import SimpleNamespace


def uncached_approach_frame(self, station):
    station = np.asarray(station)
    clipped = np.clip(station, 0., self.path.length)
    if getattr(self.path, 'outward', False):
        p, tangent, normal = self.path.frame(clipped)
    else:
        global_station = clipped if self.end == 0 else self.path.length - clipped
        p, tangent, normal = self.path.frame(global_station)
        if self.end == 1:
            tangent, normal = -tangent, -normal
    p = p + np.minimum(station, 0.)[..., None] * tangent
    return p, tangent, normal


def sorted_list_order(tasks):
    graph = {task.key: task for task in tasks}
    reverse = {}
    for task in tasks:
        for key in task.dependencies:
            reverse.setdefault(key, set()).add(task.key)
    return sorted_list_topological_order(SimpleNamespace(tasks=graph, reverse=reverse))


def sorted_list_topological_order(self):
    pending = {key: set(t.dependencies) & self.tasks.keys() for key, t in self.tasks.items()}
    ready = sorted(key for key, deps in pending.items() if not deps)
    order = []
    while ready:
        key = ready.pop(0)
        order.append(key)
        for child in sorted(self.reverse.get(key, ())):
            pending[child].discard(key)
            if not pending[child] and child not in order and child not in ready:
                ready.append(child)
                ready.sort()
    if len(order) != len(self.tasks):
        raise ValueError('Feature dependency cycle')
    return tuple(order)
