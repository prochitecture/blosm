"""Explicit ownership across UI turns; synchronous legacy allocation capture.

capture() MUST enclose only a synchronous call. Never suspend/yield inside it.
New producers can register their allocations with own() immediately instead.
"""
from contextlib import contextmanager
import bpy

PENDING = '_blosm_pending_import'
TRANSIENT = '_blosm_transient_import'

KINDS = ('scenes', 'objects', 'collections', 'meshes', 'curves', 'cameras', 'lights',
         'worlds', 'texts', 'materials', 'images', 'node_groups', 'particles')


def alive(item):
    try:
        return bool(item.as_pointer())
    except ReferenceError:
        return False


class Ownership:
    def __init__(self):
        self.pending = {}

    def own(self, kind, item):
        self.pending[item.as_pointer()] = (kind, item)
        item[PENDING] = True
        return item

    @contextmanager
    def capture(self):
        before = {kind: {item.as_pointer() for item in getattr(bpy.data, kind)} for kind in KINDS}
        allocations = []
        try:
            yield allocations
        finally:
            for kind in KINDS:
                for item in getattr(bpy.data, kind):
                    if item.as_pointer() not in before[kind]:
                        self.own(kind, item)
                        allocations.append((kind, item))

    def accept(self, allocations=None):
        if allocations is None:
            for _, item in self.pending.values():
                if alive(item) and PENDING in item:
                    del item[PENDING]
            self.pending.clear()
        else:
            for _, item in allocations:
                if alive(item):
                    if PENDING in item:
                        del item[PENDING]
                    self.pending.pop(item.as_pointer(), None)

    def rollback(self):
        entries = list(self.pending.values())
        self.pending.clear()
        # Remove consumers first, then revisit their dependencies. A single
        # pass leaks child node groups and persistence Texts when a later
        # removal releases their final user. Never purge unrelated orphans.
        while entries:
            retained = []
            for kind, item in entries:
                if alive(item):
                    if (kind not in ('objects', 'collections', 'scenes') and
                            item.users > int(item.use_fake_user)):
                        retained.append((kind, item))
                    else:
                        getattr(bpy.data, kind).remove(item, do_unlink=True)
            if len(retained) == len(entries):
                break
            entries = retained
        # Remaining assets are shared with retained results or user work.
        for _, item in entries:
            if PENDING in item:
                del item[PENDING]
            if TRANSIENT in item:
                del item[TRANSIENT]


@bpy.app.handlers.persistent
def discard_pending_history(*_):
    """A user undo/save snapshot may include hidden staged IDs, never accepted output.

Remove only explicitly tagged allocations on restore. Committing strips tags;
canonical objects and user-created IDs are never inferred from global deltas.
"""
    ownership = Ownership()
    for kind in KINDS:
        for item in getattr(bpy.data, kind):
            if item.get(PENDING) or item.get(TRANSIENT):
                ownership.pending[item.as_pointer()] = (kind, item)
    ownership.rollback()
