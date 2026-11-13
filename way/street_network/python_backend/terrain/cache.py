"""Bounded reuse of terrain constraints, never fitted output coordinates.

Each entry contains local control floors and terrain heights at source vertices.
Every fit still rebuilds contacts, propagation, curves, clearance/lift validation
and output positions. A lower terrain or removed constraint can lower the fit.
"""
from collections import OrderedDict
from dataclasses import dataclass
import hashlib
import struct
import numpy as np


@dataclass(frozen=True)
class TerrainConstraints:
    floors: np.ndarray
    ground: np.ndarray
    samples: int


class TerrainConstraintCache:
    def __init__(self, max_bytes=8 * 1024 * 1024, max_entries=512):
        if type(max_bytes) is not int or max_bytes < 0 or type(max_entries) is not int or max_entries < 0:
            raise ValueError('Invalid terrain constraint cache budget')
        self.max_bytes, self.max_entries = max_bytes, max_entries
        self.bytes = self.hits = self.misses = self.evictions = 0
        self.entries = OrderedDict()

    @staticmethod
    def key(component, reference, settings, terrain_token):
        # No token means that an arbitrary callable has no verifiable terrain
        # revision. Native Blender supplies its evaluated-geometry signature.
        if terrain_token is None:
            return None
        digest = hashlib.sha256()
        for text in (str(terrain_token), component.identity, component.reference):
            value = text.encode('utf8')
            digest.update(struct.pack('<Q', len(value)) + value)
        digest.update(struct.pack('<dddd', component.datum, component.clearance,
                                  settings.clearance, settings.sample_spacing))
        for value in (component.xyz, component.quads, component.stations, reference.stations):
            digest.update(str((value.dtype.str, value.shape)).encode('ascii'))
            digest.update(value.tobytes())
        return digest.digest()

    def get(self, key):
        if key is None:
            return None
        saved = self.entries.get(key)
        if saved is None:
            self.misses += 1
            return None
        self.hits += 1
        self.entries.move_to_end(key)
        return saved

    @staticmethod
    def size(entry):
        # Include a conservative per-entry allowance for Python metadata.
        return entry.floors.nbytes + entry.ground.nbytes + 512

    def put(self, key, floors, ground, samples):
        size = floors.nbytes + ground.nbytes + 512
        if key is None or size > self.max_bytes or self.max_entries == 0:
            return
        previous = self.entries.pop(key, None)
        if previous is not None:
            self.bytes -= self.size(previous)
        while self.entries and (self.bytes + size > self.max_bytes or len(self.entries) >= self.max_entries):
            _, old = self.entries.popitem(last=False)
            self.bytes -= self.size(old)
            self.evictions += 1
        floors, ground = floors.copy(), ground.copy()
        floors.flags.writeable = ground.flags.writeable = False
        self.entries[key] = TerrainConstraints(floors, ground, int(samples))
        self.bytes += size

    def stats(self):
        return dict(bytes=self.bytes, limit=self.max_bytes, entries=len(self.entries),
                    hits=self.hits, misses=self.misses, evictions=self.evictions)
