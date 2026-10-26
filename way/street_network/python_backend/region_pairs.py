"""Bounded spatial candidates for exact sampled-boundary intersection tests."""
import math
import numpy as np


class PairBudget:
    def __init__(self, limit):
        self.limit, self.used = limit, 0

    def charge(self, count):
        if count > self.limit - self.used:
            raise ValueError('Compound boundary exceeds its intersection-test budget')
        self.used += count


class SegmentPairIndex:
    """Conservative edge-box grid; long edges have a bounded fallback list.

    This only removes pairs whose boxes cannot meet. Callers retain their
    original narrow intersection arithmetic and curve samples. The expanded
    boxes deliberately overselect at floating-point cell boundaries.
    """
    MAX_EDGE_CELLS = 64
    BATCH_SIZE = 4096

    def __init__(self, minimum, maximum, epsilon):
        self.minimum, self.maximum, self.epsilon = minimum, maximum, epsilon
        self.origin = minimum.min(axis=0)
        span = maximum.max(axis=0) - self.origin
        self.cell = max(4. * epsilon, float(span.max()) / math.sqrt(len(minimum)))
        if not math.isfinite(self.cell):
            raise ValueError('Compound boundary has nonfinite spatial extent')
        self.buckets, self.long_edges = {}, []
        for index, (low, high) in enumerate(zip(minimum, maximum)):
            cells = self._cells(low, high)
            if cells is None:
                self.long_edges.append(index)
            else:
                for cell in cells:
                    self.buckets.setdefault(cell, []).append(index)

    def _cells(self, low, high):
        low = np.nextafter(low - self.epsilon, -np.inf)
        high = np.nextafter(high + self.epsilon, np.inf)
        lower, upper = (low - self.origin) / self.cell, (high - self.origin) / self.cell
        if not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper)):
            return None
        x0, y0 = (math.floor(float(value)) for value in lower)
        x1, y1 = (math.floor(float(value)) for value in upper)
        if (x1 - x0 + 1) * (y1 - y0 + 1) > self.MAX_EDGE_CELLS:
            return None
        return ((x, y) for x in range(x0, x1 + 1) for y in range(y0, y1 + 1))

    def batches(self, minimum, maximum, budget, *, same=False):
        # Charge queries before allocating the grid's result scratch. Each
        # unique candidate box is charged before its NumPy comparison too.
        budget.charge(len(minimum))
        rows, columns = [], []
        for row, (low, high) in enumerate(zip(minimum, maximum)):
            cells = self._cells(low, high)
            if cells is None:
                candidates = range(len(self.minimum))
            else:
                budget.charge(len(self.long_edges))
                found = set(self.long_edges)
                for cell in cells:
                    bucket = self.buckets.get(cell, ())
                    # Count repeated bucket entries too, so a dense grid
                    # cannot hide excessive set-building work from the cap.
                    budget.charge(len(bucket))
                    found.update(bucket)
                candidates = sorted(found)
            if same:
                candidates = [column for column in candidates
                              if column > row + 1 and not (row == 0 and column == len(minimum) - 1)]
            budget.charge(len(candidates))
            if not len(candidates):
                continue
            candidates = np.asarray(candidates, dtype=np.intp)
            nearby = np.all(high + self.epsilon >= self.minimum[candidates], axis=-1)
            nearby &= np.all(self.maximum[candidates] + self.epsilon >= low, axis=-1)
            for column in candidates[nearby]:
                rows.append(row)
                columns.append(int(column))
                if len(rows) == self.BATCH_SIZE:
                    yield np.asarray(rows, dtype=np.intp), np.asarray(columns, dtype=np.intp)
                    rows, columns = [], []
        if rows:
            yield np.asarray(rows, dtype=np.intp), np.asarray(columns, dtype=np.intp)
