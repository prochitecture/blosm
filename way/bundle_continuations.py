"""Deterministic, source-backed continuation choices for adjoining Bundles."""

def bundle_merge_pairs(involvedBundles):
    """Select physical continuations without relying on equal-position order.

    A Bundle can contain several head entries at one branching source node.
    Treating the first/last entry as a lateral member then chooses a different
    road when the input set happens to iterate differently. At such a node,
    continue a unique original OSM way. Keep the other branches as original
    attached Streets; ambiguous branches must retain their separate Bundles.
    """
    records = list(involvedBundles.values())
    if len(records) != 2:
        return None
    sides = [{row['end'] for row in rows} for rows in records]
    if sides[0] != sides[1] or len(sides[0]) != 2:
        return None
    result = []
    for point in sorted(sides[0], key=lambda p: tuple(p)):
        candidates = []
        for a in {row['street'] for row in records[0] if row['end'] == point}:
            for b in {row['street'] for row in records[1] if row['end'] == point}:
                if a.dst == point and b.src == point:
                    candidates.append((a, b))
                elif b.dst == point and a.src == point:
                    candidates.append((b, a))
        if len(candidates) > 1:
            candidates = [(a, b) for a, b in candidates
                          if set(a.tail.sourceWayIds).intersection(b.head.sourceWayIds)]
        if len(candidates) != 1:
            return None
        result.append((point, *candidates[0]))
    if len({street for _, a, b in result for street in (a, b)}) != 4:
        return None
    return tuple(result)

