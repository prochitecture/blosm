"""Existing whole-curve separation proof on actual visible source intervals.

Junction-owned prefixes and wholly owned query spans are outside these visible
road interiors. Their source incidence/ownership is certified by preparation;
this proof neither invents connectivity nor compares displaced mesh samples.
"""
from .bundles import CenterlineRange, certify_disjoint
from .model import Key


def certify_visible_centerlines(context, spans_by_role):
    if len(spans_by_role) != 2 or any(not spans for spans in spans_by_role):
        raise ValueError('Gapped centerline separation needs both actual visible source roles')
    paths, inventory = [], []
    for role, spans in enumerate(spans_by_role):
        ranges, records = [], []
        for span in spans:
            if span.role != role or len(span.visits) != len(span.source_intervals):
                raise ValueError('Gapped centerline proof lost its actual original source intervals')
            for visit, row in zip(span.visits, span.source_intervals):
                if (row['section'] != visit.section_id or row['forward'] != visit.forward):
                    raise ValueError('Gapped centerline proof changed its ordered original source traversal')
                path = context[Key('path', visit.section_id)]
                first, last = row['first'], row['last']
                ranges.append(CenterlineRange.build(path, first, last))
                records.append(dict(section=visit.section_id, forward=visit.forward,
                                    first=first, last=last))
        if not ranges:
            raise ValueError('Gapped centerline proof needs a positive visible interval on each role')
        paths.append(tuple(ranges));inventory.append(records)
    checks = certify_disjoint(*paths)
    return dict(construction='existing continuous original centerline hull separation',
        proof_checks=checks, source_intervals=inventory,
        actual_junction_owned_prefixes_excluded=True,
        wholly_owned_query_sections_excluded=True, mesh_samples_used=False,
        new_source_connectivity=False)
