"""Authored Bundle regression families retained from the GN numerical oracle."""
import numpy as np


def fixture_records():
    x = np.arange(7) * 12.
    bend = np.array((0., .4, 1.4, 2.2, 1.4, .4, 0.))

    def record(name, distances, *, widths=(6., 6.), offsets=(0., 0.), curved=True, reverse=False, interpolation=0, islands=None, valid=True):
        distances = np.broadcast_to(np.asarray(distances), x.shape)
        center = bend if curved else np.zeros_like(x)
        paths = []
        for sign in (1., -1.):
            points = np.column_stack((x, center + .5 * distances * sign, np.zeros_like(x)))
            first, last = points[0].copy(), points[-1].copy()
            first[0] -= 12.
            last[0] += 12.
            paths.append(np.vstack((first, points, last)))
        if reverse:
            paths[1] = paths[1][::-1]
        return dict(name=name, points_a=paths[0], points_b=paths[1], widths=widths,
                    offsets=offsets,
                    reverse=reverse, interpolation=interpolation, islands=islands, valid=valid)

    return (
        record('curved_separated', 9., islands=1),
        record('curved_overlap', 5., islands=0),
        record('three_islands', (5., 9., 5., 9., 5., 9., 5.), islands=3),
        record('unequal_widths', (5., 9., 5., 9., 5., 9., 5.), widths=(5., 7.), islands=3),
        record('reversed_member', (5., 9., 5., 9., 5., 9., 5.), reverse=True, islands=3),
        record('smooth_islands', (5., 9., 5., 9., 5., 9., 5.), interpolation=3, islands=3),
        record('contained_member', 4., widths=(4., 14.), curved=False, islands=0),
        record('exact_line_contact', 6., curved=False, islands=0),
        record('near_contact', 6.0002, curved=False),
        record('isolated_contact_probe', (8., 7., 6., 7., 8., 8., 8.), curved=False),
        record('offset_reordered', 6., curved=False, offsets=(20., 0.), islands=1),
        record('centerline_contact', (9., 7., 3., 0., 3., 7., 9.), curved=False, valid=False),
        record('side_order_reversal', (9., 7., 3., -2., 3., 7., 9.), curved=False, valid=False),
    )


def make_network(model, fixture, *, interpolation=None):
    mode = fixture['interpolation'] if interpolation is None else interpolation
    sections = (model.Section(1, 1, fixture['points_a'], width=fixture['widths'][0], offset=fixture['offsets'][0], rise=.12, interpolation=mode),
                model.Section(2, 2, fixture['points_b'], width=fixture['widths'][1], offset=fixture['offsets'][1], rise=.12, interpolation=mode))
    bundle = model.Bundle(7, (model.BundleMember((model.Traversal(1),)),
                             model.BundleMember((model.Traversal(2, not fixture['reverse']),))))
    return model.Network('Bundle / ' + fixture['name'], sections, bundles=(bundle,))


def bundle_networks(model, *, interpolation=None):
    return tuple(make_network(model, f, interpolation=interpolation) for f in fixture_records() if f['valid'])
