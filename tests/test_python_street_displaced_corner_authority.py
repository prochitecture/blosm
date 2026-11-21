"""Keep logical corner authority distinct from later physical source contacts."""
import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest
import numpy as np
import _python_street_package
from test_python_street_displaced_source_contacts import fixture, context
from tested_python_streets import model
from tested_python_streets.compound_displaced_contacts import retry
from tested_python_streets.compounds import SharedClaim
from tested_python_streets.displaced_corner_authority import (
    validate_displaced_corner, original_parent_corner, DisplacedCornerAuthorityError)
from tested_python_streets.junction_bank_transitions import original_corner


def case(network=None):
    network = fixture() if network is None else network
    network = replace(network, intersections=tuple(replace(site, radius=10.)
                      if site.id == 10 else site for site in network.intersections))
    ctx, a, b = context(network)
    members = tuple(model.Key('intersection', i) for i in (10, 11))
    corner = retry(ctx, a, b, 10., 'STYLE', set(members), {1})
    if corner is None:
        raise AssertionError('The actual source-contact producer did not return a corner')
    group = SimpleNamespace(members=members,
        claims=(SharedClaim(1, *members, 0., ctx[model.Key('path', 1)].length),),
        physical_layout=SimpleNamespace(approaches=(corner.first, corner.second),
                          outward_stations=(corner.first.length, corner.second.length)))
    node = next(site for site in network.intersections if site.id == 10)
    return ctx, group, node, corner


class DisplacedCornerAuthorityTests(unittest.TestCase):
    def test_later_contact_is_valid_without_inventing_a_logical_node_incidence(self):
        ctx, group, node, corner = case()
        before = ctx.network.to_dict()
        resolver = SimpleNamespace(context=ctx, group=group)
        self.assertFalse(original_corner(resolver, node, corner))
        proof = validate_displaced_corner(ctx, group, node, corner)
        self.assertEqual(proof.logical_ends, (model.End(2, 0), model.End(1, 0)))
        self.assertEqual(proof.contact_ends, (model.End(2, 0), model.End(4, 0)))
        self.assertNotIn(proof.contact_ends[1], node.ends)
        self.assertTrue(original_parent_corner(resolver, node, corner))
        self.assertEqual(ctx.network.to_dict(), before)

    def test_reversed_storage_and_rigid_motion_keep_real_path_authority(self):
        network = fixture()
        rotation = np.array(((np.cos(.7), -np.sin(.7)), (np.sin(.7), np.cos(.7))))
        transformed = []
        for section in network.sections:
            points = np.asarray(section.points).copy()
            points[:, :2] = points[:, :2] @ rotation.T + (1200., -800.)
            transformed.append(replace(section, points=tuple(map(tuple, points))))
        reversed_network = replace(network,
            sections=tuple(replace(s, points=s.points[::-1]) for s in network.sections),
            intersections=tuple(replace(site, ends=tuple(model.End(e.section_id, 1-e.end)
                                for e in site.ends)) for site in network.intersections))
        for variant in (replace(network, sections=tuple(transformed)), reversed_network):
            ctx, group, node, corner = case(variant)
            proof = validate_displaced_corner(ctx, group, node, corner)
            self.assertEqual(proof.junction, node.key)

    def test_wrong_prefix_or_borrowed_path_cannot_supply_contact_authority(self):
        ctx, group, node, corner = case()
        resolver = SimpleNamespace(context=ctx, group=group)
        invalid = copy.copy(corner)
        first, second = corner.source_contact_walks
        invalid.source_contact_walks = (first, replace(second, prefix=()))
        self.assertFalse(original_parent_corner(resolver, node, invalid))
        arm = copy.copy(second.approach)
        arm.path = copy.copy(arm.path)
        pieces = list(arm.path.pieces)
        pieces[0] = replace(pieces[0], path=copy.copy(pieces[0].path))
        arm.path.pieces = tuple(pieces)
        invalid = copy.copy(corner)
        invalid.second = arm
        invalid.source_contact_walks = (first, replace(second, approach=arm))
        self.assertFalse(original_parent_corner(resolver, node, invalid))

    def test_unowned_prefix_and_contact_are_not_accepted_as_original_geometry(self):
        ctx, group, node, corner = case()
        unowned = copy.copy(group)
        unowned.claims = ()
        self.assertFalse(original_parent_corner(SimpleNamespace(context=ctx, group=unowned), node, corner))
        truncated = copy.copy(group)
        truncated.physical_layout = SimpleNamespace(approaches=(corner.first, corner.second),
            outward_stations=(corner.stations[0], corner.stations[1]-.001))
        self.assertFalse(original_parent_corner(SimpleNamespace(context=ctx, group=truncated), node, corner))

    def test_radius_policy_and_actual_geometry_cannot_be_relabelled(self):
        ctx, group, node, corner = case()
        resolver = SimpleNamespace(context=ctx, group=group)
        for change in (dict(radius_mode='STRICT'), dict(requested_radius=9.),
                       dict(radius_attempts=(corner.radius,)),
                       dict(center=corner.center+(.001, 0.)),
                       dict(start_angle=corner.start_angle+.001),
                       dict(stations=(corner.first.length+.001, corner.stations[1]))):
            invalid = copy.copy(corner)
            for key, value in change.items():
                setattr(invalid, key, value)
            with self.subTest(change=tuple(change)):
                self.assertFalse(original_parent_corner(resolver, node, invalid))

    def test_real_original_node_and_group_membership_remain_required(self):
        ctx, group, node, corner = case()
        copied_node = replace(node)
        with self.assertRaises(DisplacedCornerAuthorityError):
            validate_displaced_corner(ctx, group, copied_node, corner)
        outside = copy.copy(group)
        outside.members = (model.Key('intersection', 11),)
        self.assertFalse(original_parent_corner(SimpleNamespace(context=ctx, group=outside), node, corner))


if __name__ == '__main__':
    unittest.main()
