"""A generated connector arc uses its actual sealed native contact profiles."""
from copy import copy
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from _python_street_owned_connector_fixture import fixture
from tested_python_streets.model import Key
from tested_python_streets.compound_cycles import _port_footprint
from tested_python_streets.source_owned_connector_fillets import required,insert,bind_contacts
from tested_python_streets.source_curb_conditions import (SourceCurbConditions,
    SourceCurbConditionUnavailable,_native_gradient)


def contact_fixture(reverse=False,angle=0.):
    """Finite source-contact fixture, not a complete Product qualification."""
    planner,group,grids=fixture(reverse=reverse,angle=angle)
    record=required(planner,group,());added=[]
    insert(planner,group.members,{12,13},added,{},(record,),.7,.00025)
    layout=copy(group.physical_layout)
    actual={(arm.source.id,arm.end):station+.25
            for arm,station in zip((record.corner.first,record.corner.second),record.corner.stations)}
    layout.outward_stations=tuple(max(old,actual[arm.source.id,arm.end])
        for arm,old in zip(layout.approaches,layout.outward_stations))
    # Regenerate real finite port prefixes at these actual contacts. The
    # profiles, original logical corners and consumed source graph stay exact.
    ports=[]
    for arm,station in zip(layout.approaches,layout.outward_stations):
        ports.append(_port_footprint(Key('cycle_port',2*arm.source.id+arm.end),
            arm,station,8,.7,.00025,grids))
    group=replace(group,physical_layout=layout,footprints=tuple(
        source for source in group.footprints if source.region.owner.kind!='cycle_port')+tuple(ports)+tuple(added))
    group=bind_contacts(planner.context,group)
    source=next(source for source in group.footprints
                if getattr(source.layout,'source_owned_connector_fillet',None))
    return planner.context,group,source


class OwnedConnectorConditionTests(unittest.TestCase):
    def test_sealed_connector_dependencies_and_metadata_survive_serialization(self):
        import json
        from tested_python_streets.model import fingerprint
        from tested_python_streets.source_owned_connector_fillets import signature
        _,group,_ = contact_fixture()
        value = signature(group)
        self.assertEqual(fingerprint(value),fingerprint(json.loads(json.dumps(value))))
        self.assertTrue(json.loads(json.dumps(group.physical_layout.source_owned_connector_reconciliation))['records'])

    def test_actual_wide_contact_is_not_the_logical_narrow_connector_profile(self):
        for reverse,angle in ((False,0.),(True,.713)):
            context,group,source=contact_fixture(reverse,angle)
            record=source.layout.source_owned_connector_fillet
            label=next(value for value in source.segments if value.kind=='curb')
            rule=SourceCurbConditions(context,group).source_rule(source.region.owner,label.index,0.,1.)
            self.assertTrue(rule.corner_fragment);self.assertFalse(rule.source_anchor)
            self.assertEqual(record.original_corner.last.path.section.id,12)
            self.assertEqual(record.contacts[1].path.section.id,14)
            self.assertNotEqual(record.original_corner.last.path.section.width,record.contacts[1].path.section.width)
            heights,_=rule.sample(np.array((0.,1.)))
            expected=[]
            for contact in record.contacts:
                profile=context[Key('profile',contact.path.section.id)]
                expected.append(profile.datum+float(profile.value(contact.side)))
                self.assertTrue(contact.owned_intervals)
                self.assertTrue(any(first<=contact.station<=last for first,last in contact.owned_intervals))
            # The prepared arc may be stored in the opposite boundary order.
            self.assertEqual(sorted(heights),sorted(expected))
            expected_slopes=sorted(float(np.linalg.norm(_native_gradient(
                contact.path,context[Key('profile',contact.path.section.id)],contact.station,contact.side)))
                for contact in record.contacts)
            np.testing.assert_allclose(sorted((rule.k0,rule.k1)),expected_slopes,rtol=0.,atol=1e-14)

    def test_unsealed_wrong_side_and_changed_native_contacts_are_refused(self):
        context,group,source=contact_fixture()
        record=source.layout.source_owned_connector_fillet
        label=next(value for value in source.segments if value.kind=='curb')
        first=record.contacts[0]
        for contact in (replace(first,owned_intervals=()),replace(first,side=-first.side),
                        replace(first,station=first.station+.01)):
            changed=replace(record,contacts=(contact,record.contacts[1]))
            layout=copy(source.layout);layout.generated_junction_area=changed
            layout.source_owned_connector_fillet=changed
            replacement=replace(source,layout=layout)
            current=replace(group,footprints=tuple(replacement if value is source else value
                                                  for value in group.footprints))
            with self.assertRaisesRegex(SourceCurbConditionUnavailable,'sealed finite'):
                SourceCurbConditions(context,current).source_rule(source.region.owner,label.index,0.,1.)


if __name__=='__main__':unittest.main()
