"""Semantic cross-link classification retains original input and incidences."""
from copy import deepcopy
from pathlib import Path
import importlib.util
import sys
from types import SimpleNamespace as NS
import unittest

path=Path(__file__).resolve().parents[1]/'way/street_network/junction_interiors.py'
spec=importlib.util.spec_from_file_location('tested_junction_interiors',path)
review=importlib.util.module_from_spec(spec);sys.modules[spec.name]=review;spec.loader.exec_module(review)


def sample(reverse=False):
    points=((-6.,0.,0.),(6.,0.,0.))
    if reverse:points=points[::-1]
    source=dict(id=30,endpoints=(100,100),status='SUPPORTED',items=[dict(kind='Section',id=300,
        points=points,width=4.,offset=0.,category='residential',source_way_ids=['original'])])
    snapshot=NS(streets={30:source});intersection=NS(id=100)
    bundle=NS(id=7,streetsHead=[1,2],streetsTail=[1,2],headLocs=((-6.,0.),(6.,0.)),tailLocs=((-6.,100.),(6.,100.)))
    street=NS(id=30)
    connectors=[NS(id=10,item=bundle,leaving=True),NS(id=12,item=street,leaving=False),NS(id=11,item=street,leaving=True)]
    return snapshot,intersection,connectors


class InteriorLinkTests(unittest.TestCase):
    def setUp(self):self.review=review.PreparedJunctionInteriorReview()

    def test_two_authoring_directions_and_source_preservation(self):
        for reverse in (False,True):
            s,j,c=sample(reverse);before=deepcopy(s)
            links,retained=self.review.read(s,j,c)
            self.assertEqual(s,before);self.assertFalse(retained);self.assertEqual(len(links),1)
            link=links[0]
            self.assertEqual(link.connector_ids,(11,12));self.assertEqual(link.section_ids,(300,))
            self.assertEqual(link.endpoint_members,(1,0) if reverse else (0,1))
            self.assertEqual(link.report()['status'],'PENDING_FOOTPRINT')

    def test_multisection_street_retains_every_section(self):
        s,j,c=sample();source=s.streets[30]['items'][0]
        s.streets[30]['items']=[dict(source,points=((-6.,0.,0.),(0.,0.,0.))),
                              dict(source,id=301,points=((0.,0.,0.),(6.,0.,0.)))]
        links,_=self.review.read(s,j,c);self.assertEqual(links[0].section_ids,(300,301))

    def test_true_returning_path_is_not_silently_absorbed(self):
        s,j,c=sample();s.streets[30]['items'][0]['points']=((-6.,0.,0.),(-6.,-30.,0.),(6.,-30.,0.),(6.,0.,0.))
        links,retained=self.review.read(s,j,c)
        self.assertFalse(links);self.assertEqual(retained[0]['status'],'RETAINED_APPROACH')

    def test_wrong_endpoint_or_duplicate_incidence_retained(self):
        for change in ('external','duplicate'):
            s,j,c=sample()
            if change=='external':s.streets[30]['endpoints']=(100,101)
            else:c[1].leaving=True
            links,retained=self.review.read(s,j,c);self.assertFalse(links);self.assertTrue(retained)

    def test_unsupported_feature_not_absorbed(self):
        s,j,c=sample();s.streets[30]['status']='UNSUPPORTED'
        links,retained=self.review.read(s,j,c);self.assertFalse(links);self.assertTrue(retained)

    def test_ambiguous_bundle_end_not_selected(self):
        s,j,c=sample();other=deepcopy(c[0]);other.item.id=8;c.append(other)
        links,retained=self.review.read(s,j,c);self.assertFalse(links);self.assertTrue(retained)

    def test_side_and_origin_not_hardcoded(self):
        s,j,c=sample();c[0].leaving=False
        s.streets[30]['items'][0]['points']=((-6.,100.,0.),(6.,100.,0.))
        links,_=self.review.read(s,j,c);self.assertEqual(links[0].bundle_side,1)

    def test_empty_degenerate_nonfinite_and_reversing_sources_retained(self):
        chains=([],[(-6.,0.,0.)],
                [(-6.,0.,0.),(-6.,0.,0.),(6.,0.,0.)],
                [(-6.,0.,0.),(float('nan'),0.,0.),(6.,0.,0.)],
                [(-6.,0.,0.),(0.,0.,float('inf')),(6.,0.,0.)],
                [(-6.,0.,0.),(2.,0.,0.),(-2.,0.,0.),(6.,0.,0.)])
        for points in chains:
            with self.subTest(points=points):
                s,j,c=sample();s.streets[30]['items'][0]['points']=points
                links,retained=self.review.read(s,j,c)
                self.assertFalse(links);self.assertEqual(len(retained),1)
        s,j,c=sample();s.streets[30]['items']=[]
        links,retained=self.review.read(s,j,c)
        self.assertFalse(links);self.assertEqual(len(retained),1)

    def test_disconnected_or_repeated_sections_retained(self):
        for duplicate in (False,True):
            s,j,c=sample();section=s.streets[30]['items'][0]
            s.streets[30]['items']=[dict(section,points=((-6.,0.,0.),(-1.,0.,0.))),
                dict(section,id=300 if duplicate else 301,
                     points=((-1. if duplicate else 1.,0.,0.),(6.,0.,0.)))]
            links,retained=self.review.read(s,j,c)
            self.assertFalse(links);self.assertEqual(len(retained),1)


class InteriorAdapterTests(unittest.TestCase):
    def adapter(self,snapshot,connectors,reader=None):
        self.calls=[]
        def external(intersection,**options):
            self.calls.append((intersection,options))
            return NS(source_id=intersection.id,connectors=options['source_connectors'])
        source_adapter=NS(_connectors=lambda j:tuple(connectors),from_blosm_intersection=external)
        contract=NS(StreetInput=lambda sid,sections:NS(source_id=sid,sections=sections))
        if reader is None:
            reader=lambda street:tuple(NS(source_id=s['id'],centerline=s['points'],width=s['width'],offset=s['offset'] or 0.)
                for s in snapshot.streets[street.id]['items'])
        return review.PreparedJunctionInteriorAdapter(contract,source_adapter,reader)

    def test_original_sources_and_external_selection_both_directions(self):
        for reverse in (False,True):
            s,j,c=sample(reverse);before=deepcopy(s);packet=self.adapter(s,c).read(s,j,profiles_by_connector={})
            self.assertEqual(s,before);self.assertEqual(packet.original_connector_ids,(10,12,11))
            self.assertEqual(packet.junction.connectors,(c[0],))
            self.assertEqual(packet.internal_streets[0].sections[0].centerline,s.streets[30]['items'][0]['points'])
            self.assertIs(self.calls[0][1]['source_connectors'][0],c[0])
            self.assertEqual(packet.links[0].endpoint_members,(1,0) if reverse else (0,1))

    def test_complete_multisection_source_and_options(self):
        s,j,c=sample();source=s.streets[30]['items'][0]
        s.streets[30]['items']=[dict(source,points=((-6.,0.,0.),(0.,0.,0.))),
                              dict(source,id=301,points=((0.,0.,0.),(6.,0.,0.)))]
        packet=self.adapter(s,c).read(s,j,profiles_by_connector={},right_hand_traffic=False)
        self.assertEqual(tuple(v.source_id for v in packet.internal_streets[0].sections),(300,301))
        self.assertFalse(self.calls[0][1]['right_hand_traffic'])

    def test_changed_source_reader_is_rejected(self):
        for field,value in (('source_id',999),('width',2.),('offset',1.),('centerline',((-6.,0.,0.),(7.,0.,0.)))):
            s,j,c=sample();source=s.streets[30]['items'][0]
            section=dict(source_id=300,centerline=source['points'],width=4.,offset=0.);section[field]=value
            adapter=self.adapter(s,c,lambda street:(NS(**section),))
            with self.assertRaises(ValueError):adapter.read(s,j,profiles_by_connector={})
            self.assertFalse(self.calls)

    def test_returning_source_stays_on_existing_path(self):
        s,j,c=sample();s.streets[30]['items'][0]['points']=((-6.,0.,0.),(-6.,-30.,0.),(6.,0.,0.))
        adapter=self.adapter(s,c)
        self.assertIsNone(adapter.read(s,j,profiles_by_connector={}));self.assertFalse(self.calls)


if __name__=='__main__':unittest.main()
