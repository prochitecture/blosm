"""Local merged mouths retain complete multi-Section source-run identity."""
from dataclasses import replace
import unittest
import numpy as np
import _python_street_package
from test_python_street_merged_entrances import network
from tested_python_streets import model
from tested_python_streets.backend import PythonMeshBackend
from tested_python_streets.bundle_merged_entrances import MergedEntrancePreparation, PairedAngularEntrance, preparation_key
from tested_python_streets.conditions import BoundaryConditions
from tested_python_streets.updates import BuildContext, DependencyIndex


def split_network(reverse=False, state='overlap'):
    source=network(3,state)
    sections=[]
    for s in source.sections:
        if s.id==17:
            sections.extend((replace(s,id=18,points=s.points[:3]),replace(s,points=s.points[2:])))
        elif s.id==129:
            sections.extend((replace(s,points=s.points[:-2]),replace(s,id=130,points=s.points[-3:])))
        else:sections.append(s)
    def endpoint(e):
        if e==model.End(17,0):return model.End(18,0)
        if e==model.End(129,1):return model.End(130,1)
        return e
    sites=tuple(replace(j,ends=tuple(map(endpoint,j.ends))) for j in source.intersections)
    connections=(model.Connection(900,(model.End(18,1),model.End(17,0))),
                 model.Connection(901,(model.End(129,1),model.End(130,0))))
    bundle=replace(source.bundles[0],members=(model.BundleMember(tuple(model.Traversal(sid) for sid in (18,17,117))),
        model.BundleMember(tuple(model.Traversal(sid) for sid in (29,129,130)))))
    if reverse:
        sections=[replace(s,points=s.points[::-1],offset=-s.offset) for s in sections]
        sites=tuple(replace(j,ends=tuple(model.End(e.section_id,1-e.end) for e in j.ends)) for j in sites)
        connections=tuple(replace(j,ends=tuple(model.End(e.section_id,1-e.end) for e in j.ends)) for j in connections)
        bundle=replace(bundle,members=tuple(model.BundleMember(tuple(model.Traversal(v.section_id,False) for v in m.path))
                                            for m in bundle.members))
    return replace(source,sections=tuple(sections),intersections=sites,connections=connections,bundles=(bundle,))


def prepare(source):
    backend=PythonMeshBackend(source.settings)
    graph=DependencyIndex(tuple(task for provider in backend.providers for task in provider.tasks(source)))
    target=preparation_key(source.bundles[0]);wanted=set();pending=[target]
    while pending:
        key=pending.pop()
        if key in wanted or key not in graph.tasks:continue
        wanted.add(key);pending.extend(graph.tasks[key].dependencies)
    context=BuildContext(source,{},backend.patterns)
    for key in graph.order:
        if key in wanted:context.products[key]=graph.tasks[key].build(context)
    return context,context[target],graph.tasks[target]


class MergedSourceRunTests(unittest.TestCase):
    def test_optional_contact_knots_refine_original_curbs_without_replacing_their_error_budget(self):
        source=split_network();context,spec,_=prepare(source)
        original=next(a for a in spec.layout.approaches if isinstance(a,PairedAngularEntrance))
        supplied={p.section.id:np.array((0.,.337*p.length,p.length)) for p in original.paths}
        result=PairedAngularEntrance(context,original.run,original.end_index,spec.frame,
            source_tolerance=.00025,source_parameters=supplied)
        for path,visit,stations in zip(result.paths,result.visits,result.source_curve_stations):
            ordered=stations if visit.forward else stations[::-1]
            self.assertGreater(len(ordered),3)
            self.assertEqual(ordered[0],0.);self.assertEqual(ordered[-1],path.length)
            self.assertTrue(np.any(ordered==supplied[path.section.id][1]))
        for invalid in (0.,float('nan')):
            with self.assertRaises(ValueError):PairedAngularEntrance(context,original.run,original.end_index,spec.frame,
                source_tolerance=invalid)

    def test_actual_adjacent_sections_define_aliases_without_losing_any_complete_run_source(self):
        for reverse in (False,True):
            with self.subTest(reverse=reverse):
                source=split_network(reverse);before=source.to_dict()
                candidates=MergedEntrancePreparation.candidates(source)
                self.assertEqual(len(candidates),1)
                context,spec,task=prepare(source)
                self.assertIsNotNone(spec)
                paired=[a for a in spec.layout.approaches if isinstance(a,PairedAngularEntrance)]
                self.assertEqual(len(paired),2)
                member_sections={v.section_id for member in source.bundles[0].members for v in member.path}
                self.assertEqual({sid for a in paired for sid in a.source_run_sections},member_sections)
                self.assertEqual({path.section.id for a in paired for path in a.paths},{17,29,117,129})
                self.assertEqual({v.section_id for a in paired for m in a.run.members for v in m.visits},
                                 member_sections)
                for sid in source.bundles[0].sections:
                    for kind in ('path','profile'):
                        self.assertIn(model.Key(kind,sid),task.dependencies)
                self.assertTrue(all(d.kind not in ('bundle_plan','junction_surface','road','separator') for d in task.dependencies))
                conditions=BoundaryConditions().build(spec.group,spec.group.loops[0],
                    {f.region.owner:f.layout for f in spec.group.footprints})
                ports=[p for p in conditions.ports if len(p.aliases)==2]
                self.assertEqual(len(ports),2)
                for port in ports:
                    a=spec.layout.approaches[port.approach_index]
                    station=spec.layout.outward_stations[port.approach_index]
                    self.assertLess(station,a.length)
                    self.assertEqual({alias.endpoint for alias in port.aliases},set(a.original_ends))
                    for alias in port.aliases:
                        path=context[model.Key('path',alias.endpoint.section_id)]
                        self.assertAlmostEqual(float(spec.frame.local(path.frame(alias.station)[0])[0]),
                                               float(a.chart_station(station)),places=8)
                        profile=context[model.Key('profile',alias.endpoint.section_id)].oriented(alias.endpoint.end)
                        u=np.asarray(port.outward_u)[alias.contains(port.outward_u)]
                        np.testing.assert_allclose(profile.value(2.*alias.outward_parameter(u)-1.),
                                                   a.profile.value(2.*u-1.),atol=1e-12,rtol=0.)
                self.assertEqual(source.to_dict(),before)

    def test_separated_multi_section_source_pair_retains_its_existing_provider(self):
        source=split_network(state='multiple');before=source.to_dict()
        self.assertEqual(len(MergedEntrancePreparation.candidates(source)),1)
        _,spec,_=prepare(source)
        self.assertIsNone(spec)
        self.assertEqual(source.to_dict(),before)


if __name__=='__main__':unittest.main()
