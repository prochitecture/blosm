"""Source-incidence partition of a Bundle at its declared interior crossings.

This descriptor changes no source geometry or membership. Original Section
incidences define the run boundaries; physical mouths are resolved separately.
"""
from dataclasses import dataclass
from . import model
from .bundle_axis import AxialChartUnavailable


class RunPartitionUnavailable(AxialChartUnavailable):
    pass


@dataclass(frozen=True)
class MemberRun:
    visits: tuple
    ends: tuple
    boundaries: tuple


@dataclass(frozen=True)
class BundleRun:
    members: tuple
    boundaries: tuple

    @property
    def source_sections(self):
        return tuple(visit.section_id for member in self.members for visit in member.visits)


class VisibleBundleRuns:
    @staticmethod
    def build(network, bundle, crossing_inventory):
        if len(bundle.members) != 2 or not all(member.complete for member in bundle.members):
            raise RunPartitionUnavailable('A paired run partition needs two complete ordered members')
        sites = {end:site for site in (*network.intersections,*network.connections) for end in site.ends}
        by_site = {}
        for crossing in crossing_inventory:
            for role,site in enumerate(crossing.sites):
                token=(role,site)
                if token in by_site:
                    raise RunPartitionUnavailable('A member site needs a general interior-region partition')
                by_site[token]=('crossing',crossing.street_id)
        sides=[]
        for role,member in enumerate(bundle.members):
            runs=[];visits=[];first_token=('bundle_end',0)
            for index,visit in enumerate(member.path):
                entry=model.End(visit.section_id,0 if visit.forward else 1)
                exit=model.End(visit.section_id,1 if visit.forward else 0)
                if index:
                    previous=member.path[index-1]
                    previous_end=model.End(previous.section_id,1 if previous.forward else 0)
                    site=sites.get(previous_end)
                    if site is None or entry not in site.ends:
                        raise RunPartitionUnavailable('Member runs require explicit original Section incidences')
                visits.append(visit)
                site=sites.get(exit)
                last_token=by_site.get((role,site.key)) if site is not None else None
                if index==len(member.path)-1:
                    if last_token is not None:
                        raise RunPartitionUnavailable('An interior crossing at a member terminal needs its continuation region')
                    last_token=('bundle_end',1)
                if last_token is None:
                    continue
                first=visits[0]
                start=model.End(first.section_id,0 if first.forward else 1)
                runs.append(MemberRun(tuple(visits),(start,exit),(first_token,last_token)))
                visits=[];first_token=last_token
            if visits:
                raise RunPartitionUnavailable('Incomplete member-run partition')
            sides.append(tuple(runs))
        if len(sides[0])!=len(sides[1]) or any(a.boundaries!=b.boundaries for a,b in zip(*sides)):
            raise RunPartitionUnavailable('Interior crossings need a region graph beyond a common ordered run chain')
        result=tuple(BundleRun((a,b),a.boundaries) for a,b in zip(*sides))
        used=[sid for run in result for sid in run.source_sections]
        original=[visit.section_id for member in bundle.members for visit in member.path]
        if len(used)!=len(set(used)) or set(used)!=set(original):
            raise RunPartitionUnavailable('Run partition must retain every member Section exactly once')
        return result
