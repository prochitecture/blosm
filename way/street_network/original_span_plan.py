"""Input-independent import selection for the native original-span provider.

Selection follows original incidence identities. It never invents free ends,
fits surfaces, or chooses a different radius to make a connection succeed.
"""
from dataclasses import dataclass, replace
import math
from .physical_ports import PreparedPhysicalPorts
from .port_connections import PreparedPortConnections
from .port_compounds import PreparedPortCompound
from .port_network import PreparedPortNetwork
from .port_clusters import PreparedLogicalPortClusters


@dataclass(frozen=True)
class OriginalSpanComponent:
    ports: object
    network: object
    frontiers: tuple
    compound_decisions: tuple = ()
    free_bundles: tuple = ()


@dataclass(frozen=True)
class OriginalSpanImportPlan:
    components: tuple
    issues: tuple


class OriginalSpanPlanner:
    def __init__(self, categories, *, category_properties=None, radius_for_end=None, defer_physical_categories=True, section_scopes=False,
                 internal_logical_ends=False, shared_logical_owners=False, seam_connections=False,
                 category_connections=False, boundary_runs=False, dimension_runs=False, free_bundles=False):
        self.reader = PreparedPhysicalPorts(categories, defer_physical_categories=defer_physical_categories,
                                           section_scopes=section_scopes, seam_connections=seam_connections,
                                           category_connections=category_connections)
        self.category_properties=category_properties
        # Authored inputs may have a different radius at each original port.
        # Reuse the ownership policy without replacing those values by a
        # category default. The callable receives (inventory, SpanEnd).
        self.radius_for_end=radius_for_end
        self.internal_logical_ends=internal_logical_ends
        self.shared_logical_owners=shared_logical_owners
        self.boundary_runs=boundary_runs
        self.dimension_runs=dimension_runs
        self.free_bundles=free_bundles
        if boundary_runs and not section_scopes:
            raise ValueError('Boundary runs require explicit Section category scopes')
        if shared_logical_owners and not internal_logical_ends:
            raise ValueError('Shared logical owners require retained internal ends')

    def short_compounds(self, snapshot, ports, logical, explicit):
        """Bounded import-time ownership heuristic; GN still proves geometry.

        A span shorter than its width plus both prescribed terminal radii is
        a candidate for the existing two-site construction. This is not a
        fitted mouth or a geometric feasibility certificate. Ownership remains
        stable during coordinate edits; no source points/radii are modified.
        """
        if self.category_properties is None and self.radius_for_end is None:return explicit,()
        ordinary=PreparedPortNetwork(PreparedPortCompound()).prepare(ports,
            compound_spans=explicit,logical_clusters=logical)
        occupied={s for c in ordinary.compounds for s in c.site_ids}
        occupied.update(s for c in logical for s in c.physical_junction_ids)
        candidates=[]
        for route in ordinary.routes:
            if route.kind!='CORRIDOR':continue
            span=route.span
            sections=[s for s in snapshot.streets[span.street_id]['items']
                      if s['kind']=='Section' and s['id'] in span.section_ids]
            length=sum(math.dist(a,b) for s in sections for a,b in zip(s['points'],s['points'][1:]))
            width=max(float(s['width']) for s in sections)
            try:
                radii=([float(self.radius_for_end(snapshot,end)) for end in (span.start,span.end)]
                       if self.radius_for_end is not None else
                       [float(self.category_properties[s['category']]['radius']) for s in (sections[0],sections[-1])])
            except (KeyError,TypeError,ValueError):continue
            limit=width+sum(radii)
            if all(math.isfinite(v) and v>0 for v in (length,width,*radii)) and length<limit:
                candidates.append((length/limit,span.identity,length,limit))
        selected=list(explicit);decisions=[]
        for _ratio,identity,length,limit in sorted(candidates):
            try:compound=PreparedPortCompound().prepare(ports,identity)
            except ValueError as exc:
                decisions.append(dict(span=identity,status='REJECTED',reason=str(exc)));continue
            if occupied.intersection(compound.site_ids):
                decisions.append(dict(span=identity,status='DEFERRED',reason='Sites already have a surface owner'));continue
            occupied.update(compound.site_ids);selected.append(identity)
            decisions.append(dict(span=identity,status='SELECTED',source_length_m=length,
                candidate_limit_m=limit,rule='source length < own width + terminal radii'))
        return tuple(selected),tuple(decisions)

    @staticmethod
    def site_kind(site):
        if site.ports and all(p.source_kind=='SEAM' for p in site.ports):return 'SEAM'
        return 'PHYSICAL' if all(p.source_kind=='PHYSICAL' for p in site.ports) else 'INTERIOR'

    @staticmethod
    def issues_for(snapshot, component, issues):
        sites={s.junction_id for s in component.ports.sites}
        logical={jid for row in snapshot.junction_origins.get('components',())
            if sites.intersection(row['physical_junction_ids']) for jid in row['logical_junction_ids']}
        return tuple(row for row in issues
            if (row['kind']=='LOGICAL' and row.get('id') in logical) or
               (row['kind']=='LOGICAL_COMPONENT' and sites.intersection(row.get('ids',()))))

    def prepare(self, snapshot, *, physical_ids=None, interior_ids=(), seam_ids=(), expand=True,
                compound_spans=(), standalone=False):
        """None selects all original sites; an explicit tuple selects seeds.

        Compound declarations are explicit topology policy while automatic
        short-span ownership is being completed. Invalid declarations reject;
        unsupported automatic seeds remain named issues in the import report.
        """
        automatic = physical_ids is None
        if automatic:
            category_seams = self.reader.category_seam_ids(snapshot)
            if self.reader.seam_connections:
                seam_ids = tuple(sorted(s['junction_id'] for s in snapshot.junction_origins.get('seams', ())
                                        if s['junction_id'] not in category_seams))
            physical_ids = tuple(sorted({r['junction_id'] for r in snapshot.junction_origins.get('physical', ())}
                                         -set(seam_ids)-category_seams))
            interior_ids = tuple(sorted({i['id'] for s in snapshot.streets.values() for i in s['items']
                if i['kind']=='Intersection' and i.get('surface_treatment')=='NATIVE_JUNCTION'}))
        groups=[];issues=[];seen=set()
        seeds = [('PHYSICAL',i) for i in physical_ids]+[('INTERIOR',i) for i in interior_ids]+[('SEAM',i) for i in seam_ids]
        # Explicit selections are validated together; automatic selection
        # keeps an unsupported complete ring from hiding independent components.
        selections = [[seed] for seed in seeds] if automatic else [seeds]
        for selection in selections:
            if not selection or all(seed in seen for seed in selection):
                continue
            physical=tuple(i for kind,i in selection if kind=='PHYSICAL')
            interior=tuple(i for kind,i in selection if kind=='INTERIOR')
            seams=tuple(i for kind,i in selection if kind=='SEAM')
            try:
                if expand:
                    found=PreparedPortConnections(self.reader).expand(snapshot,physical_ids=physical,interior_ids=interior,seam_ids=seams)
                    ports=found.plan;frontiers=found.frontiers
                else:
                    ports=self.reader.prepare(snapshot,physical_ids=physical,interior_ids=interior,seam_ids=seams);frontiers=()
            except ValueError as exc:
                if not automatic:
                    raise
                issues.append(dict(kind=selection[0][0],id=selection[0][1],reason=str(exc)))
                seen.update(selection);continue
            seen.update((self.site_kind(s),s.junction_id) for s in ports.sites)
            groups.append((ports,frontiers))

        # A complete authored Street must not be copied into two independently
        # editable components merely because a remote junction is unresolved.
        merged=[]
        for ports,frontiers in groups:
            members=set(ports.source_streets);matches=[]
            for at,(existing,_) in enumerate(merged):
                if members.intersection(existing.source_streets):matches.append(at)
            if matches:
                sites={s.junction_id:s for s in ports.sites};all_frontiers=list(frontiers)
                for at in reversed(matches):
                    old,old_frontiers=merged.pop(at)
                    sites.update((s.junction_id,s) for s in old.sites);all_frontiers.extend(old_frontiers)
                ports=self.reader.prepare(snapshot,
                    physical_ids=tuple(sorted(s.junction_id for s in sites.values() if self.site_kind(s)=='PHYSICAL')),
                    interior_ids=tuple(sorted(s.junction_id for s in sites.values() if self.site_kind(s)=='INTERIOR')),
                    seam_ids=tuple(sorted(s.junction_id for s in sites.values() if self.site_kind(s)=='SEAM')))
                frontiers=tuple(dict.fromkeys(all_frontiers))
            merged.append((ports,frontiers))

        components=[];used_compounds=set()
        for ports,frontiers in merged:
            selected={s.junction_id for s in ports.sites};logical=[]
            for row in snapshot.junction_origins.get('components', ()):
                ids=set(row['physical_junction_ids'])
                if len(ids)<2 or not ids.intersection(selected):continue
                if not ids<=selected:
                    issues.append(dict(kind='LOGICAL_COMPONENT',ids=sorted(ids),reason='Incomplete physical component remains unresolved'))
                    continue
                owners=row['logical_junction_ids']
                # One physical component has one surface, even when several
                # original logical junctions share it. Keep all IDs in the
                # declaration; do not create duplicate patch providers.
                for jid in ((min(owners),) if self.shared_logical_owners else owners):
                    try:logical.append(PreparedLogicalPortClusters().prepare(snapshot,ports,jid,
                        allow_internal_connectors=self.internal_logical_ends,
                        allow_shared_owners=self.shared_logical_owners))
                    except ValueError as exc:
                        issues.append(dict(kind='LOGICAL',id=jid,reason=str(exc)))
            known={s.identity for s in ports.spans}
            compounds=tuple(tuple(s) for s in compound_spans if tuple(s) in known)
            used_compounds.update(compounds)
            compounds,decisions=self.short_compounds(snapshot,ports,logical,compounds)
            network=PreparedPortNetwork(PreparedPortCompound()).prepare(ports,
                compound_spans=compounds,logical_clusters=logical)
            components.append(OriginalSpanComponent(ports,network,frontiers,decisions))
        if used_compounds!={tuple(s) for s in compound_spans}:
            raise ValueError('A declared compound span is absent from the imported components')
        if automatic and standalone:
            # Sources are selected by declared ownership, never their distance
            # from the generated network. Recheck physical/interior incidences
            # even when obsolete logical endpoint fields look unattached.
            used={sid for component in components for sid in component.ports.source_streets}
            if self.free_bundles:
                from .original_span_free_bundles import PreparedFreeBundles
                for bid, bundle in sorted(snapshot.bundles.items()):
                    members=set(bundle.get('head', ()))
                    if not members or members <= used:continue
                    try:
                        ports,network,surface=PreparedFreeBundles(self.reader).prepare(snapshot,bid,occupied=used)
                        components.append(OriginalSpanComponent(ports,network,(),free_bundles=(surface,)))
                        used.update(ports.source_streets)
                    except ValueError as exc:
                        issues.append(dict(kind='BUNDLE',id=bid,reason=str(exc)))
            for sid,source in sorted(snapshot.streets.items()):
                if sid in used or source.get('status')!='SUPPORTED':continue
                if tuple(source.get('endpoints',()))!=(None,None):continue
                sections=[item for item in source['items'] if item['kind']=='Section']
                if not sections or any(s['category'] not in self.reader.categories for s in sections):continue
                try:
                    if source.get('bundle_id') is not None:
                        raise ValueError('Unattached Bundle members need the Bundle surface branch')
                    if not self.dimension_runs and len(sections)!=len(source['items']):
                        raise ValueError('Unattached Street items need their explicit surface treatment')
                    if not self.dimension_runs and any((s['width'],s.get('offset') or 0.)!=(sections[0]['width'],sections[0].get('offset') or 0.) for s in sections):
                        raise ValueError('Unattached Street dimensional transitions need the tapered branch')
                    ports=self.reader.prepare(snapshot,street_ids=(sid,))
                    if self.dimension_runs:
                        from .physical_dimensions import PreparedPhysicalDimensions
                        for item in source['items']:
                            if item['kind']!='Section' and (item['kind'] not in ('SymLane','SideLane','AuthoredTransition')
                                    or item.get('surface_treatment') not in ('WIDTH_TRANSITION','COMPATIBLE_SEAM')
                                    or item.get('points') or item.get('transition_reason')):
                                raise ValueError('Unattached Street item has no dimensional run treatment')
                        PreparedPhysicalDimensions().read_spans(snapshot,ports.spans,source_plan=ports)
                    network=PreparedPortNetwork(PreparedPortCompound()).prepare(ports)
                    allowed=('STANDALONE','BOUNDARY_ONLY') if self.boundary_runs else ('STANDALONE',)
                    if not network.routes or any(route.kind not in allowed for route in network.routes):
                        raise ValueError('Street has a recorded junction incidence; it is not unattached')
                    components.append(OriginalSpanComponent(ports,network,()))
                except ValueError as exc:
                    issues.append(dict(kind='STREET',id=sid,reason=str(exc)))
        if self.boundary_runs:
            from .original_boundary_runs import OriginalBoundaryRuns
            if automatic:
                used={sid for component in components for sid in component.ports.source_streets}
                for sid,source in sorted(snapshot.streets.items()):
                    if sid in used or not any(s['kind']=='Section' and s['category'] in self.reader.categories
                                              for s in source['items']):continue
                    try:
                        ports=self.reader.prepare(snapshot,street_ids=(sid,))
                        network=PreparedPortNetwork(PreparedPortCompound()).prepare(ports)
                        # Keep unselected junctions pending; never replace them
                        # with an independent run merely to obtain a surface.
                        if not any(r.kind=='BOUNDARY_ONLY' for r in network.routes):continue
                        components.append(OriginalSpanComponent(ports,network,()))
                    except ValueError as exc:
                        issues.append(dict(kind='STREET',id=sid,reason=str(exc)))
            admitted=[]
            for component in components:
                runs,pending=OriginalBoundaryRuns().prepare(snapshot,component.ports,component.network,
                                                          dimensions=self.dimension_runs)
                issues.extend(pending)
                admitted.append(replace(component,network=replace(component.network,boundary_runs=runs)))
            components=admitted
        return OriginalSpanImportPlan(tuple(components),tuple(issues))
