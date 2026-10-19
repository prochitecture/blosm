"""Account for accepted native surface owners, including partial components.

The presence of some faces is not evidence that every declared route succeeded.
These records map the solver's exact FACE owner keys back to prepared Streets,
physical/logical junctions and Bundles without inspecting geometric proximity.
"""
from dataclasses import dataclass
from numbers import Integral


def matching_connection_records(components, site, mouth_ids):
    """Require one checked seam status and its two mouths across components.

    Mesh and point-cloud storage share the same role/identity contract. Never
    accept each component independently: duplicate or conflicting records
    across the components must reject, just as duplicates in a legacy mesh do.
    Components without any relevant records need no record-specific columns.
    """
    schemas = {3: {'_cn_cache_site': 'INT', 'valid': 'BOOLEAN'},
               2: {'st_cj_input_index': 'INT', 'st_cj_seam': 'BOOLEAN', 'st_cj_valid': 'BOOLEAN'}}
    records = []
    for component in components:
        if component is None: continue
        attrs = component.attributes
        role = attrs.get('_cn_cache_part')
        if role is None: continue
        if role.domain != 'POINT' or role.data_type != 'INT': return False
        selected = [i for i, item in enumerate(role.data) if item.value in (2, 3)]
        if not selected: continue
        for i in selected:
            part = role.data[i].value
            schema = schemas[part]
            if any(name not in attrs or attrs[name].domain != 'POINT' or
                   attrs[name].data_type != kind or len(attrs[name].data) != len(role.data)
                   for name, kind in schema.items()):
                return False
            records.append(dict(_cn_cache_part=part,
                                **{name: attrs[name].data[i].value for name in schema}))
    status = [row for row in records if row['_cn_cache_part'] == 3]
    mouths = [row for row in records if row['_cn_cache_part'] == 2]
    # The provider writes the site on its status row only. Mouths carry global
    # input indices; their padded _cn_cache_site value has no meaning.
    return (len(status) == 1 and len(mouths) == 2 and len(mouth_ids) == 2 and len(set(mouth_ids)) == 2
            and bool(status[0]['valid'])
            and status[0]['_cn_cache_site'] == site
            and all(row['st_cj_seam'] and row['st_cj_valid'] for row in mouths)
            and sorted(row['st_cj_input_index'] for row in mouths) == sorted(mouth_ids))


@dataclass(frozen=True)
class OriginalSurfaceOwner:
    kind: int
    index: int
    streets: tuple = ()
    junctions: tuple = ()
    bundles: tuple = ()
    connection_mouths: tuple = ()
    connection_owners: tuple = ()

    @property
    def key(self):
        return self.kind, self.index


@dataclass(frozen=True)
class SeparatorRegionOwner:
    bundle_id: int
    region_id: int
    streets: tuple = ()

    @property
    def key(self):return self.bundle_id,self.region_id


class OriginalSurfaceCoverage:
    @staticmethod
    def free_bundle_owners(rows):
        """One native union owns both full members, including any grass islands.

        Rows retain (Bundle, Street, Section, chain, order, owner). Unlike a
        separator-only owner, a missing kind-8 owner means both pavements failed.
        """
        groups={};seen=set();chain_streets={};street_bundles={}
        for bid,sid,section,chain,order,owner in rows:
            if any(not isinstance(v,Integral) or isinstance(v,bool) or v<0
                   for v in (bid,sid,section,chain,order)) or type(owner) is not bool:
                raise ValueError('Invalid free Bundle ownership identity')
            if section in seen:raise ValueError('Repeated free Bundle Section ownership')
            seen.add(section)
            if chain in chain_streets and chain_streets[chain]!=(bid,sid):
                raise ValueError('Bundle chain has conflicting Street ownership')
            if sid in street_bundles and street_bundles[sid]!=bid:
                raise ValueError('Street has several whole-Bundle owners')
            chain_streets[chain]=(bid,sid);street_bundles[sid]=bid
            groups.setdefault(bid,{}).setdefault((sid,chain),[]).append((order,owner))
        result=[]
        for bid,chains in sorted(groups.items()):
            streets={sid for sid,_ in chains}
            if len(chains)!=2 or len(streets)!=2:
                raise ValueError('Free Bundle surface requires exactly two complete member chains')
            owner_count=0
            for chain in chains.values():
                chain.sort()
                if [row[0] for row in chain]!=list(range(len(chain))) or len({row[1] for row in chain})!=1:
                    raise ValueError('Free Bundle chain order or owner is inconsistent')
                owner_count+=int(chain[0][1])
            if owner_count!=1:raise ValueError('Free Bundle needs one union-generating member')
            result.append(OriginalSurfaceOwner(8,bid,tuple(sorted(streets)),bundles=(bid,)))
        return tuple(result)

    @staticmethod
    def shared_separator_owners(rows):
        """One region owns faces jointly for its complete original Bundles.

        Each row is (region namespace, patch index, Bundle ID, Street ID).
        Surface kind 7 keeps this relationship separate from kind 4, whose
        owner is one original Bundle. Neither original Bundle ID is replaced.
        """
        groups = {}
        for row in rows:
            if len(row) != 4 or any(isinstance(v, bool) or not isinstance(v, Integral)
                                    or not 0 <= v < 2**20 for v in row):
                raise ValueError('Shared separator needs bounded integer original identities')
            region, patch, bundle, street = map(int, row)
            groups.setdefault(region, []).append((patch, bundle, street))
        result, patches = [], set()
        for region, members in sorted(groups.items()):
            owners = {p for p, _, _ in members}
            bundles = {b for _, b, _ in members}
            streets = {s for _, _, s in members}
            if (not 3 <= len(members) <= 64 or len(owners) != 1 or
                    owners.intersection(patches) or len(streets) != len(members) or
                    region <= max(bundles) or len(bundles) < 2 or
                    any(sum(b == bid for _, b, _ in members) != 2 for bid in bundles)):
                raise ValueError('Shared region must retain complete, unique two-member Bundle ownership')
            patches.update(owners)
            result.append(OriginalSurfaceOwner(7, region, tuple(sorted(streets)), bundles=tuple(sorted(bundles))))
        return tuple(result)

    def __init__(self, owners, regions=()):
        self.owners = {owner.key: owner for owner in owners}
        if len(self.owners) != len(owners):
            raise ValueError('Repeated original surface owner declaration')
        self.regions={region.key:region for region in regions}
        if len(self.regions)!=len(regions):raise ValueError('Repeated separator region declaration')
        if any((4,r.bundle_id) not in self.owners for r in regions):
            raise ValueError('Separator region has no original Bundle surface owner')

    def evaluate(self, accepted, *, matching_connections=(), accepted_regions=()):
        accepted = {(int(kind), int(index)) for kind, index in accepted}
        unknown = accepted.difference(self.owners)
        if unknown:
            raise ValueError('Undeclared native surface owners: '+str(sorted(unknown)))
        connections = set()
        for key in matching_connections:
            key = tuple(key)
            owner = self.owners.get(key)
            if owner is None or len(owner.connection_mouths) != 2:
                raise ValueError('Undeclared native matching connection: '+str(key))
            # A matching seam has no area. Its native provider must be valid,
            # and both declared incident surfaces must survive assembly. An
            # unequal-profile transition still needs its own pavement faces.
            if len(owner.connection_owners) == 2 and set(owner.connection_owners) <= accepted:
                connections.add(key)
        accepted.update(connections)
        missing = set(self.owners).difference(accepted)
        retained_streets = {sid for key in accepted for sid in self.owners[key].streets}
        affected = {sid for key in missing for sid in self.owners[key].streets}
        regions={tuple(map(int,key)) for key in accepted_regions}
        if regions.difference(self.regions):raise ValueError('Undeclared native separator regions')
        if any((4,bundle) not in accepted for bundle,_ in regions):
            raise ValueError('Region evidence lacks its Bundle surface')
        missing_regions=set(self.regions).difference(regions)
        affected.update(sid for key in missing_regions for sid in self.regions[key].streets)
        result=dict(complete=not missing and not missing_regions,
            accepted_owners=[list(key) for key in sorted(accepted)],
            missing_owners=[list(key) for key in sorted(missing)],
            affected_streets=sorted(affected),
            failed_streets=sorted(affected.difference(retained_streets)),
            missing_junctions=sorted({jid for key in missing for jid in self.owners[key].junctions}),
            missing_bundles=sorted({bid for key in missing for bid in self.owners[key].bundles}|
                                   {bid for bid,_ in missing_regions}))
        if self.regions:
            result.update(accepted_regions=[list(key) for key in sorted(regions)],
                          missing_regions=[list(key) for key in sorted(missing_regions)])
        return result

    def matching_connections(self, obj):
        """Read checked face-free joins from evaluated native patch providers."""
        from .evidence import evaluated_record_components
        from ..gn.editable_input.curved_patch_provider_bindings import input_values
        providers = input_values(obj)['Patch Results'].objects
        wanted = {owner.index: owner for owner in self.owners.values()
                  if owner.kind == 0 and len(owner.connection_mouths) == 2}
        completed = []
        for provider in providers:
            site = provider.get('st_cached_patch_owner')
            if site not in wanted:
                continue
            with evaluated_record_components(provider) as (mesh, points):
                if not mesh.polygons and matching_connection_records(
                        (mesh, points), site, wanted[site].connection_mouths):
                    completed.append((0, site))
        return tuple(completed)

    @classmethod
    def from_object(cls, obj):
        from ..gn.editable_input.curved_patch_provider_bindings import input_values, values, owner_sites
        tables = input_values(obj)
        def data(name, field):
            return values(tables[name].data, field)
        junctions = data('Junctions', 'st_cj_junction_id')
        starts = data('Junctions', 'st_cj_range_start')
        counts = data('Junctions', 'st_cj_range_count')
        streets = data('Approaches', 'st_cj_street_id')
        logical = dict(zip(data('Logical Clusters', 'st_lg_site'), data('Logical Clusters', 'st_lg_id')))
        import json
        logical_ids=json.loads(obj.get('blosm_logical_owner_ids','{}'))
        incident = {}
        for kind, name, prefix, fields in ((1, 'Corridors', 'st_cs_', ('start_index', 'end_index')),
                                           (3, 'External Approaches', 'st_ca_', ('start_index',))):
            for field in fields:
                for index, mouth in enumerate(data(name, prefix+field)):
                    if mouth in incident:
                        raise ValueError('Repeated original mouth surface ownership')
                    incident[mouth] = (kind, index)
        kinds = tables['Approaches'].data.attributes.get('st_pp_kind')
        owners = []
        for owner, sites in owner_sites(tables).items():
            ids = {sid for site in sites for sid in streets[starts[site]:starts[site]+counts[site]]}
            jids = {junctions[site] for site in sites}
            if owner in logical:
                jids.update(logical_ids.get(str(logical[owner]),(logical[owner],)))
            mouths = tuple(range(starts[owner], starts[owner]+counts[owner]))
            matching = (len(sites) == 1 and len(mouths) == 2 and kinds is not None
                        and all(kinds.data[i].value == 3 for i in mouths))
            owners.append(OriginalSurfaceOwner(0, owner, tuple(sorted(ids)), tuple(sorted(jids)),
                connection_mouths=mouths if matching else (),
                connection_owners=tuple(incident[i] for i in mouths if i in incident) if matching else ()))
        for kind, name, prefix in ((1, 'Corridors', 'st_cs_'), (3, 'External Approaches', 'st_ca_')):
            owners.extend(OriginalSurfaceOwner(kind, index, (sid,))
                          for index, sid in enumerate(data(name, prefix+'street_id')))
        if tables.get('Separators'):
            for bid, a, b in zip(data('Separators', 'st_sb_bundle_id'),
                                 data('Separators', 'st_sb_member_a'), data('Separators', 'st_sb_member_b')):
                owners.append(OriginalSurfaceOwner(4, bid, (a, b), bundles=(bid,)))
        if tables.get('Standalone Streets'):
            owners.extend(OriginalSurfaceOwner(5,index,(sid,))
                          for index,sid in enumerate(data('Standalone Streets','st_fs_street_id')))
        if tables.get('Boundary Runs'):
            owners.extend(OriginalSurfaceOwner(6,index,(sid,))
                          for index,sid in enumerate(data('Boundary Runs','st_br_street_id')))
        if tables.get('Internal Cycles'):
            owners.extend(cls.shared_separator_owners(zip(*(data('Internal Cycles', 'st_iv_'+key)
                for key in ('region_id', 'patch', 'bundle_id', 'street_id')))))
        if tables.get('Bundle Paths'):
            rows=list(zip(*(data('Bundle Paths','st_path_'+key) for key in
                ('bundle_id','street_id','id','chain_index','chain_order','chain_owner'))))
            declarations=json.loads(obj.get('blosm_free_bundle_preparation','[]'))
            expected={(row['bundle_id'],member['street_id'],section['start']['section_id'])
                for row in declarations for member in row['members'] for span in member['spans'] for section in span['sections']}
            if {tuple(row[:3]) for row in rows}!=expected:
                raise ValueError('Free Bundle paths lost original prepared Section ownership')
            owners.extend(cls.free_bundle_owners(rows))
        regions=[]
        if tables.get('Continuation Regions') and len(tables['Continuation Regions'].data.vertices):
            regions.extend(SeparatorRegionOwner(owner.index,owner.index,owner.streets) for owner in owners if owner.kind==4)
            regions.extend(SeparatorRegionOwner(bid,rid,(a,b)) for bid,rid,a,b in zip(
                data('Continuation Regions','st_sb_bundle_id'),data('Continuation Regions','st_sb_region_id'),
                data('Continuation Regions','st_sb_member_a'),data('Continuation Regions','st_sb_member_b')))
        return cls(owners,regions)
