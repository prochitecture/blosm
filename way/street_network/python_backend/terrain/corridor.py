"""Translate generated products into the complete corridor placement contract.

An extension may append CorridorComponents on an existing reference: raised
side strips, curbs and other surfaces then constrain the same fit as pavement.
No sidewalk meshes or material policy are hard-coded into the fitter.
"""
import numpy as np
from .alignment import CorridorComponent, ElevationReference
from ..model import Key
from ..providers import Part, RoadSampling


def boundary_vertices(quads):
    edges = np.sort(np.concatenate([quads[:,[i,(i+1)%4]] for i in range(4)]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    return np.unique(unique[counts == 1])


def prepare(network, products, settings):
    components, controls, membership = [], {}, {}
    for bundle in network.bundles:
        for identity in bundle.sections:
            membership[identity] = bundle.id
    section_index = network.section_index
    protected = set()
    for key, part in products.items():
        if not isinstance(part, Part) or part.mesh is None:
            continue
        owner = part.owner
        extension = part.metadata.get('terrain_corridor')
        if extension:
            sid = int(extension['section_id'])
            if sid not in section_index:
                raise ValueError('Unknown Section in component terrain contract')
            ref = 'bundle:'+str(membership[sid]) if sid in membership else Key('road',sid).text()
            stations = np.zeros(len(part.mesh.xyz)) if sid in membership else np.asarray(extension['stations'])
            components.append(CorridorComponent(key.text(),ref,part.mesh.xyz,part.mesh.quads,stations,
                float(extension.get('datum',section_index[sid].datum)),float(extension.get('clearance',0.))))
            continue
        bundle = owner.id if owner.kind == 'bundle' else membership.get(owner.id) if owner.kind == 'section' else None
        datum = float(part.metadata.get('terrain_datum', section_index[owner.id].datum
            if owner.kind == 'section' else min((s.datum for s in network.sections), default=0.)))
        if bundle is not None:
            ref = 'bundle:'+str(bundle)
            stations = np.zeros(len(part.mesh.xyz))
            protected.add(bundle)
        elif key.kind == 'road' and part.metadata.get('transverse_intervals'):
            ref = key.text()
            samples = part.metadata.get('terrain_stations')
            if samples is None:
                samples = RoadSampling.from_metadata(part.metadata).samples(products[Key('path', key.id)], network.settings)
            samples = np.asarray(samples)
            stations = np.repeat(samples, part.metadata['transverse_intervals']+1)
            if len(stations) != len(part.mesh.xyz):
                raise ValueError('Saved road station layout is obsolete; regenerate the street mesh')
            count = max(1, int(np.ceil((samples[-1]-samples[0])/settings.control_spacing)))
            # Controls finer than the published rows would have no incident
            # face constraints and could dip between rows. In that case use
            # the actual rows, also bounding allocation by existing geometry.
            controls[ref] = set(samples if (samples[-1]-samples[0])/count < np.max(np.diff(samples))
                                else np.linspace(samples[0], samples[-1], count+1))
        else:
            ref = 'rigid:'+owner.text()
            stations = np.zeros(len(part.mesh.xyz))
        if ref not in controls:
            controls[ref] = {0.}
        components.append(CorridorComponent(key.text(), ref, part.mesh.xyz, part.mesh.quads, stations, datum))
    # Connect real mesh seams, including composite junction ownership. The
    # numerical backend has already validated that these ground surfaces do
    # not overlap. Grade-separated ways are excluded during source preparation.
    buckets, contacts, serial = {}, {}, 0
    tolerance = 2e-5
    for component in components:
        if component.reference not in controls:
            raise ValueError('Component needs a surviving road reference; split it at junction-owned spans: '+component.identity)
        for vertex in boundary_vertices(component.quads):
            xy = component.xyz[vertex,:2]
            cell = tuple(np.floor(xy/tolerance).astype(np.int64))
            found = None
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    for other, index, point, token in buckets.get((cell[0]+dx,cell[1]+dy), ()):
                        if other is component or np.linalg.norm(point-xy) > tolerance:
                            continue
                        if abs((other.xyz[index,2]-other.datum)-(component.xyz[vertex,2]-component.datum)) > 1e-4:
                            raise ValueError('Terrain placement requires compatible prepared relative seam heights: '
                                             +other.identity+' / '+component.identity)
                        found = token
                        contacts.setdefault(other.identity,set()).add((float(other.stations[index]),token))
                        controls[other.reference].add(float(other.stations[index]))
                        break
                    if found is not None: break
                if found is not None: break
            if found is None:
                found = serial
                serial += 1
            else:
                contacts.setdefault(component.identity,set()).add((float(component.stations[vertex]),found))
                controls[component.reference].add(float(component.stations[vertex]))
            buckets.setdefault(cell,[]).append((component,vertex,xy,found))
    for component in components:
        component.contacts = tuple(sorted(contacts.get(component.identity,())))
    references = [ElevationReference(key, np.array(sorted(values))) for key,values in controls.items()]
    return references, components, dict(protected_bundles=sorted(protected),
        bundle_policy='Complete Bundle regions share one rigid elevation in this initial version')
