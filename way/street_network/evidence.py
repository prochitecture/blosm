"""Read native diagnostic/ownership attributes for per-entity import coverage."""
from contextlib import contextmanager
import bpy
import numpy as np


@contextmanager
def evaluated_mesh(obj, *, depsgraph=None):
    """Borrow evaluated Mesh data during synchronous, read-only inspection.

    Keep the evaluated object and dependency graph alive for the entire read.
    Callers must not edit or reevaluate the scene inside this context or retain
    RNA references after leaving it. Other object types retain the temporary
    conversion path and its cleanup. Borrowed mesh data is never cleared.
    """
    graph = depsgraph if depsgraph is not None else bpy.context.evaluated_depsgraph_get()
    evaluated = obj.evaluated_get(graph)
    if isinstance(evaluated.data, bpy.types.Mesh):
        yield evaluated.data
    else:
        mesh = evaluated.to_mesh()
        try:
            yield mesh
        finally:
            evaluated.to_mesh_clear()


@contextmanager
def evaluated_record_components(obj):
    """Keep native mesh and point-cloud record storage alive during one read.

    Providers may store records as legacy loose vertices, native points, or a
    mixture. Retain the GeometrySet owner: its component RNA must not escape
    this context or be used across another dependency-graph update.
    """
    graph = bpy.context.evaluated_depsgraph_get()
    evaluated = obj.evaluated_get(graph)
    geometry = evaluated.evaluated_geometry()
    with evaluated_mesh(obj, depsgraph=graph) as mesh:
        yield mesh, geometry.pointcloud


class GeometryEvidence:
    def read(self, obj):
        with evaluated_mesh(obj) as mesh:
            attrs = mesh.attributes
            count = len(mesh.polygons)
            def values(key, default=0):
                attr = attrs.get(key)
                if attr is None or attr.domain != 'FACE':
                    return np.full(count, default, dtype=np.int32)
                result = np.empty(count, dtype=np.bool_ if attr.data_type=='BOOLEAN' else np.int32)
                attr.data.foreach_get('value', result)
                return result
            # Fetch each native attribute once. Per-face RNA lookups made the
            # coverage report itself expensive on otherwise valid large meshes.
            diagnostic = values('st_diagnostic').astype(bool)
            surface = ~diagnostic
            paths = values('st_path_id', -3)
            kinds = values('st_component_kind')
            junction_ids = values('st_intersection_id', -1)
            sizes = np.empty(count, dtype=np.int32)
            mesh.polygons.foreach_get('loop_total', sizes)
            sections = np.unique(paths[surface & (paths>=0)]).tolist()
            junctions = np.unique(junction_ids[surface & (paths==-1) & (kinds==0)]).tolist()
            transition_ids = np.unique(values('st_transition_id', -1)[surface & values('st_transition').astype(bool)]).tolist()
            diagnostics = set(zip(*(values(key, -1)[diagnostic].tolist() for key in
                ('st_section_id','st_intersection_id','st_bundle_id','st_diagnostic_code'))))
            return dict(section_ids=sections, intersection_ids=junctions, transition_ids=transition_ids,
                diagnostics=[dict(section_id=s,intersection_id=j,bundle_id=b,code=c)
                    for s,j,b,c in sorted(diagnostics)], quads=int(np.sum(surface & (sizes==4))),
                other_faces=int(np.sum(surface & (sizes!=4))),
                vertices=len(mesh.vertices))
