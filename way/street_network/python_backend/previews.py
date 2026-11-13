"""Explicitly owned, provisional street surfaces for a first publication.

Only completed top-level Parts enter this collection. Their geometry has local
guards, but complete-network and terrain validation still precede acceptance.
The final publisher reuses the same object and mesh IDs, never duplicate output.
"""
import json
import bpy
import numpy as np
from .mesh import Mesh, publication_coordinates
from .providers import Part
from .uv import PlanarMapping, apply_mapping, publication_uv_coordinates
from .materials import resolve, assign
from ....progressive.ownership import TRANSIENT, alive


class PendingPreviews:
    batch_size = 8

    def __init__(self, controller, writer, ownership, on_preview, scale):
        self.controller = controller
        self.writer = writer
        self.ownership = ownership
        self.on_preview = on_preview
        self.scale = scale
        self.collection = None
        self.parent = controller.users_collection[0]
        self.objects = {}
        self.meshes = {}
        self.states = {}
        self.modifiers = {}
        self.waiting = []
        self.visible = 0
        self.changed = 0
        self.finished = False

    @staticmethod
    def _state(obj):
        return (tuple(v for row in obj.matrix_basis for v in row),
                tuple(v for row in obj.matrix_parent_inverse for v in row),
                obj.hide_viewport, obj.hide_render, obj.hide_get(), obj.mode)

    @staticmethod
    def _modifiers(obj):
        result = []
        for modifier in obj.modifiers:
            group = getattr(modifier, 'node_group', None)
            inputs = []
            if group is not None:
                for socket in group.interface.items_tree:
                    if socket.item_type == 'SOCKET' and socket.in_out == 'INPUT':
                        item = getattr(modifier.properties.inputs, socket.identifier, None)
                        if item is not None and hasattr(item, 'value'):
                            value = item.value
                            if not isinstance(value, (str, int, float, bool, type(None), bpy.types.ID)):
                                value = tuple(value)
                            inputs.append((socket.identifier, getattr(item, 'type', None), value))
            result.append((modifier.as_pointer(), modifier.show_viewport, modifier.show_render,
                           group, tuple(inputs)))
        return tuple(result)

    def check(self, *, full=False):
        from .scene import StreetInputChanged
        if self.collection is not None:
            if (not alive(self.collection) or not alive(self.parent)
                    or self.collection.name not in self.parent.children):
                raise StreetInputChanged('Pending street preview collection was deleted or unlinked')
            if len(self.collection.objects) != self.visible:
                raise StreetInputChanged('Pending street preview collection membership changed')
        if not full:
            return
        waiting = set(self.waiting)
        for key, obj in self.objects.items():
            self.check_object(key)
            if obj not in waiting and self.collection not in obj.users_collection:
                raise StreetInputChanged('Pending street preview was unlinked')

    def check_object(self, key):
        from .scene import StreetInputChanged
        obj = self.objects[key]
        if (not alive(obj) or not alive(self.meshes[key]) or obj.data != self.meshes[key]
                or obj.parent != self.controller or obj.get('st_py_controller') != self.controller
                or obj.get('st_py_key') != key.text() or not obj.get('st_py_preview')
                or self._state(obj) != self.states[key] or self._modifiers(obj) != self.modifiers[key]):
            raise StreetInputChanged('Pending street preview ownership changed')

    def add(self, key, product):
        """Allocate one complete Part; link only at a recorded visible batch."""
        if not isinstance(product, Part) or product.mesh is None:
            return False
        if key in self.objects:
            raise ValueError('Repeated pending street product: ' + key.text())
        mesh = apply_mapping(product.mesh, PlanarMapping())
        stored = Mesh(publication_coordinates(mesh.xyz, self.scale, physical=True), mesh.quads)
        stored.validate()
        data = self.writer._mesh(key.text(), mesh, self.scale, ownership=self.ownership)
        obj = self.ownership.own('objects', bpy.data.objects.new(key.text(), data))
        self.objects[key], self.meshes[key] = obj, data
        assign(data, resolve(self.controller, mesh.material, self.writer.materials))
        obj.parent = self.controller
        obj['st_py_controller'] = self.controller
        obj['st_py_key'] = key.text()
        obj['st_py_owner'] = product.owner.text()
        obj['st_py_preview'] = True
        self.states[key] = self._state(obj)
        self.modifiers[key] = self._modifiers(obj)
        self.waiting.append(obj)
        return self.flush() if self.visible == 0 or len(self.waiting) >= self.batch_size else False

    def flush(self):
        if not self.waiting and not self.changed:
            return False
        if self.waiting:
            self.check(full=True)
            if self.collection is None:
                self.collection = self.ownership.own('collections',
                    bpy.data.collections.new(self.controller.name + ' / pending streets'))
                self.collection['st_py_preview'] = True
                self.collection[TRANSIENT] = True
                self.parent.children.link(self.collection)
            for obj in self.waiting:
                self.collection.objects.link(obj)
                self.visible += 1
            self.waiting.clear()
        self.changed = 0
        self.on_preview(self.visible)
        return True

    def refreshed(self, key):
        """Notify the session before a turn can expose an in-place refinement."""
        self.modifiers[key] = self._modifiers(self.objects[key])
        self.changed += 1
        return self.changed >= self.batch_size

    def validate_final(self, products, raw_products):
        """Match live buffers to the already validated publication payloads.

        Preview objects are accessible to user edits during earlier UI turns.
        This readback closes that gap after their last yield, before acceptance.
        """
        from .scene import StreetInputChanged
        from .terrain.blender import RAW_Z
        self.check(full=True)
        for key, obj in self.objects.items():
            payload = products[key].mesh
            mesh = obj.data
            xyz = np.empty(len(mesh.vertices) * 3, dtype=np.float32)
            ids = np.empty(len(mesh.loops), dtype=np.int32)
            totals = np.empty(len(mesh.polygons), dtype=np.int32)
            mesh.vertices.foreach_get('co', xyz)
            mesh.loops.foreach_get('vertex_index', ids)
            mesh.polygons.foreach_get('loop_total', totals)
            material = resolve(self.controller, payload.material, self.writer.materials)
            if (not np.array_equal(xyz, publication_coordinates(payload.xyz, self.scale).ravel())
                    or not np.array_equal(ids, payload.quads.ravel())
                    or len(totals) != len(payload.quads) or np.any(totals != 4)
                    or set(layer.name for layer in mesh.uv_layers) != set(payload.uv_layers)
                    or tuple(mesh.materials) != (() if material is None else (material,))
                    or obj.get('st_py_owner') != products[key].owner.text()
                    or obj.get('st_py_metadata') != json.dumps(products[key].metadata,
                        separators=(',', ':'), allow_nan=False)):
                raise StreetInputChanged('Pending street preview geometry or materials changed')
            for name, expected in payload.uv_layers.items():
                values = np.empty(expected.size, dtype=np.float32)
                mesh.uv_layers[name].uv.foreach_get('vector', values)
                if not np.array_equal(values, publication_uv_coordinates(expected).ravel()):
                    raise StreetInputChanged('Pending street preview UV coordinates changed')
            raw = mesh.attributes.get(RAW_Z)
            if raw is None or raw.domain != 'POINT' or raw.data_type != 'FLOAT':
                raise StreetInputChanged('Pending street preview terrain coordinates changed')
            values = np.empty(len(mesh.vertices), dtype=np.float32)
            raw.data.foreach_get('value', values)
            if not np.array_equal(values, np.asarray(raw_products[key].mesh.xyz[:, 2], dtype=np.float32)):
                raise StreetInputChanged('Pending street preview terrain coordinates changed')

    def promote(self):
        """Remove provisional links/tags inside the non-yielding final commit."""
        for obj in self.objects.values():
            if 'st_py_preview' in obj:
                del obj['st_py_preview']
        if self.collection is not None:
            pointer = self.collection.as_pointer()
            bpy.data.collections.remove(self.collection, do_unlink=True)
            self.ownership.pending.pop(pointer, None)
            parent = self.ownership.parent
            if parent is not None:
                parent.pending.pop(pointer, None)
            self.collection = None

    def accept(self):
        self.finished = True
        if self.visible:
            self.visible = 0
            self.on_preview(0)

    def rollback(self):
        if self.finished:
            return
        self.finished = True
        visible = self.visible
        try:
            self.ownership.rollback()
        finally:
            self.visible = 0
            if visible:
                self.on_preview(0)
