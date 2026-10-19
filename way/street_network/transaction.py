"""Own only data created synchronously by one import/component transaction."""
import bpy


class ImportTransaction:
    # Node groups/materials come from the shared asset service and are reusable
    # after a failed import. Source/generated geometry and reports are owned here.
    kinds = ('objects', 'collections', 'meshes', 'curves', 'texts')

    def __enter__(self):
        self.before = {kind: {value.as_pointer() for value in getattr(bpy.data, kind)}
                       for kind in self.kinds}
        return self

    def __exit__(self, error_type, error, traceback):
        if error_type is not None:
            for kind in self.kinds:
                collection = getattr(bpy.data, kind)
                for value in tuple(collection):
                    if value.as_pointer() not in self.before[kind]:
                        collection.remove(value, do_unlink=True)
        return False
