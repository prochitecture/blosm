"""Recovery controls shared by every retained Blender asset feature."""
from types import SimpleNamespace
from pathlib import Path

import bpy
from ..app.asset_library import AssetError, BlenderAssetLoader
from ..app.asset_recovery import inspect_library, restore_library
from ..app.library_store import AssetStoreError


def retained_libraries():
    for library in bpy.data.libraries:
        reference = inspect_library(library)
        if reference is not None:
            yield library, reference


class BLOSM_OT_RestoreRetainedLibrary(bpy.types.Operator):
    bl_idname = 'blosm.restore_retained_library'
    bl_label = 'Restore retained library'
    bl_description = 'Choose a copy of the exact saved revision, verify it and reload the linked assets'
    bl_options = {'INTERNAL'}

    library_name: bpy.props.StringProperty(options={'HIDDEN', 'SKIP_SAVE'})
    expected_digest: bpy.props.StringProperty(options={'HIDDEN', 'SKIP_SAVE'})
    filepath: bpy.props.StringProperty(name='Library copy', subtype='FILE_PATH')
    filter_glob: bpy.props.StringProperty(default='*.blend', options={'HIDDEN'})
    use_current_directory: bpy.props.BoolProperty(name='Use current retained directory', default=False,
        description='Restore into the directory from the import settings and relink this project; otherwise keep the saved location')

    def invoke(self, context, event):
        library = bpy.data.libraries.get(self.library_name)
        reference = inspect_library(library) if library is not None else None
        if reference is None:
            self.report({'ERROR'}, 'The retained library reference is no longer available')
            return {'CANCELLED'}
        self.expected_digest = reference.digest
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def draw(self, context):
        self.layout.label(text='Choose the exact saved library revision')
        self.layout.prop(self, 'use_current_directory')
        if self.use_current_directory:
            self.layout.prop(context.scene.blosm, 'assetLibraryDirectory')

    def execute(self, context):
        library = bpy.data.libraries.get(self.library_name)
        if library is None:
            self.report({'ERROR'}, 'The retained library reference is no longer available')
            return {'CANCELLED'}
        try:
            root = None
            if self.use_current_directory:
                settings = SimpleNamespace(basePath=str(Path(__file__).resolve().parents[1]), assetLibraryMode='RETAINED',
                    assetLibraryDirectory=context.scene.blosm.assetLibraryDirectory)
                root = BlenderAssetLoader.for_app(settings).linked_store.root
            restore_library(library, bpy.path.abspath(self.filepath), root=root,
                            expected_digest=self.expected_digest or None)
        except (AssetError, AssetStoreError, OSError, ValueError) as exc:
            self.report({'ERROR'}, str(exc))
            return {'CANCELLED'}
        self.report({'INFO'}, 'Exact library revision restored and reloaded; save the project to retain any changed link')
        return {'FINISHED'}


class BLOSM_PT_AssetLibraries(bpy.types.Panel):
    bl_label = 'Asset libraries'
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'Blosm'

    @classmethod
    def poll(cls, context):
        return any(retained_libraries())

    def draw(self, context):
        for library, reference in retained_libraries():
            box = self.layout.box()
            box.label(text=library.name)
            # Never hash large files while drawing a panel. Verify on demand.
            missing = not reference.path.is_file()
            box.label(text='Library file is missing' if missing else 'Library file is present',
                      icon='ERROR' if missing else 'FILE_BLEND')
            action = box.operator('blosm.restore_retained_library',
                                  text='Restore exact revision' if missing else 'Verify or restore revision')
            action.library_name = library.name


classes = (BLOSM_OT_RestoreRetainedLibrary, BLOSM_PT_AssetLibraries)
