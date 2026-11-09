"""Compatibility boundary for Blender's internal glTF importer.

Both entry paths create Blender data synchronously on Blender's main thread.
A source path supplies naming and relative resources; it need not exist when
content supplies a self-contained GLB/glTF.
"""

from dataclasses import dataclass
from os import fspath
import logging
import struct
import sys

import bpy


@dataclass(frozen=True)
class ImportResult:
    copyright: str
    messages: tuple


def _import_settings():
    # Read the installed operator's defaults without executing it. This also
    # includes settings added by newer Blender releases.
    import io_scene_gltf2

    settings = {}
    for prop in bpy.ops.import_scene.gltf.get_rna_type().properties:
        if prop.type in {'BOOLEAN', 'INT', 'FLOAT', 'STRING', 'ENUM'}:
            settings[prop.identifier] = (
                tuple(prop.default_array) if getattr(prop, 'is_array', False)
                else set(prop.default_flag) if getattr(prop, 'is_enum_flag', False)
                else prop.default
            )
    if 'import_scene_as_collection' in settings:
        settings['import_scene_as_collection'] = True
    if hasattr(io_scene_gltf2, 'set_debug_log'):
        settings['loglevel'] = io_scene_gltf2.set_debug_log()
    else:
        # Blender 3.6-4.2 kept this method on the operator itself.
        settings['loglevel'] = {
            0: logging.CRITICAL, 1: logging.ERROR,
            2: logging.WARNING, 3: logging.INFO,
        }.get(bpy.app.debug_value, logging.NOTSET)
    # Preserve the user-extension discovery performed by the glTF operator.
    settings['import_user_extensions'] = []
    for addon_name in bpy.context.preferences.addons.keys():
        module = sys.modules.get(addon_name)
        extension = getattr(module, 'glTF2ImportUserExtension', None)
        if extension is not None:
            settings['import_user_extensions'].append(extension())
    return settings


def import_gltf(filepath, *, content=None, report=None, coordinate_offset=None):
    """Import a file, or bytes with filepath as their resource base path.

    Return attribution and importer messages without retaining decoded buffers.
    Errors propagate to the owning Blosm operation without per-tile UI reports.
    Optional report forwards messages from successful imports using the
    report({level}, message) signature of a Blender operator. Like Blender's
    importer, this does not roll back partially created objects.
    External resources retain ordinary file-resolution behavior.
    """
    # Blender 4.3 renamed the scene
    # construction modules. Keep both call paths explicit at this boundary.
    version = bpy.app.version
    if (version[0] == 3 and version[1] == 6) or (version[0] == 4 and version[1] <= 2):
        from io_scene_gltf2.blender.imp.gltf2_blender_gltf import BlenderGlTF
    elif (version[0] == 4 and version[1] >= 3) or version[0] > 4:
        from io_scene_gltf2.blender.imp.blender_gltf import BlenderGlTF
    else:
        raise RuntimeError('Unsupported Blender version for direct 3D Tiles import')

    from io_scene_gltf2.io.imp.gltf2_io_gltf import glTFImporter, ImportError
    from io_scene_gltf2.io.com.gltf2_io import gltf_from_dict

    filepath = fspath(filepath)
    settings = _import_settings()
    if coordinate_offset is not None:
        settings['import_select_created_objects'] = False
    settings['filepath'] = filepath
    importer = glTFImporter(filepath, settings)
    try:
        if content is None:
            importer.read()
        else:
            content = memoryview(content).cast('B')
            try:
                if content[:4] == b'glTF':
                    document, importer.glb_buffer = importer.load_glb(content)
                else:
                    document = importer.load_json(content)
                importer.check_version(document)
                importer.data = gltf_from_dict(document)
            except (AssertionError, struct.error) as exc:
                raise ImportError('Could not parse the glTF tile: invalid or truncated data') from exc

        importer.checks()
        if coordinate_offset is not None and bpy.app.debug_value != 100:
            # Rebase the parsed, private node transforms in double precision,
            # before Blender converts them to single-precision mathutils values.
            # glTF is Y-up: (x, y, z) becomes Blender's (x, -z, y).
            # This replaces the old process-wide glTF class monkey-patches.
            x, y, z = coordinate_offset
            for node in importer.data.nodes or ():
                if node.matrix is not None:
                    node.matrix[12] -= x
                    node.matrix[13] -= z
                    node.matrix[14] += y
                else:
                    translation = node.translation or [0., 0., 0.]
                    node.translation = [translation[0] - x, translation[1] - z, translation[2] + y]
        BlenderGlTF.create(importer)
        # Evaluate newly created objects before Blosm reads matrix_world for
        # RTC and tile transforms; otherwise these matrices can be stale.
        bpy.context.view_layer.update()
        messages = tuple(importer.log.messages()) if hasattr(importer.log, 'messages') else ()
        if report is not None:
            for level, message in messages:
                report({level}, message)
        return ImportResult(importer.data.asset.copyright or '', messages)
    finally:
        # Failed tiles are handled by Blosm's manager: details go to the System
        # Console, with one aggregate UI report from the outer import operator.
        # Per-import handlers use shared loggers; remove them on failure too.
        if hasattr(importer, 'log_handler'):
            importer.log.removeHandler(importer.log_handler)
        else:
            importer.log.flush()
