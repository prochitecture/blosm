"""Respect Blender's network preference without requiring bpy in the worker."""
import sys


class OnlineAccessError(RuntimeError):
    """A download was prevented by Blender's online-access preference."""


def require_online_access():
    bpy = sys.modules.get('bpy')
    if bpy is not None and not bpy.app.online_access:
        raise OnlineAccessError(
            'Online access is disabled in Blender. Enable Allow Online Access '
            'in Preferences > System to download data.')
