"""Shared cooperative import scheduling and background I/O.

core and transport have no Blender dependencies. blender owns the UI lifecycle;
data-type adapters own traversal, parsing, scene changes and result reporting.
"""
