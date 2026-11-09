"""Street materials loaded through the feature-independent asset service."""
from dataclasses import replace
from ...app.asset_library import AssetCatalog, AssetRef, LibrarySpec

MATERIAL_LIBRARY = LibrarySpec('streets.materials', 'prochitecture_street_materials.blend', revision='2.0')
PAVEMENT = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Pavement')
GRASS = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Grass')
DIAGNOSTICS = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Diagnostics')
DIAGNOSTIC_RED = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Diagnostic red')
ASPHALTS = tuple(AssetRef(MATERIAL_LIBRARY.key, 'materials',
    f'Prochitecture Streets / Asphalt / {wear} / {variant}')
    for wear in ('New', 'Worn', 'Cracked') for variant in ('PBR', 'Diffuse'))
TERRAIN_LIBRARY = LibrarySpec('streets.terrain', 'prochitecture_street_terrain.blend', revision='1.0')
TERRAIN_PROJECTION = AssetRef(TERRAIN_LIBRARY.key, 'node_groups', 'Prochitecture Streets / Direct terrain projection')


def terrain_catalog():
    return AssetCatalog((TERRAIN_LIBRARY,))


def material_catalog(external_path=None):
    library = replace(MATERIAL_LIBRARY, path=external_path, source='external') if external_path else MATERIAL_LIBRARY
    return AssetCatalog((library,))
