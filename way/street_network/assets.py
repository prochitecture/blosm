"""Street materials loaded through the feature-independent asset service."""
from dataclasses import replace
from ...app.asset_library import AssetCatalog, AssetRef, LibrarySpec

MATERIAL_LIBRARY = LibrarySpec('streets.materials', 'prochitecture_street_materials.blend', revision='1.0')
PAVEMENT = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Pavement')
GRASS = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Grass')
DIAGNOSTICS = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Diagnostics')
DIAGNOSTIC_RED = AssetRef(MATERIAL_LIBRARY.key, 'materials', 'Prochitecture Streets / Diagnostic red')


def material_catalog(external_path=None):
    library = replace(MATERIAL_LIBRARY, path=external_path, source='external') if external_path else MATERIAL_LIBRARY
    return AssetCatalog((library,))
