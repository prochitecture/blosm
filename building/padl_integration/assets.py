"""Versioned building-program roles, independent of Blender group names."""
import json
from pathlib import Path

from ...app.asset_library import AssetCatalog, AssetRef, BlenderAssetLoader, LibrarySpec
from ...app.node_group_library import VerifiedNodeGroupLibrary
from ...padl.errors import PadlError
from ...padl.parser import PadlParser


class BuildingProgramLibrary:
    def __init__(self, assets_root=None, *, app=None, catalog='prochitecture_buildings.runtime.json'):
        self.root = Path(assets_root) if assets_root else Path(__file__).resolve().parents[2] / 'assets'
        self.catalog = json.loads((self.root / catalog).read_text(encoding='utf-8'))
        if self.catalog.get('schema') != 'blosm.buildings.v1':
            raise PadlError('Unsupported building asset catalog', code='incompatible_asset')
        self.app = app
        self._loader = None

    def record(self, role):
        try:
            record = self.catalog['programs'][role]
        except KeyError as exc:
            raise PadlError(f'No building program for {role}', code='missing_asset') from exc
        if record.get('edition', 'BASE') not in {'BASE', 'PRO'}:
            raise PadlError('Building program edition must be BASE or PRO', code='incompatible_asset')
        if record.get('edition', 'BASE') == 'PRO' and not getattr(self.app, 'isPro', False):
            raise PadlError(f'Building program {role} requires Blosm Pro', code='unavailable_program')
        return record

    def program(self, role):
        record = self.record(role)
        path = (self.root / record['document']).resolve()
        if not path.is_relative_to(self.root.resolve()):
            raise PadlError('Building document path escapes assets', code='incompatible_asset')
        program = PadlParser().load(path)
        if program.fingerprint != record['program_sha256']:
            raise PadlError('Building document differs from the published program', code='program_mismatch')
        return program

    def available_roofs(self, role):
        roofs = self.record(role).get('available_roofs')
        if not isinstance(roofs, list) or any(roof not in {'flat', 'skillion', 'gabled', 'hipped'} for roof in roofs):
            raise PadlError('Republish the runtime catalog with its roof capabilities', code='incompatible_asset')
        return frozenset(roofs)

    def tree(self, role):
        import bpy
        if bpy.app.version < (5, 2, 0):
            raise PadlError('Geometry Nodes buildings require Blender 5.2 or later', code='unsupported_runtime')
        record = self.record(role)
        library = self.catalog['library']
        spec = LibrarySpec('prochitecture_buildings', library['file'], revision=library['sha256'],
                           sha256=library['sha256'], self_contained=True,
                           edition=record.get('edition', 'BASE'))
        if self._loader is None:
            catalog = AssetCatalog([spec])
            self._loader = (BlenderAssetLoader.for_app(self.app, catalog) if self.app else
                            BlenderAssetLoader(self.root, catalog))
        program = self.program(role)

        def validate(reference, tree):
            from ...padl.persistence import load_program
            saved = load_program(tree)
            if (tree.get('padl_compiler_version') != '2' or
                    tree.get('padl_schema_version') != program.schema_version or
                    tree.get('padl_program_sha256') != program.fingerprint or
                    saved.fingerprint != program.fingerprint):
                raise PadlError('Compiled building asset identity is incompatible', code='program_mismatch')
            sockets = [item for item in tree.interface.items_tree if item.item_type == 'SOCKET']
            inputs = {item.name for item in sockets if item.in_out == 'INPUT'}
            outputs = {item.name for item in sockets if item.in_out == 'OUTPUT'}
            if inputs != {'Geometry', 'Seed', 'Terrain', 'Terrain limits'} or not {'Geometry', 'Status'} <= outputs:
                raise PadlError('Compiled building asset has an incompatible public interface', code='incompatible_asset')

        return VerifiedNodeGroupLibrary(self._loader).load(
            AssetRef(spec.key, 'node_groups', record['tree']), record['fingerprint'], validate=validate)
