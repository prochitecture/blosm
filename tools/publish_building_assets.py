"""Publish the runtime catalog and Python document from a precompiled .blend.

Run in Blender: --background LIBRARY --python this_file -- [--tree NAME].
Only sidecars are written; the compiled library is not modified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import ModuleType


def main():
    import bpy
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default='PADL simple_3d')
    parser.add_argument('--role', default='simple_3d')
    args = parser.parse_args(sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else [])
    root = Path(__file__).resolve().parents[1]
    # Import services without registering the addon or its UI.
    for name, path in (('_blosm_publish', root), ('_blosm_publish.app', root / 'app'),
                       ('_blosm_publish.padl', root / 'padl')):
        module = ModuleType(name)
        module.__path__ = [str(path)]
        sys.modules[name] = module
    from _blosm_publish.app.node_group_identity import graph_fingerprint
    from _blosm_publish.padl.persistence import load_program
    library = Path(bpy.data.filepath)
    if library.parent.resolve() != (root / 'assets').resolve():
        parser.error('Open the compiled library in blosm/assets before publishing')
    tree = bpy.data.node_groups[args.tree]
    program = load_program(tree)
    if (tree.get('padl_compiler_version') != '1' or tree.get('padl_schema_version') != program.schema_version or
            tree.get('padl_program_sha256') != program.fingerprint):
        parser.error('Compiled document fingerprint differs from the saved document')
    digest = hashlib.sha256(library.read_bytes()).hexdigest()
    authored = json.loads(library.with_suffix('.assets.json').read_text(encoding='utf-8'))
    if authored.get('schema') != 'blosm-padl-roofs-v1' or authored.get('sha256') != digest:
        parser.error('Rebuild the compiler asset catalog for this library before publishing')
    closure, pending = set(), [tree]
    while pending:
        group = pending.pop()
        if group.name in closure:
            continue
        closure.add(group.name)
        pending.extend(node.node_tree for node in group.nodes if getattr(node, 'node_tree', None))
    available = sorted(roof for roof, record in authored['groups'].items() if record['name'] in closure)
    fingerprint = graph_fingerprint(tree)
    editions = [record.get('edition', 'BASE') for record in authored.get('programs', [])
                if record.get('fingerprint') == fingerprint]
    if not editions:
        parser.error('Compiled program is absent from the matching compiler catalog')
    document = library.with_name(library.stem + '.' + args.role + '.json')
    document.write_text(json.dumps(program.document, indent=2) + '\n', encoding='utf-8')
    # Keep this separate from the compiler's roof-provider .assets.json catalog.
    catalog_path = library.with_suffix('.runtime.json')
    prior = json.loads(catalog_path.read_text()) if catalog_path.exists() else {'programs': {}}
    if not isinstance(prior.get('programs'), dict):
        parser.error('Existing runtime catalog has an incompatible program layout')
    prior.update(schema='blosm.buildings.v1', library={'file': library.name, 'sha256': digest})
    prior['programs'][args.role] = {'tree': tree.name, 'document': document.name,
        'program_sha256': program.fingerprint, 'fingerprint': fingerprint,
        'available_roofs': available, 'edition': 'PRO' if 'PRO' in editions else 'BASE'}
    catalog_path.write_text(json.dumps(prior, indent=2) + '\n', encoding='utf-8')
    print('Published:', catalog_path)


if __name__ == '__main__':
    main()
