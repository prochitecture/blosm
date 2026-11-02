"""Verified GN library loading, with shared dependencies across entry points.

Blender's native append reuse covers one source library. Reconciliation also
covers groups built in the current file and identical helpers in other published
libraries. Original names come from Blender's weak references and the asset
loader's existing provenance property, never from stripping a numeric suffix.
No additional PADL persistence fields are needed.
"""
from pathlib import Path
import hashlib
import json

from ..padl.errors import PadlError
from .node_group_identity import graph_fingerprint


def original_name(group):
    # Blender clears a weak reference when an ID is renamed. The ordinary
    # loader provenance remains stable through renaming and saving/reopening.
    try:
        reference = json.loads(group.get('blosm_asset_reference', 'null'))
        if isinstance(reference, list) and len(reference) in (5, 6) and reference[2] == 'node_groups' and isinstance(reference[3], str):
            return reference[3]
    except (TypeError, ValueError):
        pass
    reference = group.library_weak_reference
    return reference.id_name[2:] if reference else group.name


def _storage(group):
    import bpy
    if group.library:
        return str(Path(bpy.path.abspath(group.library.filepath)).resolve())
    return None


def _closure(roots):
    found = {}
    pending = list(roots)
    while pending:
        group = pending.pop()
        if group.as_pointer() in found:
            continue
        found[group.as_pointer()] = group
        pending.extend(node.node_tree for node in group.nodes
                       if getattr(node, 'node_tree', None) is not None)
    return found


def _metadata(value):
    """Conservative equality for behavior beyond the published graph digest."""
    import bpy
    if isinstance(value, bpy.types.Text):
        return ('TEXT', value.as_string())
    if isinstance(value, bpy.types.ID):
        # Never alias different scene objects/materials merely by their names.
        return ('ID', value.as_pointer())
    if hasattr(value, 'items'):
        return tuple(sorted((key, _metadata(item)) for key, item in value.items()
                            if key not in {'blosm_asset_reference', '_RNA_UI'}))
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return tuple(_metadata(item) for item in value)


class VerifiedNodeGroupLibrary:
    """Share unchanged groups; refuse conflicts without changing existing data.

    A load is synchronous. All graph and metadata checks finish before any
    remapping. Rollback removes only datablocks allocated by this attempt.
    Linked groups retain their library ownership; explicit linked/local copies
    are intentionally separate storage choices.
    """

    def __init__(self, loader, *, error_type=PadlError):
        self.loader = loader
        self.error_type = error_type

    def load(self, reference, fingerprint, *, validate=None):
        return self.load_many({reference: fingerprint}, validate=validate)[reference]

    def load_many(self, expected, *, validate=None):
        import bpy
        if bpy.app.version < (5, 2, 0):
            raise self.error_type('PADL GN loading requires Blender 5.2+', code='unsupported_runtime')
        expected = dict(expected)
        if any(reference.kind != 'node_groups' for reference in expected):
            raise ValueError('VerifiedNodeGroupLibrary only loads node groups')

        # Even an already loaded root must pass edition and library-byte checks.
        storage, sources = {}, {}
        for reference in expected:
            linked = self.loader.uses_link(reference.library)
            spec, path = self.loader.resolve(reference.library, link=linked)
            storage[reference] = str(path) if linked else None
            sources[reference] = (spec, path)

        before = set(bpy.data.user_map())
        existing = list(bpy.data.node_groups)
        candidates = {}
        for group in existing:
            candidates.setdefault((_storage(group), original_name(group)), []).append(group)
        memo = {}
        result = {}

        def check_root(reference, group):
            if graph_fingerprint(group, memo) != expected[reference]:
                raise self.error_type(f'Node group {reference.name!r} or a dependency differs from its published fingerprint',
                                code='asset_conflict')
            if validate is not None:
                validate(reference, group)

        # Cached roots and groups imported with a wrapper are equally eligible.
        # Reading live datablocks also handles new providers, renames and reopen.
        for reference in expected:
            matches = candidates.get((storage[reference], reference.name), [])
            if matches:
                for group in matches:
                    check_root(reference, group)
                result[reference] = matches[0]

        try:
            missing = [reference for reference in expected if reference not in result]
            if missing:
                result.update(self.loader.load_many(missing, reuse_local_id=True))
            # Native reuse can return edited dependencies. Verify the whole root,
            # even when the root itself was just appended from trusted bytes.
            memo.clear()
            for reference, group in result.items():
                check_root(reference, group)

            roots = _closure(result.values())
            signatures = {}

            def signature(group):
                pointer = group.as_pointer()
                if pointer in signatures:
                    return signatures[pointer]
                if group.animation_data is not None or group.override_library is not None:
                    raise self.error_type(f'Animated/overridden GN asset {group.name!r} cannot be reused as an immutable generator',
                                    code='asset_conflict')
                # Graph fingerprints intentionally omit editor layout and custom
                # metadata. Include metadata and external ID references here so
                # structurally similar but behaviorally different groups cannot
                # become shared accidentally.
                pointers = []
                children = []
                node_properties = []
                for node in group.nodes:
                    node_properties.append((node.name, _metadata(node)))
                    child = getattr(node, 'node_tree', None)
                    if child is not None:
                        children.append(signature(child))
                    for prop in node.bl_rna.properties:
                        if prop.type == 'POINTER' and prop.identifier not in {'rna_type', 'node_tree', 'parent'}:
                            value = getattr(node, prop.identifier, None)
                            if isinstance(value, bpy.types.ID):
                                pointers.append((node.name, prop.identifier, _metadata(value)))
                    for socket in (*node.inputs, *node.outputs):
                        value = getattr(socket, 'default_value', None)
                        if isinstance(value, bpy.types.ID):
                            pointers.append((node.name, socket.identifier, _metadata(value)))
                for socket in group.interface.items_tree:
                    value = getattr(socket, 'default_value', None)
                    if isinstance(value, bpy.types.ID):
                        pointers.append(('interface', socket.identifier, _metadata(value)))
                value = (graph_fingerprint(group, memo), _metadata(dict(group.items())), tuple(node_properties),
                         tuple(pointers), tuple(children))
                value = hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()
                signatures[pointer] = value
                return value

            replacements = {}
            for group in roots.values():
                signature(group)
                if group.library:
                    continue
                key = (None, original_name(group))
                matches = candidates.setdefault(key, [])
                if group not in matches:
                    matches.append(group)
                for candidate in matches:
                    if signature(candidate) != signature(group):
                        raise self.error_type(f'Conflicting shared GN group {key[1]!r}; existing graph was edited or uses a different revision',
                                        code='asset_conflict')
                # Prefer an existing canonical ID. Also repair verified copies
                # left by the previous append-only loading path. Hand-made
                # copies without source identity are not inferred from suffixes.
                canonical = min(matches, key=lambda candidate: (candidate not in before, candidate.name != key[1]))
                for candidate in matches:
                    if candidate != canonical:
                        replacements[candidate] = canonical

            # Use the existing Blosm provenance layout on dependencies too.
            # This is needed because native weak references disappear on rename.
            provenance = {}
            for reference, root in result.items():
                spec, path = sources[reference]
                for group in _closure([root]).values():
                    if group.library:
                        continue
                    name = original_name(group)
                    identity = (str(path), spec.revision, 'node_groups', name, False)
                    provenance[group] = json.dumps(identity + (spec.sha256,) if spec.sha256 else identity)

            # No validation after this point can reject a graph. Only verified
            # equivalents of the same source group are replaced. Remapping all
            # users preserves socket bindings and also repairs older duplicates.
            for duplicate, canonical in replacements.items():
                canonical.use_fake_user |= duplicate.use_fake_user
                duplicate.user_remap(canonical)
            for group, key in provenance.items():
                group = replacements.get(group, group)
                if not group.get(self.loader.provenance_key):
                    group[self.loader.provenance_key] = key
            result = {reference: replacements.get(group, group) for reference, group in result.items()}
            if replacements:
                bpy.data.batch_remove(ids=list(replacements))
            for group in result.values():
                if not group.library:
                    group.use_fake_user = True
            return result
        except Exception:
            created = set(bpy.data.user_map()) - before
            if created:
                bpy.data.batch_remove(ids=created)
            raise
