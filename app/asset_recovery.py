"""Restore a saved retained-library reference without choosing a new revision.

The digest directory is part of ImmutableLibraryStore's persistent contract.
Inspection does not hash files or load Blender data. Restoration streams bytes
through the same verified, atomic, no-replacement store as normal imports.
"""
from dataclasses import dataclass
from pathlib import Path

from .library_store import AssetStoreError, ImmutableLibraryStore, validate_digest


@dataclass(frozen=True)
class RetainedLibraryReference:
    name: str
    path: Path
    digest: str

    @classmethod
    def inspect(cls, name, absolute_path):
        path = Path(absolute_path)
        if not path.is_absolute() or path.name != 'library.blend':
            return None
        try:
            digest = validate_digest(path.parent.name)
        except ValueError:
            return None
        return cls(name, path, digest)

    def restore(self, source, *, root=None):
        """Keep the original location by default; allow an explicit new store."""
        store = ImmutableLibraryStore(self.path.parent.parent if root is None else root)
        return store.install(source, self.digest)


def inspect_library(library):
    import bpy
    return RetainedLibraryReference.inspect(library.name,
        bpy.path.abspath(library.filepath, library=library.parent))


def restore_library(library, source, *, root=None, expected_digest=None):
    """Restore exact saved bytes and reload their existing Blender references.

    No append, catalog upgrade, scene reopening or save is performed. Existing
    corrupt files are never overwritten. A reload error leaves the verified
    copy on disk and restores the previous filepath; the caller reports it.
    """
    reference = inspect_library(library)
    if reference is None:
        raise AssetStoreError('This is not a retained library reference')
    if expected_digest is not None and expected_digest != reference.digest:
        raise AssetStoreError('The saved library revision changed while choosing its replacement')
    target = reference.restore(source, root=root)
    previous_path = library.filepath
    try:
        library.filepath = str(target)
        library.reload()
        missing = [block.name for block in library.users_id if block.is_missing]
        if missing:
            raise RuntimeError('Linked datablocks are still missing: '+', '.join(missing[:5]))
    except (RuntimeError, OSError) as exc:
        library.filepath = previous_path
        raise AssetStoreError(f'Exact revision restored to {target}, but Blender could not reload it: {exc}') from exc
    return target
