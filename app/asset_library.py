"""Feature-independent access to bundled and external Blender libraries.

Asset selection belongs to a feature's catalog, not this loader. Loading an
object/collection does not instantiate it in a scene. Direct callers default to
embedding; app callers use the selected storage policy. Reuse is by source,
revision, kind and original name, never a conflicting local datablock name.
"""
from dataclasses import dataclass
from pathlib import Path
import json
from .library_store import AssetStoreError, ImmutableLibraryStore, file_digest, validate_digest
from .asset_policy import AssetLoadPolicy


class AssetError(RuntimeError):
    pass


class AssetUnavailable(AssetError):
    pass


@dataclass(frozen=True)
class LibrarySpec:
    key: str
    path: str
    source: str = 'bundled'
    edition: str = 'BASE'
    revision: str = '1'
    sha256: str | None = None
    self_contained: bool = False

    def __post_init__(self):
        if not self.key or self.source not in ('bundled', 'external'):
            raise ValueError('Library needs a key and bundled/external source')
        if self.edition not in ('BASE', 'PRO'):
            raise ValueError('Library edition must be BASE or PRO')
        if self.sha256 is not None:
            validate_digest(self.sha256)


@dataclass(frozen=True)
class AssetRef:
    library: str
    kind: str
    name: str


class AssetCatalog:
    def __init__(self, libraries=()):
        self._libraries = {}
        for library in libraries:
            self.register(library)

    def register(self, library, *, replace=False):
        if library.key in self._libraries and not replace:
            raise ValueError('Duplicate asset library: ' + library.key)
        self._libraries[library.key] = library

    def get(self, key):
        try:
            return self._libraries[key]
        except KeyError:
            raise AssetUnavailable('Unknown asset library: ' + key) from None


class BlenderAssetLoader:
    """One service usable by streets, buildings, nature and future features."""
    provenance_key = 'blosm_asset_reference'

    def __init__(self, bundled_root, catalog=None, *, is_pro=False, linked_root=None, policy='APPEND'):
        self.bundled_root = Path(bundled_root).resolve()
        self.catalog = catalog if catalog is not None else AssetCatalog()
        self.is_pro = bool(is_pro)
        self._cache = {}
        self.policy = AssetLoadPolicy(policy)
        self.linked_store = ImmutableLibraryStore(linked_root) if linked_root is not None else None
        if self.policy.mode != 'APPEND' and self.linked_store is None:
            raise AssetError('Automatic or retained library loading requires a durable library directory')

    @classmethod
    def for_app(cls, app, catalog=None, *, linked_root=None):
        import bpy
        policy = AssetLoadPolicy(getattr(app, 'assetLibraryMode', 'AUTO'))
        if linked_root is None and policy.mode != 'APPEND':
            directory = getattr(app, 'assetLibraryDirectory', '')
            if directory:
                if directory.startswith('//') and not bpy.data.filepath:
                    raise AssetError('Save the scene before using a project-relative library directory, '
                                     'or choose an absolute directory')
                linked_root = Path(bpy.path.abspath(directory))
                if not linked_root.is_absolute():
                    raise AssetError('Choose an absolute or saved-project-relative library directory')
            else:
                # Resolve only. The store creates directories on first use;
                # an embedded/uncertified asset causes no disk allocation.
                resource_root = bpy.utils.resource_path('USER')
                if not resource_root:
                    raise AssetError('Choose a durable Blender library directory in the import settings')
                linked_root = Path(resource_root)/'datafiles'/'blosm'/'libraries'
        return cls(Path(app.basePath) / 'assets', catalog,
                   is_pro=getattr(app, 'isPro', False), linked_root=linked_root, policy=policy.mode)

    def uses_link(self, key, override=None):
        return self.policy.uses_link(self.catalog.get(key), override)

    def resolve(self, key, *, link=None):
        spec = self.catalog.get(key)
        link = self.uses_link(key, link)
        if spec.edition == 'PRO' and not self.is_pro:
            raise AssetUnavailable(f'Library {key!r} requires Blosm Pro')
        if link and spec.sha256 is not None and self.linked_store is None:
            raise AssetError('Fingerprint-pinned linking requires a retained library root: ' + key)
        if spec.source == 'bundled':
            path = (self.bundled_root / spec.path).resolve()
            if not path.is_relative_to(self.bundled_root):
                raise AssetError('Bundled asset path escapes the asset directory')
        else:
            import bpy
            path = Path(bpy.path.abspath(spec.path)).resolve()
        if link and self.linked_store is not None:
            if not spec.sha256 or not spec.self_contained:
                raise AssetError('Managed linking needs a published SHA-256 and a '
                                 'certified self-contained library: ' + key)
            try:
                # Edition and bundled-path rules still apply even when the
                # original distribution is absent after an addon update.
                retained = self.linked_store.existing(spec.sha256)
                if retained is not None:
                    return spec, retained
                if path.suffix.lower() != '.blend' or not path.is_file():
                    raise AssetUnavailable(f'Blender library is unavailable: {path}; '
                                           f'restore revision {spec.revision} ({spec.sha256})')
                return spec, self.linked_store.install(path, spec.sha256)
            except AssetStoreError as exc:
                raise AssetError(str(exc)) from exc
        if path.suffix.lower() != '.blend' or not path.is_file():
            raise AssetUnavailable(f'Blender library is unavailable: {path}')
        if spec.sha256 is not None and file_digest(path) != spec.sha256:
            raise AssetError(f'Library fingerprint mismatch: {path}; expected {spec.sha256}')
        return spec, path

    def load(self, reference, *, link=None):
        return self.load_many((reference,), link=link)[reference]

    def contains(self, reference, *, link=None):
        """Inspect a library catalog without allocating any datablocks.

        Edition, path and fingerprint errors remain errors. Only an absent
        named block returns False, allowing a feature to batch optional assets
        with its primary assets and avoid duplicate appended dependencies.
        """
        import bpy
        link = self.uses_link(reference.library, link)
        _spec, path = self.resolve(reference.library, link=link)
        collection = getattr(bpy.data, reference.kind, None)
        if not reference.name or collection is None or not hasattr(collection, 'get'):
            raise AssetError('Invalid Blender datablock reference: '+str(reference))
        with bpy.data.libraries.load(str(path), link=link) as (available, _requested):
            names = getattr(available, reference.kind, None)
            if not isinstance(names, list):
                raise AssetError('Invalid library datablock kind: '+reference.kind)
            return reference.name in names

    def load_many(self, references, *, link=None):
        import bpy
        references = tuple(dict.fromkeys(references))
        result, pending, resolved = {}, {}, {}
        # Resolve every library/edition before mutating Blender data.
        for ref in references:
            use_link = self.uses_link(ref.library, link)
            if ref.library not in resolved:
                resolved[ref.library] = self.resolve(ref.library, link=use_link)
            spec, path = resolved[ref.library]
            identity = (str(path), spec.revision, ref.kind, ref.name, use_link)
            # Keep legacy provenance keys readable. Pinned libraries also need
            # byte identity when a publisher mistakenly reuses a revision label.
            key = json.dumps(identity + (spec.sha256,) if spec.sha256 else identity)
            collection = getattr(bpy.data, ref.kind, None)
            if not ref.name or collection is None or not hasattr(collection, 'get'):
                raise AssetError('Invalid Blender datablock reference: ' + str(ref))
            value = self._cache.get(key)
            if value is not None:
                try:
                    if collection.get(value.name) != value:
                        value = None
                except ReferenceError:
                    value = None
            if value is None:
                for candidate in collection:
                    if not use_link and not candidate.library and candidate.get(self.provenance_key) == key:
                        value = candidate
                        break
                    if use_link and candidate.library and candidate.name == ref.name:
                        source = Path(bpy.path.abspath(candidate.library.filepath)).resolve()
                        if source == path:
                            value = candidate
                            break
            if value is not None:
                result[ref] = self._cache[key] = value
            else:
                pending.setdefault((path, use_link), []).append((ref, key))
        # Check all requested names first: missing assets never load a partial batch.
        for (path, use_link), entries in pending.items():
            with bpy.data.libraries.load(str(path), link=use_link) as (available, _requested):
                for ref, _key in entries:
                    names = getattr(available, ref.kind, None)
                    if not isinstance(names, list) or ref.name not in names:
                        raise AssetUnavailable(f'{path}: no {ref.kind} named {ref.name!r}')
        for (path, use_link), entries in pending.items():
            kinds = {}
            for ref, key in entries:
                kinds.setdefault(ref.kind, []).append((ref, key))
            with bpy.data.libraries.load(str(path), link=use_link) as (_available, requested):
                for kind, values in kinds.items():
                    setattr(requested, kind, [ref.name for ref, _key in values])
            for kind, values in kinds.items():
                for (ref, key), value in zip(values, getattr(requested, kind)):
                    if value is None:
                        raise AssetError(f'Blender could not load {ref} from {path}')
                    if not use_link:
                        value[self.provenance_key] = key
                    result[ref] = self._cache[key] = value
        return result
