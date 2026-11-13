"""Numerical backend assembly; feature policy is injected through providers."""
import time
from dataclasses import replace
from .intersections import IntersectionProvider
from .patterns import PatternLibrary
from .providers import PathProvider, CarriagewayProvider, FurnitureProvider
from .updates import UpdateCoordinator, drain_steps, _PreparationProgress
from .validation import validate_overlaps_steps
from .profiles import ProfileProvider
from .bundles import BundleProvider

DEFAULT_FEATURES = ('paths', 'profiles', 'intersections', 'carriageways', 'bundles', 'furniture')
FEATURES = {'paths': PathProvider, 'intersections': IntersectionProvider,
            'carriageways': CarriagewayProvider, 'furniture': FurnitureProvider, 'profiles': ProfileProvider,
            'bundles': BundleProvider}


def register_feature(name, factory):
    """Register trusted addon extension code; saved files contain only IDs.

    Reopening never imports or executes Python module paths from scene data.
    An addon extension must register its provider before its scenes are updated.
    """
    if not isinstance(name, str) or not name or name in FEATURES or not callable(factory):
        raise ValueError('Feature provider needs a unique stable name and factory')
    FEATURES[name] = factory


class PythonMeshBackend:
    name = 'python'

    def __init__(self, settings, providers=None, *, feature_names=None):
        self.patterns = PatternLibrary(settings.pattern_cache_bytes)
        from .terrain.cache import TerrainConstraintCache
        self.terrain_constraints = TerrainConstraintCache()
        self.feature_names = tuple(DEFAULT_FEATURES if feature_names is None else feature_names)
        if providers is not None and feature_names is not None:
            raise ValueError('Choose explicit provider objects or registered feature names')
        if providers is None:
            unknown = set(self.feature_names) - FEATURES.keys()
            if unknown:
                raise ValueError('Required feature providers are not registered: ' + ', '.join(sorted(unknown)))
            providers = tuple(FEATURES[name]() for name in self.feature_names)
        self.providers = tuple(providers)
        self.updates = UpdateCoordinator(self.providers, self.patterns)

    def prepare(self, network, *, force=False):
        return drain_steps(self.prepare_steps(network, force=force))

    def prepare_steps(self, network, *, force=False, on_progress=None):
        """Prepare numerical products cooperatively without accepting them.

        UV fallback yields between products, and overlap validation yields at
        its existing numerical batch boundaries. Individual builders and UV
        mappings remain indivisible operations; the phase callback identifies
        those operations before they begin. No Blender or
        scheduler dependency is needed by this numerical producer.
        """
        progress = _PreparationProgress(on_progress)
        yield from progress.checkpoint('Checking street feature support')
        supported = set().union(*(provider.source_kinds for provider in self.providers))
        missing = {key.kind for key in network.records()} - supported
        if missing:
            raise ValueError('No registered numerical provider supports: ' + ', '.join(sorted(missing)))
        self.patterns.set_budget(network.settings.pattern_cache_bytes)
        prepared = yield from progress.forward(self.updates.prepare_steps(
            network, force=force, on_progress=on_progress))
        # Feature providers can supply their own charts and extra UV sets.
        # Every otherwise uncharted Part gets a shared physical-plan fallback,
        # after its final geometry has been assembled and before placement.
        from .providers import Part
        from .uv import DEFAULT_UV_MAP, PlanarMapping, apply_mapping
        uv_seconds = 0.
        products = dict(prepared.products)
        yield from progress.checkpoint('Checking street UV maps')
        for index, (key, part) in enumerate(products.items(), 1):
            if isinstance(part, Part) and part.mesh is not None and DEFAULT_UV_MAP not in part.mesh.uv_layers:
                yield from progress.checkpoint(f'Mapping street UVs for {key.text()}')
                before = time.perf_counter()
                products[key] = Part(part.owner, apply_mapping(part.mesh, PlanarMapping()),
                                     dict(part.metadata, uv_mapping='planar_xy_metres'))
                uv_seconds += time.perf_counter() - before
            elif index % 32 == 0:
                yield from progress.checkpoint()
        prepared = replace(prepared, products=products)
        prepared.timings['uv_seconds'] = uv_seconds
        yield from progress.checkpoint('Validating street overlaps')
        meshes = {key: part.mesh for key, part in prepared.products.items()
                  if hasattr(part, 'mesh') and part.mesh is not None}
        prepared.timings['overlap_check'] = yield from progress.forward(validate_overlaps_steps(
            meshes, changed=prepared.rebuilt, cell_size=max(1., network.settings.spacing * 2.),
            on_progress=on_progress))
        yield from progress.checkpoint('Street numerical validation complete')
        prepared.timings['geometry_seconds'] = prepared.timings['total']
        prepared.timings['total'] = progress.elapsed()
        return prepared

    def commit(self, prepared):
        self.updates.commit(prepared)

    def manifest(self):
        return self.updates.manifest()

    def release_published_meshes(self):
        """Blender owns published geometry; avoid a second whole-scene cache.

        Keep the small source-derived path/layout artifacts and output metadata.
        A changed mesh task always builds a fresh payload before publication.
        Standalone numerical callers can retain buffers for comparison instead.
        """
        from .providers import Part
        for product in self.updates.products.values():
            if isinstance(product, Part) and product.mesh is not None:
                product.metadata['bounds'] = product.mesh.bounds
                product.mesh = None
            elif callable(getattr(product, 'release_geometry', None)):
                product.release_geometry()
