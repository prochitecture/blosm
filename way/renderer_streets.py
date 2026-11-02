"""Stable importer facade for registered street rendering backends."""
from .street_network.backends import DEFAULT_BACKEND, registry


class StreetRenderer:
    """Delegate an import lifecycle without coupling setup to a mesh technology."""

    def __init__(self, app, *, backend=None, design=None, features=None, asset_loader=None):
        self.backend = backend if backend is not None else getattr(app, 'streetBackend', DEFAULT_BACKEND)
        self.backend = self.backend or DEFAULT_BACKEND
        self.implementation = registry.create(self.backend, app, design=design,
            features=features, asset_loader=asset_loader)

    def prepare(self):
        return self.implementation.prepare()

    def render(self, manager, data):
        return self.implementation.render(manager, data)

    def finalize(self):
        return self.implementation.finalize()

    def cleanup(self):
        return self.implementation.cleanup()
