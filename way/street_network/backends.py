"""Street rendering lifecycle and lazy backend registration, without bpy.

All backends consume Blosm's shared prepared Streets/Sections/Bundles. Keep
geometry generation, assets and backend-specific state in their implementations.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Callable

DEFAULT_BACKEND = 'python'
# Saved RNA values must never acquire a new meaning when a backend is added.
RETIRED_ENUM_VALUES = frozenset((1, 2))
RETIRED_IDENTIFIERS = frozenset(('geometry_nodes', 'legacy'))


class StreetRenderBackend(ABC):
    """Lifecycle used by the importer, independently of the mesh technology."""

    @abstractmethod
    def prepare(self):
        """Check dependencies and initialize this import's backend resources."""

    @abstractmethod
    def render(self, manager, data):
        """Generate output from the shared prepared network; return its context."""

    @abstractmethod
    def finalize(self):
        """Finish publication of a successful import."""

    @abstractmethod
    def cleanup(self):
        """Release temporary backend resources."""


@dataclass(frozen=True)
class StreetBackendDefinition:
    key: str
    label: str
    description: str
    enum_value: int
    factory: Callable
    draw_settings: Callable | None = None


class StreetBackendRegistry:
    def __init__(self):
        self._definitions = {}
        self._enum_items = ()
        # Blender borrows strings returned by EnumProperty callbacks. Retain
        # earlier tuples too when a provider registers/unregisters at runtime.
        self._enum_history = []

    def __len__(self):
        return len(self._definitions)

    def __contains__(self, key):
        return key in self._definitions

    def register(self, definition):
        if not isinstance(definition, StreetBackendDefinition):
            raise TypeError('Expected a StreetBackendDefinition')
        if not definition.key or not definition.label or not callable(definition.factory):
            raise ValueError('A street backend needs a key, label and callable factory')
        if type(definition.enum_value) is not int or definition.enum_value < 0:
            raise ValueError('Street backend enum values must be nonnegative integers')
        if definition.enum_value in RETIRED_ENUM_VALUES:
            raise ValueError('Street backend enum values 1 and 2 are retired')
        if definition.draw_settings is not None and not callable(definition.draw_settings):
            raise ValueError('Backend settings drawing must be callable')
        if definition.key in self._definitions:
            raise ValueError('Street backend already registered: '+definition.key)
        if any(d.enum_value == definition.enum_value for d in self._definitions.values()):
            raise ValueError('Street backend enum value already registered')
        self._definitions[definition.key] = definition
        self._refresh_enum_items()

    def unregister(self, key):
        if key == DEFAULT_BACKEND:
            raise ValueError('The default Python street backend must remain registered')
        self.get(key)
        del self._definitions[key]
        self._refresh_enum_items()

    def get(self, key):
        try:
            return self._definitions[key]
        except KeyError:
            raise ValueError('Unknown street backend: '+str(key)) from None

    def create(self, key, app, *, design=None, features=None, asset_loader=None):
        backend = self.get(key).factory(app, design=design, features=features, asset_loader=asset_loader)
        if not isinstance(backend, StreetRenderBackend):
            raise TypeError('Street backend factory must return a StreetRenderBackend')
        return backend

    def enum_items(self):
        return self._enum_items

    def _refresh_enum_items(self):
        self._enum_items = tuple((d.key, d.label, d.description, d.enum_value)
                                 for d in self._definitions.values())
        self._enum_history.append(self._enum_items)


def _create_python_backend(app, **options):
    from .python_backend.renderer import PythonStreetRenderer
    return PythonStreetRenderer(app, **options)


registry = StreetBackendRegistry()
registry.register(StreetBackendDefinition(
    DEFAULT_BACKEND, 'Python (experimental)',
    'Direct mesh generation; edit centerlines and press Update', 0, _create_python_backend))


def migrate_saved_selection(settings):
    """Normalize retired selections without overriding a newly registered GN backend."""
    raw = settings.get('streetBackend')
    if isinstance(raw, int) and raw in RETIRED_ENUM_VALUES:
        settings.streetBackend = DEFAULT_BACKEND
    elif isinstance(raw, str):
        if raw in registry:
            settings.streetBackend = raw
        elif raw in RETIRED_IDENTIFIERS:
            settings.streetBackend = DEFAULT_BACKEND
