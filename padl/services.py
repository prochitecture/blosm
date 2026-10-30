"""Import shared addon services without registering the addon in a build process."""
import importlib
from pathlib import Path
import sys
from types import ModuleType

from .errors import PadlError

_root = None
_package = None


def configure_blosm(root=None, *, package=None):
    """Opt into Blosm services by checkout path or installed addon package.

    Parsing and source-based GN builds do not require this integration. Asset
    loading, modifier access and the Python hipped kernel currently reuse Blosm.
    """
    global _root, _package
    if (root is None) == (package is None):
        raise ValueError('Specify either a Blosm checkout root or its installed package')
    if package is not None:
        module = importlib.import_module(package)
        candidate = Path(module.__file__).resolve().parent
    else:
        candidate = Path(root).resolve()
    if not (candidate / '__init__.py').is_file() or not (candidate / 'lib').is_dir():
        raise PadlError('Configured directory does not contain Blosm integration services', code='missing_integration')
    if _root is not None and (_root != candidate or _package != package):
        raise PadlError('A different Blosm integration is already configured in this process', code='integration_conflict')
    _root, _package = candidate, package


def blosm_root():
    if _root is None:
        raise PadlError('Call padl.services.configure_blosm with the Blosm checkout or installed package first',
                        code='missing_integration')
    return _root


def addon_service(relative_name):
    root = blosm_root()
    if _package:
        return importlib.import_module(_package + '.' + relative_name)
    # Load shared services without running addon registration in development.
    prefix = '_blosm_padl_services'
    parent = relative_name.split('.')[0]
    for name, path in ((prefix, root), (prefix + '.' + parent, root / parent)):
        if name not in sys.modules:
            module = ModuleType(name)
            module.__path__ = [str(path)]
            module.__package__ = name
            sys.modules[name] = module
    return importlib.import_module(prefix + '.' + relative_name)
