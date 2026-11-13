"""Load the scheduler under a package shell without executing Blender startup."""
import importlib
import sys
from types import ModuleType
from release_target import ADDON_ROOT

PACKAGE = 'tested_blosm_scheduler'
if PACKAGE not in sys.modules:
    shell = ModuleType(PACKAGE)
    shell.__path__ = [str(ADDON_ROOT)]
    sys.modules[PACKAGE] = shell

core = importlib.import_module(PACKAGE + '.progressive.core')
downloads = importlib.import_module(PACKAGE + '.downloads')
transport = importlib.import_module(PACKAGE + '.downloads.transport')
worker = importlib.import_module(PACKAGE + '.downloads.worker')
