"""Load the numerical package without importing Blender's addon entry point."""
import importlib.util
from pathlib import Path
import sys

NAME = 'tested_python_streets'
if NAME not in sys.modules:
    root = Path(__file__).resolve().parents[1] / 'way/street_network/python_backend'
    spec = importlib.util.spec_from_file_location(NAME, root / '__init__.py',
                                                 submodule_search_locations=[str(root)])
    package = importlib.util.module_from_spec(spec)
    sys.modules[NAME] = package
    spec.loader.exec_module(package)
