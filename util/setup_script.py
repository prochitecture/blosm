"""Load user-selected setup scripts within this extension's namespace."""
import importlib.util
from pathlib import Path
import sys

from .. import __package__ as base_package


def load_setup_script(filepath):
    # Setup scripts use the same relative imports as modules in blosm.setup.
    # A fixed private name avoids exposing filenames as global module names.
    name = base_package + '.setup._user_script'
    spec = importlib.util.spec_from_file_location(name, Path(filepath))
    module = importlib.util.module_from_spec(spec)
    previous = sys.modules.get(name)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        if previous is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = previous
        raise
    return module.setup
