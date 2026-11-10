"""
blosm_networkx
========

blosm_networkx is a Python package for the creation, manipulation, and study of the
structure, dynamics, and functions of complex networks.

See https://blosm_networkx.org for complete documentation.
"""

__version__ = "3.3rc0.dev0"


# These are imported in order as listed
from .exception import *

from . import utils
from .utils.backends import _dispatch

from . import classes
from .classes import filters
from .classes import *

from . import convert
from .convert import *

from . import relabel
from .relabel import *

from . import generators
from .generators import *

from . import readwrite
from .readwrite import *

# Need to test with SciPy, when available
from . import algorithms
from .algorithms import *
