"""Initialize everything"""


import os

from pyfast import utilities

from .executor import Executor

ROOT = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(ROOT, "..", "VERSION")) as version_file:
    VERSION = version_file.read().strip()

__version__ = VERSION
