"""Initialize everything"""

# Import the local packages for external use
from .case_map import CASE_MAP
from .executor import Executor

import os
ROOT = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(ROOT, "..", "VERSION")) as version_file:
    VERSION = version_file.read().strip()
__version__ = VERSION
