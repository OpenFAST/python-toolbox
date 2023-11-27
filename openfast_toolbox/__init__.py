"""Initialize everything"""
import os
from openfast_toolbox.common import *
# Make main io tools available
from .io import read
from .io.fast_input_file import FASTInputFile
from .io.fast_output_file import FASTOutputFile
from .io.fast_input_deck import FASTInputDeck

# Add version to package
with open(os.path.join(os.path.dirname(__file__), "..", "VERSION")) as fid:
    __version__ = fid.read().strip()


