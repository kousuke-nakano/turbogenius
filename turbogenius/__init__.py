import os
import sys
import glob

__all__ = [
    os.path.split(os.path.splitext(file)[0])[1]
    for file in glob.glob(
        os.path.join(os.path.dirname(__file__), "[a-zA-Z0-9]*.py")
    )
]
from . import *
from .pyturbo import *

# The following line is a workaround for the compatibility between
# the old and new turbogenius versions! (i.e. name spaces in import)
# It will be removed in the next release.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
