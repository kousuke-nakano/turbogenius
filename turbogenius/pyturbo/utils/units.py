#!python -u
# -*- coding: utf-8 -*-

"""

pyturbo: units used in pyturbo

"""

from __future__ import print_function

# python modules


# set logger
from logging import getLogger

logger = getLogger("pyturbo").getChild(__name__)

"""
Units used in Turbo-Genius

Lentgh: Bohr
Energy: Ha

are defined as 1.0
"""

# Length
Bohr = 1.0
Angstrom = 1.0 / 0.529177210903  # Bohr

# Energy
Ha = 1.0
Ry = 2.0  # Ha
