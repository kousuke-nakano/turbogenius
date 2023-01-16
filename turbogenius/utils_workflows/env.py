#!python -u
# -*- coding: utf-8 -*-

from __future__ import print_function

# python modules
import os

# set logger
from logging import getLogger

logger = getLogger("turbogenius").getChild(__name__)

# turbo-genius related path lists
turbo_genius_root = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../")
)
turbo_genius_source_dir = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
)
turbo_genius_data_dir = os.path.join(turbo_genius_source_dir, "database")
turbo_genius_tmp_dir = os.path.join(
    os.path.abspath(os.environ["HOME"]), ".turbo_genius_tmp"
)

# generate pyturbo temp. dir.
os.makedirs(turbo_genius_tmp_dir, exist_ok=True)
