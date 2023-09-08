#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil

# pyturbo modules
from turbogenius.prep_genius import DFT_genius

root_dir = os.path.dirname(__file__)


# prep, WIP! The test is not comprehensive
def test_prep_genius_basic():
    os.chdir(root_dir)
    prep_genius = DFT_genius(grid_size=[0.30, 0.30, 0.30], lbox=[3.0, 3.0, 3.0])
    prep_genius.generate_input()
    prep_genius.run()
    flags = prep_genius.check_results()
    assert flags
