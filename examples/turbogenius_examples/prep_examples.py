#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys
import shutil

#set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo')
logger.setLevel("DEBUG")
stream_handler = StreamHandler()
stream_handler.setLevel("DEBUG")
handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)

#pyturbo modules
from turbogenius.prep_genius import DFT_genius

from turbogenius.utils_workflows.env import turbo_genius_root

# convertfort10mol
prefix="prep"
example_root_dir=os.path.join(turbo_genius_root, "examples", "turbogenius_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

prep_genius = DFT_genius(grid_size=[0.30, 0.30, 0.30], lbox=[3.0, 3.0, 3.0])
prep_genius.generate_input()
prep_genius.run()