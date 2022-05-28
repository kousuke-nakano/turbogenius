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
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from turbogenius.prep_genius import DFT_genius

prep_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "prep")
os.chdir(prep_test_dir)

prep_genius = DFT_genius(grid_size=[0.30, 0.30, 0.30], lbox=[3.0, 3.0, 3.0])
prep_genius.generate_input()
prep_genius.run()