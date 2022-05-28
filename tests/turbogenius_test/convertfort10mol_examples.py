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
from turbogenius.convertfort10mol_genius import Convertfort10mol_genius

convertfort10mol_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "convertfort10mol")
os.chdir(convertfort10mol_test_dir)

convertfort10mol_genius = Convertfort10mol_genius()
convertfort10mol_genius.generate_input()
convertfort10mol_genius.run()