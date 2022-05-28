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

trexio_to_turborvb_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "trexio_to_turborvb")

os.chdir(trexio_to_turborvb_test_dir)

# trexio to turborvb, check, as many cases as possible, using a workflow.