#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys

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
from turbogenius.pyturbo.utils.env import pyturbo_data_dir
from turbogenius.pyturbo.namelist import Namelist

file_name = os.path.join(pyturbo_data_dir, "makefort10", "makefort10.input")
namelist = Namelist.parse_namelist_from_file(file_name)
logger.info(namelist.parameters)
logger.info(namelist.get_parameter("posunits"))
namelist.set_parameter("posunits", "ang")  # set here
logger.info(namelist.get_parameter("posunits"))
