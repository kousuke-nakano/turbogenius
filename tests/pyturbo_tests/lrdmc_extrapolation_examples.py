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
from turbogenius.pyturbo.lrdmc_extrapolation import LRDMC_extrapolation

# LRDMC
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "lrdmc_extrapolation"))
lrdmc_ext = LRDMC_extrapolation.parse_from_default_namelist(etry=-1.10, alat_list=[-0.40, -0.50])
# if you want to edit a value
# lrdmc_ext.set_parameters()
# if you want to see a value
# lrdmc_ext.get_parameters()
lrdmc_ext.generate_input()
lrdmc_ext.run()
lrdmc_ext.get_extrapolated_energy()