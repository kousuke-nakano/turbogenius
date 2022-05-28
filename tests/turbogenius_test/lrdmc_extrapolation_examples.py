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
from turbogenius.lrdmc_genius import LRDMC_ext_genius

lrdmc_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "lrdmc_extrapolation")

os.chdir(lrdmc_test_dir)

lrdmc_ext_genius = LRDMC_ext_genius(
    fort10="fort.10",
    lrdmcsteps=200,
    bin_block=2,
    warmupblocks=2,
    correcting_factor=2,
    num_walkers=2,
    alat_list=[-0.3, -0.4],
    etry=-1.00,
    polynominal_order=2,
    twist_average=False,
    force_calc_flag=False,
    nonlocalmoves="dlatm"
)

lrdmc_ext_genius.generate_input()
lrdmc_ext_genius.run()
lrdmc_ext_genius.do_extrapolation()
logger.info(lrdmc_ext_genius.energy)
logger.info(lrdmc_ext_genius.energy_error)