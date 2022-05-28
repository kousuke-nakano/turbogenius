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
from turbogenius.lrdmc_genius import LRDMC_genius

lrdmc_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "lrdmc")

os.chdir(lrdmc_test_dir)

lrdmc_genius = LRDMC_genius(
    fort10="fort.10",
    lrdmcsteps=500,
    bin_block=10,
    warmupblocks=10,
    correcting_factor=2,
    num_walkers=20,
    alat=-0.4,
    etry=-1.00,
    twist_average=False,
    force_calc_flag=False,
    nonlocalmoves="dlatm" # tmove, dla, dlatm
)

lrdmc_genius.run()
lrdmc_genius.compute_energy_and_forces()
logger.info(lrdmc_genius.energy)
logger.info(lrdmc_genius.energy_error)