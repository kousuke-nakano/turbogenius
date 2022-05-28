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
from turbogenius.vmc_genius import VMC_genius

vmc_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "vmc")

os.chdir(vmc_test_dir)

vmc_genius=VMC_genius(
                     fort10="fort.10",
                     vmcsteps=600,
                     bin_block=10,
                     warmupblocks=10,
                     num_walkers=40,
                     twist_average=False,
                     force_calc_flag=True
                     )

vmc_genius.generate_input()
vmc_genius.run()
vmc_genius.compute_energy_and_forces()
logger.info(vmc_genius.energy)
logger.info(vmc_genius.energy_error)