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
from turbogenius.vmc_genius import VMC_genius
from turbogenius.utils_workflows.env import turbo_genius_root

# vmc
prefix="vmc"
example_root_dir=os.path.join(turbo_genius_root, "examples", "turbogenius_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

vmc_genius=VMC_genius(
                     fort10="fort.10",
                     vmcsteps=600,
                     num_walkers=40,
                     twist_average=False,
                     force_calc_flag=True
                     )

vmc_genius.generate_input()
vmc_genius.run()
vmc_genius.compute_energy_and_forces(bin_block=10, warmupblocks=10)
logger.info(vmc_genius.energy)
logger.info(vmc_genius.energy_error)