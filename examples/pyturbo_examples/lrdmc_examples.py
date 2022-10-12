#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys, shutil

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
from turbogenius.pyturbo.lrdmc import LRDMC
from turbogenius.utils_workflows.env import turbo_genius_root

# LRDMC
prefix="lrdmc"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

namelist = LRDMC.read_default_namelist()
namelist.set_parameter(parameter='ngen', value=500)
namelist.set_parameter(parameter='alat', value=0.5)
lrdmc = LRDMC(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=False,
)
lrdmc.sanity_check()
lrdmc.set_parameter(parameter="ieskin", value=1, namelist="&parameters")  # edit parameter
lrdmc.generate_input(input_name="datasfn.input")
lrdmc.run(input_name="datasfn.input", output_name="out_fn")

flags = lrdmc.check_results(output_names=["out_fn"])
logger.info(f"check_results={flags}")

energy, error = lrdmc.get_energy(init=10, correct=3, bin=5)
logger.info(f"{energy}+-{error} Ha")
force_matrix, force_error_matrix = lrdmc.get_forces(init=10, correct=3, bin=10)
logger.info(force_matrix)
logger.info(force_error_matrix)
