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
from turbogenius.pyturbo.vmc import VMC
from turbogenius.utils_workflows.env import turbo_genius_root

#VMC
prefix="vmc"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

namelist=VMC.read_default_namelist()
vmc=VMC(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=False,
)
vmc.set_parameter("ngen", 100)
vmc.set_parameter("iopt", 1)
vmc.set_parameter("ieskin", 1, "&parameters")
vmc.sanity_check()
vmc.generate_input(input_name="datasvmc0.input")
vmc.run(input_name="datasvmc0.input", output_name="out_vmc0")

vmc.set_parameter("ngen", 200)
vmc.set_parameter("iopt", 0)
vmc.sanity_check()
vmc.generate_input(input_name="datasvmc1.input")
vmc.run(input_name="datasvmc1.input", output_name="out_vmc1")

flags=vmc.check_results(output_names=["out_vmc0", "out_vmc1"])
logger.info(f"check_results={flags}")

energy,energy_error=vmc.get_energy(init=10, bin=10)
force_matrix,force_error_matrix=vmc.get_forces(init=10, bin=10)
logger.info(energy)
logger.info(energy_error)
logger.info(force_matrix)
logger.info(force_error_matrix)
