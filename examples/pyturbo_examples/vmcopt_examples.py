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
from turbogenius.pyturbo.vmcopt import VMCopt
from turbogenius.utils_workflows.env import turbo_genius_root

#VMCopt
prefix="vmcopt"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

namelist=VMCopt.read_default_namelist()
vmcopt=VMCopt(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=False,
)
vmcopt.set_parameter("ngen", 10)
vmcopt.set_parameter("nweight", 100)
vmcopt.set_parameter("iopt", 1)
#vmcopt.set_parameter("ieskin", 1, "&parameters")
vmcopt.sanity_check()
vmcopt.generate_input(input_name="datasmin0.input")
vmcopt.run(input_name="datasmin0.input", output_name="out_min0")

vmcopt.set_parameter("ngen", 30)
vmcopt.set_parameter("nbinr", 2)
vmcopt.set_parameter("nweight", 3)
vmcopt.set_parameter("iopt", 0)
vmcopt.sanity_check()
vmcopt.generate_input(input_name="datasmin1.input")
vmcopt.run(input_name="datasmin1.input", output_name="out_min1")

flags=vmcopt.check_results(output_names=["out_min0", "out_min1"])
logger.info(f"check_results={flags}")

vmcopt.average_optimized_parameters(input_file_used="datasmin1.input",  equil_steps=3, graph_plot=True)

energy, error=vmcopt.get_energy(output_names=["out_min0", "out_min1"])