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
from turbogenius.pyturbo.lrdmcopt import LRDMCopt
from turbogenius.utils_workflows.env import turbo_genius_root

# LRDMC opt
prefix="lrdmcopt"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

namelist = LRDMCopt.read_default_namelist()
namelist.set_parameter(parameter='ngen', value=500)
namelist.set_parameter(parameter='nweight', value=10)
namelist.set_parameter(parameter='alat', value=0.5)

lrdmcopt = LRDMCopt(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=False,
)
lrdmcopt.sanity_check()
lrdmcopt.generate_input(input_name="datasfn_opt.input")
lrdmcopt.run(input_name="datasfn_opt.input", output_name="out_fn_opt")

flags = lrdmcopt.check_results(output_names=["out_fn_opt"])
logger.info(f"check_results={flags}")

lrdmcopt.average_optimized_parameters(input_file_used="datasfn_opt.input", equil_steps=2, graph_plot=False)
