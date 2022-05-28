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
from turbogenius.pyturbo.lrdmcopt import LRDMCopt

# LRDMC opt
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "lrdmcopt"))
namelist = LRDMCopt.read_default_namelist()
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
