#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil

import numpy as np

# pyturbo modules
from turbogenius.pyturbo.lrdmcopt import LRDMCopt

shutil.copyfile("fort.10_ori", "fort.10")
namelist = LRDMCopt.read_default_namelist()
namelist.set_parameter(parameter="ngen", value=500)
namelist.set_parameter(parameter="nweight", value=100)
namelist.set_parameter(parameter="alat", value=0.30)
namelist.set_parameter(parameter="tpar", value=0.00001)

lrdmcopt = LRDMCopt(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=False,
)
lrdmcopt.sanity_check()

assert lrdmcopt.namelist.get_parameter(parameter="ngen") == 500
assert lrdmcopt.namelist.get_parameter(parameter="nweight") == 100
assert lrdmcopt.namelist.get_parameter(parameter="alat") == 0.30
assert lrdmcopt.namelist.get_parameter(parameter="tpar") == 0.00001
assert lrdmcopt.namelist.get_parameter(parameter="iopt") == 1

lrdmcopt.generate_input(input_name="datasfn_opt.input")
lrdmcopt.run(input_name="datasfn_opt.input", output_name="out_fn_opt")

flags = lrdmcopt.check_results(output_names=["out_fn_opt"])
assert flags

lrdmcopt.average_optimized_parameters(
    input_file_used="datasfn_opt.input", equil_steps=2, graph_plot=False
)

energy, error = lrdmcopt.get_energy(output_names=["out_fn_opt"])

np.testing.assert_array_almost_equal(
    [energy[-1], error[-1]],
    [-8.08925591853051, 0.05099428498660021],
    decimal=12,
)
