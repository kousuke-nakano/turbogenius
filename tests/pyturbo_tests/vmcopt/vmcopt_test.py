#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil

import numpy as np
from turbogenius.pyturbo.vmcopt import VMCopt

root_dir = os.path.dirname(__file__)
os.chdir(root_dir)

shutil.copyfile("fort.10_ori", "fort.10")
namelist = VMCopt.read_default_namelist()
vmcopt = VMCopt(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=False,
)
vmcopt.set_parameter("ngen", 10)
vmcopt.set_parameter("nweight", 100)
vmcopt.set_parameter("iopt", 1)
vmcopt.sanity_check()

assert vmcopt.namelist.get_parameter(parameter="ngen") == 10
assert vmcopt.namelist.get_parameter(parameter="nweight") == 100
assert vmcopt.namelist.get_parameter(parameter="iopt") == 1

vmcopt.generate_input(input_name="datasmin0.input")
vmcopt.run(input_name="datasmin0.input", output_name="out_min0")

vmcopt.set_parameter("ngen", 30)
vmcopt.set_parameter("nbinr", 2)
vmcopt.set_parameter("nweight", 3)
vmcopt.set_parameter("iopt", 0)
vmcopt.sanity_check()

assert vmcopt.namelist.get_parameter(parameter="ngen") == 30
assert vmcopt.namelist.get_parameter(parameter="nbinr") == 2
assert vmcopt.namelist.get_parameter(parameter="nweight") == 3
assert vmcopt.namelist.get_parameter(parameter="iopt") == 0

vmcopt.generate_input(input_name="datasmin1.input")
vmcopt.run(input_name="datasmin1.input", output_name="out_min1")

flags = vmcopt.check_results(output_names=["out_min0", "out_min1"])
assert flags

vmcopt.average_optimized_parameters(
    input_file_used="datasmin1.input", equil_steps=3, graph_plot=True
)

energy, error = vmcopt.get_energy(output_names=["out_min0", "out_min1"])

np.testing.assert_array_almost_equal(
    [energy[-1], error[-1]],
    [-1.16828602751887, 0.007059206136597875],
    decimal=12,
)
