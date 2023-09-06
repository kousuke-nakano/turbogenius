#!python
# -*- coding: utf-8 -*-
import os, sys
import numpy as np

# pyturbo modules
from turbogenius.pyturbo.lrdmc import LRDMC

root_dir = os.path.dirname(__file__)
os.chdir(root_dir)

namelist = LRDMC.read_default_namelist()
namelist.set_parameter(parameter="ngen", value=500)
namelist.set_parameter(parameter="alat", value=0.5)
lrdmc = LRDMC(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=False,
)
lrdmc.sanity_check()
lrdmc.set_parameter(
    parameter="ieskin", value=1, namelist="&parameters"
)  # edit parameter
lrdmc.generate_input(input_name="datasfn.input")
lrdmc.run(input_name="datasfn.input", output_name="out_fn")

flags = lrdmc.check_results(output_names=["out_fn"])
assert flags

energy, error = lrdmc.get_energy(init=10, correct=3, bin=5)
# print(f"VMC energy = {energy:.4f} +- {error:.4f}")
np.testing.assert_array_almost_equal(
    [energy, error], [-1.17234824773288, 2.457302575528722e-3], decimal=12
)
force_matrix, force_error_matrix = lrdmc.get_forces(init=10, correct=3, bin=5)
np.testing.assert_array_almost_equal(
    [force_matrix[0][2], force_error_matrix[0][2]],
    [-0.005737375445, 0.018993667436],
    decimal=12,
)
