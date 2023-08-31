#!python
# -*- coding: utf-8 -*-
import os, sys

import numpy as np

# pyturbo modules
from turbogenius.pyturbo.namelist import Namelist
from turbogenius.pyturbo.vmc import VMC

# VMC
input_parameters = {
    "&simulation": {
        "itestr4": 2,
        "ngen": 5000,
        "iopt": 1,
    },
    "&pseudo": {},
    "&vmc": {},
    "&readio": {},
    "&parameters": {},
}
fortran_namelist = Namelist(namelist=input_parameters)

# VMC run (cont.=0)
vmc = VMC(namelist=fortran_namelist)
vmc.generate_input(input_name="datasvmc0.input")
vmc.run(input_name="datasvmc0.input", output_name="out_vmc0")

# VMC run (cont.=1)
vmc.set_parameter("iopt", 0)
vmc.generate_input(input_name="datasvmc1.input")
vmc.run(input_name="datasvmc1.input", output_name="out_vmc1")

# VMC, check if the jobs finished correctly
flags = vmc.check_results(output_names=["out_vmc0", "out_vmc1"])

# Reblock the MCMC samples
energy, error = vmc.get_energy(init=10, bin=10)
# print(f"VMC energy = {energy:.4f} +- {error:.4f}")
np.testing.assert_array_almost_equal(
    [energy, error], [-1.17507593466782, 6.158026435079502e-4], decimal=12
)
