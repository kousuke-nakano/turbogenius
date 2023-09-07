#!python
# -*- coding: utf-8 -*-
import os, sys

import numpy as np

# turbogenius modules
from turbogenius.vmc_genius import VMC_genius

root_dir = os.path.dirname(__file__)


# VMC
def test_vmc_genius_basic():
    # create a vmc_genius instance
    vmc_genius = VMC_genius(
        fort10="fort.10",  # WF file of the H2 dimer
        vmcsteps=600,  # The number of MCMC steps
        num_walkers=40,  # The number of walkers
    )

    # generate an input file
    vmc_genius.generate_input(input_name="datas_vmc.input")

    # launch a VMC run
    vmc_genius.run(input_name="datas_vmc.input", output_name="out_vmc")

    # compute energy and forces
    vmc_genius.compute_energy_and_forces(bin_block=10, warmupblocks=10)

    # print vmc energy
    # print(f"VMC energy = {vmc_genius.energy:.5f} +- {vmc_genius.energy_error:.5f} Ha")
    np.testing.assert_almost_equal(vmc_genius.energy, -1.17426860723634, decimal=10)
    np.testing.assert_almost_equal(
        vmc_genius.energy_error, 2.827731334376716e-4, decimal=10
    )
