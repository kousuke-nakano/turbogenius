#!python
# -*- coding: utf-8 -*-
import os, sys
import numpy as np

# pyturbo modules
from turbogenius.pyturbo.pseudopotentials import Pseudopotentials
from turbogenius.pyturbo.utils.downloader import BFD, ccECP

# basis sets
root_dir = os.path.dirname(__file__)

# BFD pseudo potential
bfd = BFD(pseudo_potential_output_dir=root_dir)
bfd.to_file(basis_list=["vtz"], element_list=["C", "Ca"])

# ccECP pseudo potential
ccECP = ccECP(pseudo_potential_output_dir=root_dir)
ccECP.to_file(basis_list=["cc-pVDZ"], element_list=["C", "Ca"])

# pseudo potential write
file_name = "C_ccECP.pseudo"
pseudopotentials = (
    Pseudopotentials.parse_pseudopotential_from_gamess_format_files(
        files=[file_name, file_name]
    )
)
pseudopotentials.set_cutoffs()
pseudopotentials.write_pseudopotential_turborvb_file()

# pseudo potential read
turbo_pseudopotentials = (
    Pseudopotentials.parse_pseudopotential_from_turborvb_pseudo_dat(
        file="pseudo.dat"
    )
)
assert np.array_equal(turbo_pseudopotentials.element_list, [None, None])
assert np.array_equal(turbo_pseudopotentials.max_ang_mom_plus_1, [1, 1])
assert np.array_equal(turbo_pseudopotentials.z_core, [None, None])
assert np.array_equal(turbo_pseudopotentials.cutoff, [1.41, 1.41])
assert np.array_equal(
    turbo_pseudopotentials.nucleus_index, [0, 0, 0, 0, 1, 1, 1, 1]
)
assert np.array_equal(turbo_pseudopotentials.ang_mom, [0, 1, 1, 1, 0, 1, 1, 1])
assert np.array_equal(
    turbo_pseudopotentials.exponent,
    [7.76079, 14.43502, 8.39889, 7.38188, 7.76079, 14.43502, 8.39889, 7.38188],
)
assert np.array_equal(
    turbo_pseudopotentials.coefficient,
    [52.13345, 4.0, 57.74008, -25.81955, 52.13345, 4.0, 57.74008, -25.81955],
)
assert np.array_equal(
    turbo_pseudopotentials.power, [0.0, -1.0, 1.0, 0.0, 0.0, -1.0, 1.0, 0.0]
)

os.chdir(root_dir)
