#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil
import numpy as np
from turbogenius.vmc_opt_genius import VMCopt_genius

root_dir = os.path.dirname(__file__)

# WIP! not comprehensive


# jasopt: 1b,2b,3b, detopt: fix
def test_vmcopt_genius_jasopt_1b_2b_3b_basic():
    os.chdir(root_dir)
    prefix = "jasopt_1b_2b_3b"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    vmcopt_genius = VMCopt_genius(
        fort10="fort.10",
        vmcoptsteps=500,
        steps=10,
        bin_block=1,
        warmupblocks=0,
        num_walkers=20,
        optimizer="lr",
        learning_rate=0.35,
        regularization=0.001,
        opt_onebody=True,
        opt_twobody=True,
        opt_det_mat=False,
        opt_jas_mat=True,
        opt_det_basis_exp=False,
        opt_jas_basis_exp=False,
        opt_det_basis_coeff=False,
        opt_jas_basis_coeff=False,
        twist_average=False,
    )
    vmcopt_genius.generate_input()
    vmcopt_genius.run()
    energy_list, energy_error_list = vmcopt_genius.get_energy()
    np.testing.assert_almost_equal(energy_list[-1], -1.17260688110847, decimal=6)
    np.testing.assert_almost_equal(
        energy_error_list[-1], 3.236482791458927e-3, decimal=5
    )
    vmcopt_genius.average(optwarmupsteps=400)
    os.chdir(root_dir)


""" for the time being.
# jasopt: 1b, detopt: fix
def test_vmcopt_genius_jasopt_1b_basic():
    os.chdir(root_dir)
    prefix = "jasopt_1b"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    vmcopt_genius = VMCopt_genius(
        fort10="fort.10",
        vmcoptsteps=500,
        steps=10,
        bin_block=1,
        warmupblocks=0,
        num_walkers=20,
        optimizer="lr",
        learning_rate=0.35,
        regularization=0.001,
        opt_onebody=True,
        opt_twobody=False,
        opt_det_mat=False,
        opt_jas_mat=False,
        opt_det_basis_exp=False,
        opt_jas_basis_exp=False,
        opt_det_basis_coeff=False,
        opt_jas_basis_coeff=False,
        twist_average=False,
    )
    vmcopt_genius.generate_input()
    vmcopt_genius.run()
    energy_list, energy_error_list = vmcopt_genius.get_energy()
    np.testing.assert_almost_equal(energy_list[-1], -1.17061954590115, decimal=6)
    np.testing.assert_almost_equal(
        energy_error_list[-1], 3.019892562222970e-3, decimal=6
    )
    vmcopt_genius.average(optwarmupsteps=400)
    os.chdir(root_dir)
"""


# jasopt: 2b, detopt: fix
def test_vmcopt_genius_jasopt_2b_basic():
    os.chdir(root_dir)
    prefix = "jasopt_2b"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    vmcopt_genius = VMCopt_genius(
        fort10="fort.10",
        vmcoptsteps=500,
        steps=10,
        bin_block=1,
        warmupblocks=0,
        num_walkers=20,
        optimizer="lr",
        learning_rate=0.35,
        regularization=0.001,
        opt_onebody=False,
        opt_twobody=True,
        opt_det_mat=False,
        opt_jas_mat=False,
        opt_det_basis_exp=False,
        opt_jas_basis_exp=False,
        opt_det_basis_coeff=False,
        opt_jas_basis_coeff=False,
        twist_average=False,
    )
    vmcopt_genius.generate_input()
    vmcopt_genius.run()
    energy_list, energy_error_list = vmcopt_genius.get_energy()
    np.testing.assert_almost_equal(energy_list[-1], -1.17702418461567, decimal=6)
    np.testing.assert_almost_equal(
        energy_error_list[-1], 1.709183522590725e-3, decimal=6
    )
    vmcopt_genius.average(optwarmupsteps=400)
    os.chdir(root_dir)


# jasopt: fix, detopt: matrix
def test_vmcopt_genius_detopt_mat_basic():
    os.chdir(root_dir)
    prefix = "detopt_mat"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    vmcopt_genius = VMCopt_genius(
        fort10="fort.10",
        vmcoptsteps=500,
        steps=10,
        bin_block=1,
        warmupblocks=0,
        num_walkers=20,
        optimizer="lr",
        learning_rate=0.35,
        regularization=0.001,
        opt_onebody=False,
        opt_twobody=False,
        opt_det_mat=True,
        opt_jas_mat=False,
        opt_det_basis_exp=False,
        opt_jas_basis_exp=False,
        opt_det_basis_coeff=False,
        opt_jas_basis_coeff=False,
        twist_average=False,
    )
    vmcopt_genius.generate_input()
    vmcopt_genius.run()
    energy_list, energy_error_list = vmcopt_genius.get_energy()
    np.testing.assert_almost_equal(energy_list[-1], -1.17307645275826, decimal=5)
    np.testing.assert_almost_equal(
        energy_error_list[-1], 2.344883420011352e-3, decimal=6
    )
    vmcopt_genius.average(optwarmupsteps=400)
    os.chdir(root_dir)
