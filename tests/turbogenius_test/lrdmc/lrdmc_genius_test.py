#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil
import numpy as np
from turbogenius.lrdmc_genius import LRDMC_genius

root_dir = os.path.dirname(__file__)


# lrdmc, all-electrons
def test_lrdmc_genius_all_electrons_basic():
    os.chdir(root_dir)
    prefix = "ae"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10_ae"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    lrdmc_genius = LRDMC_genius(
        fort10="fort.10",
        lrdmcsteps=100,
        num_walkers=5,
        alat=-0.4,
        etry=-1.00,
        twist_average=False,
        force_calc_flag=False,
    )

    lrdmc_genius.generate_input()
    lrdmc_genius.run()
    lrdmc_genius.compute_energy_and_forces(
        bin_block=5, warmupblocks=5, correcting_factor=2
    )

    np.testing.assert_almost_equal(lrdmc_genius.energy, -1.17825782756524, decimal=10)
    np.testing.assert_almost_equal(
        lrdmc_genius.energy_error, 1.587798623416236e-3, decimal=10
    )

    os.chdir(root_dir)


# lrdmc, pp with dlatm
def test_lrdmc_genius_dlatm_basic():
    os.chdir(root_dir)
    prefix = "pp-dlatm"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10_pp"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    shutil.copy(
        os.path.join(root_dir, "pseudo.dat"),
        os.path.join(root_dir, prefix, "pseudo.dat"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    lrdmc_genius = LRDMC_genius(
        fort10="fort.10",
        lrdmcsteps=100,
        num_walkers=5,
        alat=-0.4,
        etry=-1.00,
        twist_average=False,
        force_calc_flag=False,
        nonlocalmoves="dlatm",  # tmove, dla, dlatm
    )

    lrdmc_genius.generate_input()
    lrdmc_genius.run()
    lrdmc_genius.compute_energy_and_forces(
        bin_block=5, warmupblocks=5, correcting_factor=2
    )

    np.testing.assert_almost_equal(lrdmc_genius.energy, -8.02726816434455, decimal=10)
    np.testing.assert_almost_equal(
        lrdmc_genius.energy_error, 2.278447489670330e-2, decimal=10
    )

    os.chdir(root_dir)


# lrdmc, pp with dla
def test_lrdmc_genius_dla_basic():
    os.chdir(root_dir)
    prefix = "pp-dla"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10_pp"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    shutil.copy(
        os.path.join(root_dir, "pseudo.dat"),
        os.path.join(root_dir, prefix, "pseudo.dat"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    lrdmc_genius = LRDMC_genius(
        fort10="fort.10",
        lrdmcsteps=100,
        num_walkers=5,
        alat=-0.4,
        etry=-1.00,
        twist_average=False,
        force_calc_flag=False,
        nonlocalmoves="dla",  # tmove, dla, dlatm
    )

    lrdmc_genius.generate_input()
    lrdmc_genius.run()
    lrdmc_genius.compute_energy_and_forces(
        bin_block=5, warmupblocks=5, correcting_factor=2
    )

    np.testing.assert_almost_equal(lrdmc_genius.energy, -8.10605306381933, decimal=10)
    np.testing.assert_almost_equal(
        lrdmc_genius.energy_error, 1.889584491014384e-2, decimal=10
    )

    os.chdir(root_dir)


def test_lrdmc_genius_tmove_basic():
    os.chdir(root_dir)
    prefix = "pp-tmove"
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, "fort.10_pp"),
        os.path.join(root_dir, prefix, "fort.10"),
    )
    shutil.copy(
        os.path.join(root_dir, "pseudo.dat"),
        os.path.join(root_dir, prefix, "pseudo.dat"),
    )
    os.chdir(os.path.join(root_dir, prefix))

    lrdmc_genius = LRDMC_genius(
        fort10="fort.10",
        lrdmcsteps=100,
        num_walkers=5,
        alat=-0.4,
        etry=-1.00,
        twist_average=False,
        force_calc_flag=False,
        nonlocalmoves="tmove",  # tmove, dla, dlatm
    )

    lrdmc_genius.generate_input()
    lrdmc_genius.run()
    lrdmc_genius.compute_energy_and_forces(
        bin_block=5, warmupblocks=5, correcting_factor=2
    )

    np.testing.assert_almost_equal(lrdmc_genius.energy, -8.02574492184647, decimal=10)
    np.testing.assert_almost_equal(
        lrdmc_genius.energy_error, 2.202829058518887e-2, decimal=10
    )

    os.chdir(root_dir)
