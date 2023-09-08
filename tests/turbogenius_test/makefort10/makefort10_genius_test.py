#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil

from turbogenius.makefort10_genius import Makefort10_genius

root_dir = os.path.dirname(__file__)

# to be improved!! hese tests are not comprehensive!!


# Hydrogen dimer, ccECP
def test_makefort10_hydrogen_dimer_ccECP():
    os.chdir(root_dir)
    str_input = "H2.xyz"
    prefix, _ = os.path.splitext(str_input)
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, str_input),
        os.path.join(root_dir, prefix, str_input),
    )
    os.chdir(os.path.join(root_dir, prefix))

    makefort10_genius = Makefort10_genius(
        structure_file=str_input,
        det_basis_set="cc-pVTZ",
        jas_basis_set="cc-pVTZ",
        det_contracted_flag=True,
        jas_contracted_flag=True,
        all_electron_jas_basis_set=True,
        pseudo_potential="ccECP",
        det_cut_basis_option=True,
        jas_cut_basis_option=True,
        complex=False,
    )
    makefort10_genius.generate_input()
    makefort10_genius.run()
    flags = makefort10_genius.check_results()
    assert flags
    os.chdir(root_dir)


# Benzene, all-electrons
def test_makefort10_benzene():
    os.chdir(root_dir)
    str_input = "benzene.xyz"
    prefix, _ = os.path.splitext(str_input)
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, str_input),
        os.path.join(root_dir, prefix, str_input),
    )
    os.chdir(os.path.join(root_dir, prefix))

    makefort10_genius = Makefort10_genius(
        structure_file=str_input,
        det_basis_set="cc-pVQZ",
        jas_basis_set="cc-pVDZ",
        det_contracted_flag=True,
        jas_contracted_flag=True,
        all_electron_jas_basis_set=True,
        pseudo_potential=None,
        det_cut_basis_option=True,
        jas_cut_basis_option=True,
        complex=False,
    )
    makefort10_genius.generate_input()
    makefort10_genius.run()
    flags = makefort10_genius.check_results()
    assert flags
    os.chdir(root_dir)


# Diamond, PP, orthorhombic case, BFD
def test_makefort10_diamond():
    os.chdir(root_dir)
    str_input = "Diamond.cif"
    prefix, _ = os.path.splitext(str_input)
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, str_input),
        os.path.join(root_dir, prefix, str_input),
    )
    os.chdir(os.path.join(root_dir, prefix))

    makefort10_genius = Makefort10_genius(
        structure_file=str_input,
        det_basis_set="cc-pVQZ",
        jas_basis_set="cc-pVDZ",
        det_contracted_flag=True,
        jas_contracted_flag=True,
        all_electron_jas_basis_set=True,
        pseudo_potential=None,
        det_cut_basis_option=True,
        jas_cut_basis_option=True,
        complex=False,
    )
    makefort10_genius.generate_input()
    makefort10_genius.run()
    flags = makefort10_genius.check_results()
    assert flags
    os.chdir(root_dir)
