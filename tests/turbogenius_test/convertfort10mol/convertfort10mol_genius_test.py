#!python
# -*- coding: utf-8 -*-
import os, sys

# pyturbo modules
from turbogenius.convertfort10mol_genius import Convertfort10mol_genius

root_dir = os.path.dirname(__file__)


# convertfort10mol_genius
def test_convertfort10mol_genius_basic():
    os.chdir(root_dir)

    convertfort10mol_genius = Convertfort10mol_genius(
        fort10="fort.10_in",
        add_random_mo=True,
        grid_size=0.10,
        additional_mo=0,
    )

    assert convertfort10mol_genius.convertfort10mol.get_parameter("nmol") == 1
    assert convertfort10mol_genius.convertfort10mol.get_parameter("epsdgm") < 0

    convertfort10mol_genius.generate_input()
    convertfort10mol_genius.run()
    flags = convertfort10mol_genius.check_results(output_names=["out_mol"])
    assert flags

    os.chdir(root_dir)


# convertfort10mol_genius
def test_convertfort10mol_genius_options():
    os.chdir(root_dir)

    convertfort10mol_genius = Convertfort10mol_genius(
        fort10="fort.10_in",
        add_random_mo=False,
        grid_size=0.20,
        additional_mo=2,
    )

    assert convertfort10mol_genius.convertfort10mol.get_parameter("nmol") == 3
    assert convertfort10mol_genius.convertfort10mol.get_parameter("epsdgm") > 0
    assert convertfort10mol_genius.convertfort10mol.get_parameter("ax") == 0.20
    assert convertfort10mol_genius.convertfort10mol.get_parameter("ay") == 0.20
    assert convertfort10mol_genius.convertfort10mol.get_parameter("az") == 0.20

    convertfort10mol_genius.generate_input()
    convertfort10mol_genius.run()
    flags = convertfort10mol_genius.check_results(output_names=["out_mol"])
    assert flags

    os.chdir(root_dir)
