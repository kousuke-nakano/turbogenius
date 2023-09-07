#!python
# -*- coding: utf-8 -*-
import os, sys

# pyturbo modules
from turbogenius.pyturbo.convertfort10mol import Convertfort10mol
from turbogenius.utils_workflows.env import turbo_genius_root

# convertfort10mol
root_dir = os.path.dirname(__file__)


def test_convertfort10mol_basic():
    os.chdir(root_dir)

    convertfort10mol = Convertfort10mol.parse_from_default_namelist(
        in_fort10="fort.10_in"
    )
    convertfort10mol.generate_input(input_name="convertfort10mol.input")
    convertfort10mol.run(input_name="convertfort10mol.input", output_name="out_mol")
    flags = convertfort10mol.check_results(output_names=["out_mol"])
    assert flags

    os.chdir(root_dir)
