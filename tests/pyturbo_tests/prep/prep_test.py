#!python
# -*- coding: utf-8 -*-
import os, sys
import numpy as np
from turbogenius.pyturbo.prep import Prep

root_dir = os.path.dirname(__file__)


def test_prep_basic():
    os.chdir(root_dir)

    namelist = Prep.read_default_namelist()

    # open boundary condition
    prep = Prep(
        in_fort10="fort.10",
        namelist=namelist,
    )

    np.array_equal(prep.nelocc_list, [2, 2])
    np.array_equal(prep.neloccdn_list, [])

    prep.generate_input(input_name="prep.input")
    prep.run(input_name="prep.input", output_name="out_prep")
    flags = prep.check_results(output_names=["out_prep"])
    assert flags

    # twist average case
    kpoints = [
        [
            [0.0, 0.0, 0.0, 1.0],
            [0.2, 0.2, 0.2, 0.5],
            [0.3, 0.3, 0.3, 0.4],
        ],
        [
            [-0.0, 0.0, 0.0, 1.0],
            [-0.2, 0.2, 0.2, 0.5],
            [-0.3, 0.3, 0.3, 0.4],
        ],
    ]

    prep = Prep(
        in_fort10="fort.10",
        namelist=namelist,
        twist_average=2,
    )
    assert prep.namelist.get_parameter(parameter="yes_kpoints") == ".false."
    prep.manual_kpoints = kpoints
    assert prep.namelist.get_parameter(parameter="yes_kpoints") == ".true."
    assert prep.namelist.get_parameter(parameter="yeswrite10") == ".true."
    assert prep.namelist.get_parameter(parameter="kp_type") == 2
    assert prep.namelist.get_parameter(parameter="nk1") == 3

    os.chdir(root_dir)
