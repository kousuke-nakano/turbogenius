#!python
# -*- coding: utf-8 -*-
import os
from turbogenius.pyturbo.basis_set import Basis_set, Basis_sets

data_dir = os.path.dirname(os.path.abspath(__file__))


def test_basis_set():
    basis_set = Basis_set.parse_basis_set_info_from_gamess_format_file(
        os.path.join(data_dir, "C_cc-pVDZ.bas")
    )
    shell_ang_mom_test = basis_set.shell_ang_mom
    shell_index_test = basis_set.shell_index
    exponent_list_test = basis_set.exponent_list
    coefficient_list_test = basis_set.coefficient_list
    coefficient_imag_list_test = basis_set.coefficient_imag_list

    # reference values of C_cc-pVDZ.basis
    shell_ang_mom_ref = [0, 0, 0, 1, 1, 2]
    shell_index_ref = [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        3,
        3,
        3,
        3,
        4,
        5,
    ]
    exponent_list_ref = [
        6665.0,
        1000.0,
        228.0,
        64.71,
        21.06,
        7.495,
        2.797,
        0.5215,
        0.1596,
        6665.0,
        1000.0,
        228.0,
        64.71,
        21.06,
        7.495,
        2.797,
        0.5215,
        0.1596,
        0.1596,
        9.439,
        2.002,
        0.5456,
        0.1517,
        0.1517,
        0.55,
    ]
    coefficient_list_ref = [
        0.000692,
        0.005329,
        0.027077,
        0.101718,
        0.27474,
        0.448564,
        0.285074,
        0.015204,
        -0.003191,
        -0.000146,
        -0.001154,
        -0.005725,
        -0.023312,
        -0.063955,
        -0.149981,
        -0.127262,
        0.544529,
        0.580496,
        1.0,
        0.038109,
        0.20948,
        0.508557,
        0.468842,
        1.0,
        1.0,
    ]
    coefficient_imag_list_ref = []

    # assersion tests
    assert shell_ang_mom_test == shell_ang_mom_ref
    assert shell_index_test == shell_index_ref
    assert exponent_list_test == exponent_list_ref
    assert coefficient_list_test == coefficient_list_ref
    assert coefficient_imag_list_test == coefficient_imag_list_ref


def test_basis_sets():
    basis_sets = Basis_sets.parse_basis_sets_from_gamess_format_files(
        files=[
            os.path.join(data_dir, "C_cc-pVDZ.bas"),
            os.path.join(data_dir, "H_cc-pVTZ.bas"),
        ]
    )

    nuclei_num_test = basis_sets.nuclei_num
    nucleus_index_test = basis_sets.nucleus_index
    shell_ang_mom_test = basis_sets.shell_ang_mom
    shell_ang_mom_turbo_notation_test = basis_sets.shell_ang_mom_turbo_notation
    shell_index_test = basis_sets.shell_index
    coefficient_test = basis_sets.coefficient
    exponent_test = basis_sets.exponent

    nuclei_num_test_ref = 2
    nucleus_index_ref = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
    shell_ang_mom_ref = [0, 0, 0, 1, 1, 2, 0, 0, 0, 1, 1, 2]
    shell_ang_mom_turbo_notation_ref = [
        300,
        300,
        16,
        400,
        36,
        37,
        16,
        300,
        16,
        36,
        36,
        37,
    ]
    shell_index_ref = [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        3,
        3,
        3,
        3,
        4,
        5,
        6,
        7,
        7,
        7,
        7,
        7,
        8,
        9,
        10,
        11,
    ]
    coefficient_ref = [
        0.000692,
        0.005329,
        0.027077,
        0.101718,
        0.27474,
        0.448564,
        0.285074,
        0.015204,
        -0.003191,
        -0.000146,
        -0.001154,
        -0.005725,
        -0.023312,
        -0.063955,
        -0.149981,
        -0.127262,
        0.544529,
        0.580496,
        1.0,
        0.038109,
        0.20948,
        0.508557,
        0.468842,
        1.0,
        1.0,
        1.0,
        0.006068,
        0.045308,
        0.202822,
        0.503903,
        0.383421,
        1.0,
        1.0,
        1.0,
        1.0,
    ]
    exponent_ref = [
        6665.0,
        1000.0,
        228.0,
        64.71,
        21.06,
        7.495,
        2.797,
        0.5215,
        0.1596,
        6665.0,
        1000.0,
        228.0,
        64.71,
        21.06,
        7.495,
        2.797,
        0.5215,
        0.1596,
        0.1596,
        9.439,
        2.002,
        0.5456,
        0.1517,
        0.1517,
        0.55,
        0.3258,
        33.87,
        5.095,
        1.159,
        0.3258,
        0.1027,
        0.1027,
        1.407,
        0.388,
        1.057,
    ]

    assert nuclei_num_test == nuclei_num_test_ref
    assert nucleus_index_test == nucleus_index_ref
    assert shell_ang_mom_test == shell_ang_mom_ref
    assert shell_ang_mom_turbo_notation_test == shell_ang_mom_turbo_notation_ref
    assert shell_index_test == shell_index_ref
    assert coefficient_test == coefficient_ref
    assert exponent_test == exponent_ref

    assert basis_sets.get_largest_angmom(nucleus_index=0) == 2
    assert basis_sets.get_largest_angmom(nucleus_index=1) == 2

    basis_sets.cut_orbitals(
        thr_exp=1.0,
        nucleus_index=None,
        method="larger",
    )

    nucleus_index_test = basis_sets.nucleus_index
    shell_ang_mom_test = basis_sets.shell_ang_mom
    shell_ang_mom_turbo_notation_test = basis_sets.shell_ang_mom_turbo_notation
    shell_index_test = basis_sets.shell_index
    coefficient_test = basis_sets.coefficient
    exponent_test = basis_sets.exponent

    nucleus_index_ref = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1]
    shell_ang_mom_ref = [0, 0, 0, 1, 1, 2, 0, 0, 0, 1]
    shell_ang_mom_turbo_notation_ref = [300, 300, 16, 400, 36, 37, 16, 300, 16, 36]
    shell_index_ref = [0, 0, 1, 1, 2, 3, 3, 4, 5, 6, 7, 7, 8, 9]
    coefficient_ref = [
        0.015204,
        -0.003191,
        0.544529,
        0.580496,
        1.0,
        0.508557,
        0.468842,
        1.0,
        1.0,
        1.0,
        0.503903,
        0.383421,
        1.0,
        1.0,
    ]
    exponent_ref = [
        0.5215,
        0.1596,
        0.5215,
        0.1596,
        0.1596,
        0.5456,
        0.1517,
        0.1517,
        0.55,
        0.3258,
        0.3258,
        0.1027,
        0.1027,
        0.388,
    ]

    assert shell_ang_mom_ref == shell_ang_mom_test
    assert shell_ang_mom_turbo_notation_ref == shell_ang_mom_turbo_notation_test
    assert shell_index_test == shell_index_ref
    assert coefficient_ref == coefficient_test
    assert exponent_ref == exponent_test

    basis_sets.contracted_to_uncontracted()

    shell_ang_mom_ref = [0, 0, 0, 0, 0, 1, 1, 1, 2, 0, 0, 0, 0, 1]
    shell_ang_mom_turbo_notation_ref = [
        16,
        16,
        16,
        16,
        16,
        36,
        36,
        37,
        16,
        16,
        16,
        16,
        36,
    ]
    coefficient_ref = [
        0.015204,
        -0.003191,
        0.544529,
        0.580496,
        1.0,
        0.508557,
        0.468842,
        1.0,
        1.0,
        1.0,
        0.503903,
        0.383421,
        1.0,
        1.0,
    ]
    exponent_ref = [
        0.5215,
        0.1596,
        0.5215,
        0.1596,
        0.1596,
        0.5456,
        0.1517,
        0.1517,
        0.55,
        0.3258,
        0.3258,
        0.1027,
        0.1027,
        0.388,
    ]

    shell_ang_mom_test = basis_sets.shell_ang_mom
    shell_ang_mom_turbo_notation_test = basis_sets.shell_ang_mom_turbo_notation
    shell_index_test = basis_sets.shell_index
    coefficient_test = basis_sets.coefficient
    exponent_test = basis_sets.exponent

    shell_ang_mom_ref = [0, 0, 1, 1, 2, 0, 0, 1]
    shell_ang_mom_turbo_notation_ref = [16, 16, 36, 36, 37, 16, 16, 36]
    shell_index_ref = [0, 1, 2, 3, 4, 5, 6, 7]
    coefficient_ref = [1.0] * len(shell_index_ref)
    exponent_ref = [0.5215, 0.1596, 0.5456, 0.1517, 0.55, 0.3258, 0.1027, 0.388]

    assert shell_ang_mom_ref == shell_ang_mom_test
    assert shell_ang_mom_turbo_notation_ref == shell_ang_mom_turbo_notation_test
    assert shell_index_ref == shell_index_test
    assert coefficient_ref == coefficient_test
    assert exponent_ref == exponent_test
