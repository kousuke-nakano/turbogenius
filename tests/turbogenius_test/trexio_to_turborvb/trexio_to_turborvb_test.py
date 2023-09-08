#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil

import numpy as np

from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.pyturbo.basis_set import Jas_Basis_sets
from turbogenius.trexio_to_turborvb import trexio_to_turborvb_wf

ref_BOHR = 0.529177210903

root_dir = os.path.dirname(__file__)


# H2
def test_trexio_converter_H2():
    os.chdir(root_dir)
    trexio_filename = "H2_trexio.hdf5"
    prefix, _ = os.path.splitext(trexio_filename)
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, trexio_filename),
        os.path.join(root_dir, prefix, trexio_filename),
    )
    os.chdir(os.path.join(root_dir, prefix))

    # Jastrow basis (GAMESS format)
    H_jastrow_basis = """
            S  1
            1       1.873529  1.00000000
            S  1
            1       0.343709  1.00000000
            S  1
            1       0.139013  1.00000000
            P  1
            1       0.740212  1.00000000
    """

    jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_texts(
        [H_jastrow_basis, H_jastrow_basis], format="gamess"
    )

    # convert the TREXIO file to
    trexio_to_turborvb_wf(
        trexio_file=trexio_filename,
        jas_basis_sets=jas_basis_sets,
        max_occ_conv=0,
        mo_num_conv=-1,
        only_mol=True,
        cleanup=True,
    )

    # assertions!
    fort10 = IO_fort10("fort.10")
    assert fort10.pp_flag
    assert not fort10.pbc_flag
    assert fort10.det_contraction_flag
    assert not fort10.complex_flag
    assert fort10.ansatz_type == "agpn"

    assert fort10.f10header.nelup == 1
    assert fort10.f10header.neldn == 1
    assert fort10.f10header.nel == 2
    assert fort10.f10header.natom == 2

    np.testing.assert_array_almost_equal(
        np.abs(
            fort10.f10structure.positions[0][2] - fort10.f10structure.positions[1][2]
        ),
        0.37000000 * 2 / ref_BOHR,
        decimal=10,
    )

    os.chdir(root_dir)


# wBN
def test_trexio_converter_wBN():
    os.chdir(root_dir)
    trexio_filename = "wBN_trexio.hdf5"
    prefix, _ = os.path.splitext(trexio_filename)
    if os.path.isdir(os.path.join(root_dir, prefix)):
        shutil.rmtree(os.path.join(root_dir, prefix))
    os.makedirs(os.path.join(root_dir, prefix))
    shutil.copy(
        os.path.join(root_dir, trexio_filename),
        os.path.join(root_dir, prefix, trexio_filename),
    )
    os.chdir(os.path.join(root_dir, prefix))

    # Jastrow basis (GAMESS format)
    wBN_jastrow_basis = """
            S  1
            1       1.873529  1.00000000
            S  1
            1       0.343709  1.00000000
            S  1
            1       0.139013  1.00000000
            P  1
            1       0.740212  1.00000000
    """

    jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_texts(
        [wBN_jastrow_basis, wBN_jastrow_basis], format="gamess"
    )

    # convert the TREXIO file to
    trexio_to_turborvb_wf(
        trexio_file=trexio_filename,
        jas_basis_sets=jas_basis_sets,
        max_occ_conv=0,
        mo_num_conv=-1,
        only_mol=True,
        cleanup=True,
    )

    # assertions!
    fort10 = IO_fort10("fort.10")
    assert fort10.pp_flag
    assert fort10.pbc_flag
    assert fort10.det_contraction_flag
    assert fort10.complex_flag
    assert fort10.ansatz_type == "agpn"

    assert fort10.f10header.nelup == 8
    assert fort10.f10header.neldn == 8
    assert fort10.f10header.nel == 16
    assert fort10.f10header.natom == 4

    ref_vec_a = [2.536 / ref_BOHR, 0.0, 0.0]
    ref_vec_b = [
        2.536 * np.cos(2 * np.pi / 3) / ref_BOHR,
        2.536 * np.sin(2 * np.pi / 3) / ref_BOHR,
        0.0,
    ]
    ref_vec_c = [0.0, 0.0, 4.199 / ref_BOHR]
    np.testing.assert_array_almost_equal(
        fort10.f10structure.vec_a, ref_vec_a, decimal=8
    )
    np.testing.assert_array_almost_equal(
        fort10.f10structure.vec_b, ref_vec_b, decimal=8
    )
    np.testing.assert_array_almost_equal(
        fort10.f10structure.vec_c, ref_vec_c, decimal=8
    )

    os.chdir(root_dir)
