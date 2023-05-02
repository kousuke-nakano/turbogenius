#!python
# -*- coding: utf-8 -*-
import os
import shutil
from turbogenius.pyturbo.io_fort10 import IO_fort10

data_dir = os.path.dirname(os.path.abspath(__file__))


def test_fort10_hydrogen():
    shutil.copy(
        os.path.join(data_dir, "fort.10_hydrogen"), os.path.join(data_dir, "fort.10")
    )
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    assert io_fort10.f10structure.ortho_flag
    io_fort10.f10header.io_flag = -1  # => fort.10 is automatically overwritten!!
    assert io_fort10.f10header.io_flag == -1
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    assert io_fort10.f10header.io_flag == -1

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=False)
    assert io_fort10.f10header.io_flag == -1
    io_fort10.f10header.io_flag = -1 * io_fort10.f10header.io_flag
    assert io_fort10.f10header.io_flag == 1

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    assert io_fort10.f10header.io_flag == -1

    assert io_fort10.f10structure.atomic_numbers == [1.09, 1.10]
    assert io_fort10.f10structure.positions[0][2] == -1.00003
    assert io_fort10.f10structure.positions[1][2] == 1.00003
    assert io_fort10.f10forceconstraint.atom_label == [1, -2]
    assert io_fort10.f10forceconstraint.direction == [3, 3]
    assert io_fort10.f10detbasissets.mo_coefficient[0] == [
        -0.5,
        -0.499983966350555,
        -0.389259159564972,
        0.263080060482025,
        -0.320021986961365,
        0.402878105640411,
        0.389840364456177,
        -0.112988829612732,
        -0.02405738830566406,
        0.38219290971756,
        -0.314220011234283,
        0.138799846172333,
        -0.230787396430969,
        -0.253317177295685,
        -0.36382919549942,
        -0.316983282566071,
        0.284486711025238,
        0.216519355773926,
        -0.220531404018402,
        0.118112206459045,
        -0.106018006801605,
        0.443526446819305,
        -0.06911081075668335,
        -0.41402530670166,
        -0.313779950141907,
        -0.465687811374664,
        0.408909916877747,
        0.001605749130249023,
        -0.09116250276565552,
        0.498206079006195,
        -0.403267800807953,
        -0.300212204456329,
        -0.355205833911896,
        0.291875839233398,
        -0.426920294761658,
        0.04561859369277954,
        -0.168644666671753,
        -0.117501437664032,
        0.296707689762115,
        0.305788099765778,
        0.481221675872803,
        -0.498637914657593,
        -0.421001791954041,
        -0.172527253627777,
        -0.28481650352478,
        0.009770333766937256,
        -0.170983791351318,
        0.321721255779266,
        -0.03432732820510864,
        0.0491291880607605,
        0.305401802062988,
        -0.20179158449173,
        0.46058863401413,
        0.397086381912231,
        0.3594731092453,
        0.448230683803558,
        -0.153227925300598,
        -0.298398673534393,
        -0.09608477354049683,
        -0.477181375026703,
    ]
    assert io_fort10.f10detbasissets.exponent == [
        23.843185,
        10.212443,
        4.374164,
        1.873529,
        0.802465,
        0.343709,
        0.147217,
        0.063055,
        0.029292,
        0.091791,
        0.287637,
        0.106105,
        0.393954,
        1.462694,
        0.295883,
        1.065841,
        23.843185,
        10.212443,
        4.374164,
        1.873529,
        0.802465,
        0.343709,
        0.147217,
        0.063055,
        0.029292,
        0.091791,
        0.287637,
        0.106105,
        0.393954,
        1.462694,
        0.295883,
        1.065841,
    ]
    assert io_fort10.f10detbasissets.coefficient == [
        0.004115,
        0.010464,
        0.028011,
        0.075886,
        0.182106,
        0.348521,
        0.378231,
        0.116424,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        0.004115,
        0.010464,
        0.028011,
        0.075886,
        0.182106,
        0.348521,
        0.378231,
        0.116424,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ]
    assert io_fort10.f10jasbasissets.coefficient == [
        0.004115,
        0.010464,
        0.028011,
        0.075886,
        0.182106,
        0.348521,
        0.378231,
        0.116424,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        0.004115,
        0.010464,
        0.028011,
        0.075886,
        0.182106,
        0.348521,
        0.378231,
        0.116424,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    ]

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    mo_coefficient = io_fort10.f10detbasissets.mo_coefficient
    mo_coefficient[0][0] = 1000.0
    io_fort10.f10detbasissets.mo_coefficient = mo_coefficient
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    assert io_fort10.f10detbasissets.mo_coefficient[0][0] == 1000.0

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=False)
    mo_coefficient = io_fort10.f10detbasissets.mo_coefficient
    mo_coefficient[0][0] = -1000.0
    io_fort10.f10detbasissets.mo_coefficient = mo_coefficient
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=False)
    assert io_fort10.f10detbasissets.mo_coefficient[0][0] == 1000.0


def test_fort10_hBN():
    shutil.copy(
        os.path.join(data_dir, "fort.10_hBN"), os.path.join(data_dir, "fort.10")
    )

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    assert not io_fort10.f10structure.ortho_flag
    assert io_fort10.complex_flag
    assert io_fort10.f10detmatrix.coeff_real == [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    assert io_fort10.f10detmatrix.coeff_imag == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    assert io_fort10.f10structure.phase_up == [+0.25, +0.35, +0.45]
    assert io_fort10.f10structure.phase_dn == [-0.25, -0.35, -0.45]

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    mo_coefficient = io_fort10.f10detbasissets.mo_coefficient
    mo_coefficient[0][0] = 1000.0
    io_fort10.f10detbasissets.mo_coefficient = mo_coefficient
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    assert io_fort10.f10detbasissets.mo_coefficient[0][0] == 1000.0

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    mo_coefficient_imag = io_fort10.f10detbasissets.mo_coefficient_imag
    mo_coefficient_imag[0][0] = 5000.0
    io_fort10.f10detbasissets.mo_coefficient_imag = mo_coefficient_imag
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=True)
    assert io_fort10.f10detbasissets.mo_coefficient_imag[0][0] == 5000.0

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=False)
    mo_coefficient = io_fort10.f10detbasissets.mo_coefficient
    mo_coefficient[0][0] = -1000.0
    io_fort10.f10detbasissets.mo_coefficient = mo_coefficient
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=False)
    assert io_fort10.f10detbasissets.mo_coefficient[0][0] == 1000.0

    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=False)
    mo_coefficient_imag = io_fort10.f10detbasissets.mo_coefficient_imag
    mo_coefficient_imag[0][0] = -5000.0
    io_fort10.f10detbasissets.mo_coefficient_imag = mo_coefficient_imag
    io_fort10 = IO_fort10(os.path.join(data_dir, "fort.10"), in_place=False)
    assert io_fort10.f10detbasissets.mo_coefficient_imag[0][0] == 5000.0
