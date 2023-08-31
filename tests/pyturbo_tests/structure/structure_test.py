#!python
# -*- coding: utf-8 -*-
import os, sys
import numpy as np
from ase.spacegroup import crystal

# turbogenius modules
from turbogenius.pyturbo.structure import Cell
from turbogenius.pyturbo.structure import Structure

# unit conversion
ref_BOHR = 0.529177210903

# H2.xyz
structure = Structure.parse_structure_from_file(file="H2.xyz")
ref_atomic_numbers = [1, 1]
np.testing.assert_array_equal(structure.atomic_numbers, ref_atomic_numbers)
ref_vec_a = [0.0, 0.0, 0.0]
ref_vec_b = [0.0, 0.0, 0.0]
ref_vec_c = [0.0, 0.0, 0.0]
np.testing.assert_array_equal(structure.cell.vec_a, ref_vec_a)
np.testing.assert_array_equal(structure.cell.vec_b, ref_vec_b)
np.testing.assert_array_equal(structure.cell.vec_c, ref_vec_c)
assert not structure.pbc_flag
assert not structure.tilted_flag
ref_celldm_1 = 0.0
ref_celldm_2 = 0.0
ref_celldm_3 = 0.0
ref_celldm_4 = 0.0
ref_celldm_5 = 0.0
ref_celldm_6 = 0.0
np.testing.assert_array_almost_equal(
    structure.cell.celldm_1, ref_celldm_1, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_2, ref_celldm_2, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_3, ref_celldm_3, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_4, ref_celldm_4, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_5, ref_celldm_5, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_6, ref_celldm_6, decimal=15
)

# Diamond.vasp
structure = Structure.parse_structure_from_file(file="Diamond.vasp")
ref_atomic_numbers = [6, 6, 6, 6, 6, 6, 6, 6]
np.testing.assert_array_equal(structure.atomic_numbers, ref_atomic_numbers)
ref_l = 3.5669999123
ref_vec_a = [ref_l / ref_BOHR, 0.0, 0.0]
ref_vec_b = [0.0, ref_l / ref_BOHR, 0.0]
ref_vec_c = [0.0, 0.0, ref_l / ref_BOHR]
np.testing.assert_array_equal(structure.cell.vec_a, ref_vec_a)
np.testing.assert_array_equal(structure.cell.vec_b, ref_vec_b)
np.testing.assert_array_equal(structure.cell.vec_c, ref_vec_c)
assert structure.pbc_flag
assert not structure.tilted_flag
ref_celldm_1 = ref_l / ref_BOHR
ref_celldm_2 = 1.0
ref_celldm_3 = 1.0
ref_celldm_4 = np.pi / 2  # radian
ref_celldm_5 = np.pi / 2  # radian
ref_celldm_6 = np.pi / 2  # radian
np.testing.assert_array_almost_equal(
    structure.cell.celldm_1, ref_celldm_1, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_2, ref_celldm_2, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_3, ref_celldm_3, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_4, ref_celldm_4, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_5, ref_celldm_5, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_6, ref_celldm_6, decimal=15
)

# Diamond.cif
structure = Structure.parse_structure_from_file(file="Diamond.cif")
ref_atomic_numbers = [6, 6, 6, 6, 6, 6, 6, 6]
np.testing.assert_array_equal(structure.atomic_numbers, ref_atomic_numbers)
ref_l = 3.567
ref_vec_a = [ref_l / ref_BOHR, 0.0, 0.0]
ref_vec_b = [0.0, ref_l / ref_BOHR, 0.0]
ref_vec_c = [0.0, 0.0, ref_l / ref_BOHR]
np.testing.assert_array_equal(structure.cell.vec_a, ref_vec_a)
np.testing.assert_array_equal(structure.cell.vec_b, ref_vec_b)
np.testing.assert_array_equal(structure.cell.vec_c, ref_vec_c)
assert structure.pbc_flag
assert not structure.tilted_flag
ref_celldm_1 = ref_l / ref_BOHR
ref_celldm_2 = 1.0
ref_celldm_3 = 1.0
ref_celldm_4 = np.pi / 2  # radian
ref_celldm_5 = np.pi / 2  # radian
ref_celldm_6 = np.pi / 2  # radian
np.testing.assert_array_almost_equal(
    structure.cell.celldm_1, ref_celldm_1, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_2, ref_celldm_2, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_3, ref_celldm_3, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_4, ref_celldm_4, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_5, ref_celldm_5, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_6, ref_celldm_6, decimal=15
)

# benzene.xyz
structure = Structure.parse_structure_from_file(file="benzene.xyz")
ref_atomic_numbers = [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1]
np.testing.assert_array_equal(structure.atomic_numbers, ref_atomic_numbers)
ref_vec_a = [0.0, 0.0, 0.0]
ref_vec_b = [0.0, 0.0, 0.0]
ref_vec_c = [0.0, 0.0, 0.0]
np.testing.assert_array_equal(structure.cell.vec_a, ref_vec_a)
np.testing.assert_array_equal(structure.cell.vec_b, ref_vec_b)
np.testing.assert_array_equal(structure.cell.vec_c, ref_vec_c)
assert not structure.pbc_flag
assert not structure.tilted_flag
ref_celldm_1 = 0.0
ref_celldm_2 = 0.0
ref_celldm_3 = 0.0
ref_celldm_4 = 0.0
ref_celldm_5 = 0.0
ref_celldm_6 = 0.0
np.testing.assert_array_almost_equal(
    structure.cell.celldm_1, ref_celldm_1, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_2, ref_celldm_2, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_3, ref_celldm_3, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_4, ref_celldm_4, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_5, ref_celldm_5, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_6, ref_celldm_6, decimal=15
)

# silicon_oxide.cif
structure = Structure.parse_structure_from_file(file="silicon_oxide.cif")
ref_atomic_numbers = [14, 14, 14, 8, 8, 8, 8, 8, 8]
np.testing.assert_array_equal(structure.atomic_numbers, ref_atomic_numbers)
ref_vec_a = [4.705 / ref_BOHR, 0.0, 0.0]
ref_vec_b = [
    4.705 * np.cos(2 * np.pi / 3) / ref_BOHR,
    4.705 * np.sin(2 * np.pi / 3) / ref_BOHR,
    0.0,
]
ref_vec_c = [0.0, 0.0, 5.250 / ref_BOHR]
np.testing.assert_array_equal(structure.cell.vec_a, ref_vec_a)
np.testing.assert_array_almost_equal(
    structure.cell.vec_b, ref_vec_b, decimal=15
)
np.testing.assert_array_equal(structure.cell.vec_c, ref_vec_c)
assert structure.pbc_flag
assert structure.tilted_flag
ref_celldm_1 = 4.705 / ref_BOHR
ref_celldm_2 = 1.0
ref_celldm_3 = 5.250 / 4.705
ref_celldm_4 = np.pi / 2  # radian
ref_celldm_5 = np.pi / 2  # radian
ref_celldm_6 = 2 * np.pi / 3  # radian
np.testing.assert_array_almost_equal(
    structure.cell.celldm_1, ref_celldm_1, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_2, ref_celldm_2, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_3, ref_celldm_3, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_4, ref_celldm_4, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_5, ref_celldm_5, decimal=15
)
np.testing.assert_array_almost_equal(
    structure.cell.celldm_6, ref_celldm_6, decimal=15
)

# bulk Al (from ase instance)
a = 4.05
atom = crystal(
    "Al", [(0, 0, 0)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90]
)
structure = Structure.parse_structure_from_ase_atom(atom)
ref_atomic_numbers = [13, 13, 13, 13]
np.testing.assert_array_equal(structure.atomic_numbers, ref_atomic_numbers)
ref_vec_a = [a / ref_BOHR, 0.0, 0.0]
ref_vec_b = [0.0, a / ref_BOHR, 0.0]
ref_vec_c = [0.0, 0.0, a / ref_BOHR]
np.testing.assert_array_equal(structure.cell.vec_a, ref_vec_a)
np.testing.assert_array_equal(structure.cell.vec_b, ref_vec_b)
np.testing.assert_array_equal(structure.cell.vec_c, ref_vec_c)
assert structure.pbc_flag
assert not structure.tilted_flag

# celldm test (hexagonal, silicon oxide)
celldm_1 = 4.705 / ref_BOHR
celldm_2 = 1.0
celldm_3 = 5.250 / 4.705
celldm_4 = np.pi / 2  # radian
celldm_5 = np.pi / 2  # radian
celldm_6 = 2 * np.pi / 3  # radian
cell = Cell.parse_cell_from_celldm(
    celldm_1=celldm_1,
    celldm_2=celldm_2,
    celldm_3=celldm_3,
    celldm_4=celldm_4,
    celldm_5=celldm_5,
    celldm_6=celldm_6,
)
ref_vec_a = [4.705 / ref_BOHR, 0.0, 0.0]
ref_vec_b = [
    4.705 * np.cos(2 * np.pi / 3) / ref_BOHR,
    4.705 * np.sin(2 * np.pi / 3) / ref_BOHR,
    0.0,
]
ref_vec_c = [0.0, 0.0, 5.250 / ref_BOHR]

np.testing.assert_array_almost_equal(cell.vec_a, ref_vec_a, decimal=15)
np.testing.assert_array_almost_equal(cell.vec_b, ref_vec_b, decimal=15)
np.testing.assert_array_almost_equal(cell.vec_c, ref_vec_c, decimal=15)
