#!python -u
# -*- coding: utf-8 -*-

"""

pyturbo: structure related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# python modules
import os, sys
import math
import numpy as np
from numpy import linalg as LA

# python material modules
from ase import Atoms
from ase.visualize import view
from ase.io import write, read

# set logger
from logging import config, getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)

# turbogenius module
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils.units import Angstrom


class Cell:
    def __init__(
        self,
        vec_a=np.array([0.0, 0.0, 0.0], dtype=float),  # bohr
        vec_b=np.array([0.0, 0.0, 0.0], dtype=float),  # bohr
        vec_c=np.array([0.0, 0.0, 0.0], dtype=float),  # bohr
    ):

        # cell vectors (the units are bohr)
        self.__vec_a = vec_a
        self.__vec_b = vec_b
        self.__vec_c = vec_c

        # calc norm
        self.__norm_vec_a = LA.norm(self.__vec_a)
        self.__norm_vec_b = LA.norm(self.__vec_b)
        self.__norm_vec_c = LA.norm(self.__vec_c)

        if np.sum(self.__vec_a + self.__vec_b + self.__vec_c) == 0.0:
            self.__pbc_flag = False
        else:
            self.__pbc_flag = True

        if self.__pbc_flag == False:
            self.__cos_alpha = np.nan
            self.__cos_beta = np.nan
            self.__cos_gamma = np.nan
            self.__tilted_flag = False
            self.__celldm_1 = 0.0
            self.__celldm_2 = 0.0
            self.__celldm_3 = 0.0
            self.__celldm_4 = 0.0
            self.__celldm_5 = 0.0
            self.__celldm_6 = 0.0

        else:  # self.pbc_flag==True:

            # convert vec_x,y,z to celldm_1 - celldm_6

            # compute cos
            self.__cos_alpha = np.inner(self.__vec_b, self.__vec_c) / (
                self.__norm_vec_b * self.__norm_vec_c
            )
            self.__cos_beta = np.inner(self.__vec_c, self.__vec_a) / (
                self.__norm_vec_c * self.__norm_vec_a
            )
            self.__cos_gamma = np.inner(self.__vec_a, self.__vec_b) / (
                self.__norm_vec_a * self.__norm_vec_b
            )

            # check if the unitcell is orthorhombic
            if (
                self.__cos_alpha != 0.0
                or self.__cos_beta != 0.0
                or self.__cos_gamma != 0.0
            ):
                logger.info(
                    "The unit cell is not orthorhombic. The tilted option is turned on."
                )
                self.__tilted_flag = True
            else:
                logger.info(
                    "The unit cell is orthorhombic. The tilted option is turned off."
                )
                self.__tilted_flag = False

            # set celldm(1-6)
            if self.__tilted_flag:
                self.__celldm_1 = self.__norm_vec_a
                self.__celldm_2 = self.__norm_vec_b / self.__celldm_1
                self.__celldm_3 = self.__norm_vec_c / self.__celldm_1
                self.__celldm_4 = math.acos(self.__cos_alpha)
                self.__celldm_5 = math.acos(self.__cos_beta)
                self.__celldm_6 = math.acos(self.__cos_gamma)

            else:  # self.tilted_flag:
                self.__celldm_1 = self.__norm_vec_a
                self.__celldm_2 = self.__norm_vec_b / self.__celldm_1
                self.__celldm_3 = self.__norm_vec_c / self.__celldm_1
                self.__celldm_4 = 0.0
                self.__celldm_5 = 0.0
                self.__celldm_6 = 0.0

    @property
    def vec_a(self):
        return self.__vec_a

    @property
    def vec_b(self):
        return self.__vec_b

    @property
    def vec_c(self):
        return self.__vec_c

    @property
    def celldm_1(self):
        return self.__celldm_1

    @property
    def celldm_2(self):
        return self.__celldm_2

    @property
    def celldm_3(self):
        return self.__celldm_3

    @property
    def celldm_4(self):
        return self.__celldm_4

    @property
    def celldm_5(self):
        return self.__celldm_5

    @property
    def celldm_6(self):
        return self.__celldm_6

    @property
    def pbc_flag(self):
        return self.__pbc_flag

    @property
    def tilted_flag(self):
        return self.__tilted_flag

    """
    def parse_cell_from_celldm
    Description: Quantum Espresso
    14 Triclinic                    celldm(2)= b/a,
                                    celldm(3)= c/a,
                                    celldm(4)= cos(bc),
                                    celldm(5)= cos(ac),
                                    celldm(6)= cos(ab)
    v1 = (a, 0, 0),
    v2 = (b*cos(gamma), b*sin(gamma), 0)
    v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
        c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
                     - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
    where alpha is the angle between axis b and c
          beta is the angle between axis a and c
          gamma is the angle between axis a and b
    """

    @classmethod
    def parse_cell_from_celldm(
        cls, celldm_1, celldm_2, celldm_3, celldm_4, celldm_5, celldm_6
    ):

        # note! the definition of celldm 4-6 are different from QE
        alpha = celldm_4
        beta = celldm_5
        gamma = celldm_6

        if alpha == 0.0 and beta == 0.0 and gamma == 0.0:
            vec_a = np.array([celldm_1, 0.0, 0.0], dtype=float)

            vec_b = np.array([0.0, celldm_2 * celldm_1, 0.0], dtype=float)

            vec_c = np.array([0.0, 0.0, celldm_3 * celldm_1], dtype=float)

        else:
            vec_a = np.array([celldm_1, 0.0, 0.0], dtype=float)

            b = celldm_2 * celldm_1
            vec_b = np.array(
                [b * math.cos(gamma), b * math.sin(gamma), 0.0], dtype=float
            )

            c = celldm_3 * celldm_1
            vec_c = np.array(
                [
                    c * math.cos(beta),
                    c
                    * (math.cos(alpha) - math.cos(beta) * math.cos(gamma))
                    / math.sin(gamma),
                    c
                    * np.sqrt(
                        1
                        + 2
                        * math.cos(alpha)
                        * math.cos(beta)
                        * math.cos(gamma)
                        - math.cos(alpha) ** 2
                        - math.cos(beta) ** 2
                        - math.cos(gamma) ** 2
                    )
                    / math.sin(gamma),
                ],
                dtype=float,
            )

        return cls(vec_a=vec_a, vec_b=vec_b, vec_c=vec_c)


class Structure:
    def __init__(
        self,
        cell=Cell(),
        atomic_numbers=[],
        element_symbols=[],
        positions=np.array([[]]),  # 3 * N matrix, the unit is bohr!!
    ):

        self.__cell = cell
        self.__atomic_numbers = atomic_numbers
        self.__element_symbols = element_symbols
        self.__positions = positions

    def __str__(self):
        output = [
            "Structure class",
        ]
        return "\n".join(output)

    @property
    def cell(self):
        return self.__cell

    @property
    def atomic_numbers(self):
        return self.__atomic_numbers

    @property
    def element_symbols(self):
        return self.__element_symbols

    @property
    def positions(self):
        return self.__positions

    @property
    def natom(self):
        return len(self.__atomic_numbers)

    @property
    def ntyp(self):
        return len(list(set(self.__atomic_numbers)))

    @property
    def pbc_flag(self):
        return self.__cell.pbc_flag

    @property
    def tilted_flag(self):
        return self.__cell.tilted_flag

    def view(self):
        atom = self.get_ase_atom()
        view(atom)

    def get_ase_atom(self):
        # define ASE-type structure (used inside this class)
        # Note! unit in ASE is angstrom, so one should convert bohr -> ang
        if self.__cell.pbc_flag:
            ase_atom = Atoms(
                self.__element_symbols, positions=self.__positions / Angstrom
            )
            ase_atom.set_cell(
                np.array(
                    [
                        self.__cell.vec_a / Angstrom,
                        self.__cell.vec_b / Angstrom,
                        self.__cell.vec_c / Angstrom,
                    ]
                )
            )
            ase_atom.set_pbc(True)
        else:
            ase_atom = Atoms(
                self.__element_symbols, positions=self.__positions / Angstrom
            )
            ase_atom.set_pbc(False)

        return ase_atom

    @classmethod
    def parse_structure_from_ase_atom(cls, ase_atom):
        if all(ase_atom.get_pbc()):
            vec_a = ase_atom.get_cell()[0] * Angstrom
            vec_b = ase_atom.get_cell()[1] * Angstrom
            vec_c = ase_atom.get_cell()[2] * Angstrom
            cell = Cell(vec_a=vec_a, vec_b=vec_b, vec_c=vec_c)
        else:
            cell = Cell()

        atomic_numbers = ase_atom.get_atomic_numbers()
        element_symbols = ase_atom.get_chemical_symbols()
        positions = ase_atom.get_positions() * Angstrom

        return cls(
            cell=cell,
            atomic_numbers=atomic_numbers,
            element_symbols=element_symbols,
            positions=positions,
        )

    @classmethod
    def parse_structure_from_file(cls, file):
        logger.info(
            f"Structure is read from {file} using the ASE read function."
        )
        atoms = read(file)
        return cls.parse_structure_from_ase_atom(atoms)

    def write(self, file):
        atoms = self.get_ase_atom()
        write(file, atoms)


if __name__ == "__main__":
    log = getLogger("pyturbo")
    log.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    log.addHandler(stream_handler)

    # moved to examples
