# python
# -*- coding: utf-8 -*-
"""
pyturbo: io_fort10 related classes and methods, why??
Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not
    * be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.
"""

# python modules
import os
import sys
import numpy as np
import math
import time
import shutil
import psutil
import linecache
from typing import Union
from tqdm import tqdm

# set logger
import io
import logging
from logging import getLogger, StreamHandler, Formatter

# pyturbo module
from turbogenius.pyturbo.utils.utility import (
    pygrep_lineno,
    pysed_replace,
    pygetline,
    pysed_replace_lines,
)
from turbogenius.pyturbo.utils.utility import (
    return_orb_type_chr,
    return_contraction_flag,
    return_num_twobody_and_flag_onebody,
)
from turbogenius.pyturbo.utils.utility import return_element_symbol
from turbogenius.pyturbo.structure import Structure, Cell
from turbogenius.pyturbo.basis_set import Det_Basis_sets, Jas_Basis_sets

# set logger
logger = getLogger("pyturbo").getChild(__name__)


class TqdmToLogger(io.StringIO):
    """
    Output stream for TQDM which will output to logger module instead of
    the StdOut.
    """

    logger = None
    level = None
    buf = ""

    def __init__(self, logger, level=None):
        super(TqdmToLogger, self).__init__()
        self.logger = logger
        self.level = level or logging.INFO

    def write(self, buf):
        self.buf = buf.strip("\r\n\t ")

    def flush(self):
        self.logger.log(self.level, self.buf)


class Value:
    def __init__(self, value=None, lineno=0, index=0, file=None):
        self.__value = value
        self.__lineno = lineno
        self.__index = index
        self.__file = file

    def replace(self, value: Union[str, int, float], in_place: bool = True):
        self.__value = value
        if in_place:
            pysed_replace(self.__file, value, self.__lineno, self.__index)

    def type(self, cast):
        self.__value = cast(self.__value)

    @property
    def v(self):
        return self.__value

    @property
    def l(self):
        return self.__lineno

    @property
    def i(self):
        return self.__index

    @property
    def f(self):
        return self.__file


class IO_fort10:
    """
    This class is a wrapper for python fort.10 file

    Input:
        fort.10 (str): the name of fort.10 WF file
        in_place (bool): if True, fort.10 file is updated
        whenever fort10 instance is updated.

    Attributes:
         fort.10 (str): the name of fort.10 WF file
         f10header (F10header): F10header instance
         ...
    """

    f10structure_start_keyword = "Ion coordinates"
    f10structure_end_keyword = "Constraints for forces"
    f10forceconstraint_start_keyword = "Constraints for forces"
    f10forceconstraint_end_keyword = "Parameters Jastrow two body *$"
    f10jastwobody_start_keyword = "Parameters Jastrow two body"
    f10jastwobody_end_keyword = "Parameters atomic wf"
    f10detbasis_start_keyword = "Parameters atomic wf *$"
    f10detbasis_end_keyword = "Parameters atomic Jastrow wf *$"
    f10jasbasis_start_keyword = "Parameters atomic Jastrow wf *$"
    f10jasbasis_end_keyword = "Occupation atomic orbitals *$"
    f10detocc_start_keyword = "Occupation atomic orbitals *$"
    f10detocc_end_keyword = "Occupation atomic orbitals  Jastrow *$"
    f10jasocc_start_keyword = "Occupation atomic orbitals  Jastrow *$"
    f10jasocc_end_keyword = "Nonzero values of  detmat *$"
    f10detmat_start_keyword = "Nonzero values of  detmat *$"
    f10detmat_end_keyword = "Grouped par.  in the chosen ordered basis *$"
    f10detmatsym_start_keyword = "Grouped par.  in the chosen ordered basis *$"
    f10detmatsym_end_keyword = "Nonzero values of  jasmat *$"
    f10jasmat_start_keyword = "Nonzero values of  jasmat *$"
    f10jasmat_end_keyword = "Eq. par. in the 3-body Jastrow in the chosen basis *$"
    f10jasmatsym_start_keyword = "Eq. par. in the 3-body Jastrow in the chosen basis *$"
    f10jasmatsym_end_keyword = "Eq. par. in the atomic Det par.in the chosen basis *$"
    f10detbasis_sym_start_keyword = (
        "Eq. par. in the atomic Det par.in the chosen basis *$"
    )
    f10detbasis_sym_end_keyword = (
        "Eq. par. in the atomic 3-body  par. in the chosen basis *$"
    )
    f10jasbasis_sym_start_keyword = (
        "Eq. par. in the atomic 3-body  par. in the chosen basis *$"
    )
    f10jasbasis_sym_end_keyword = "New parameters *$"

    def __init__(self, fort10: str = "fort.10", in_place: bool = True):

        self.fort10 = fort10
        self.in_place = in_place
        self.f10header = F10header(fort10=self.fort10, in_place=self.in_place)
        # logger.debug("header")
        self.f10structure = F10structure(
            fort10=self.fort10,
            Nel=self.f10header.nel,
            Ion=self.f10header.natom,
            start_keyword=self.f10structure_start_keyword,
            end_keyword=self.f10structure_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("structure")
        self.f10forceconstraint = F10forceconstraint(
            fort10=self.fort10,
            ieskinr=self.f10header.ieskinr,
            start_keyword=self.f10forceconstraint_start_keyword,
            end_keyword=self.f10forceconstraint_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10forceconstraint")
        self.f10jastwobody = F10jastwobody(
            fort10=self.fort10,
            jas_2body=self.f10header.jas_2body,
            start_keyword=self.f10jastwobody_start_keyword,
            end_keyword=self.f10jastwobody_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10jastwobody")
        self.f10detbasissets = F10detbasissets(
            fort10=self.fort10,
            shell_det=self.f10header.shell_det,
            complex_flag=self.f10header.complex_flag,
            start_keyword=self.f10detbasis_start_keyword,
            end_keyword=self.f10detbasis_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10detbasissets")
        self.f10jasbasissets = F10jasbasissets(
            fort10=self.fort10,
            shell_jas=self.f10header.shell_jas,
            start_keyword=self.f10jasbasis_start_keyword,
            end_keyword=self.f10jasbasis_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10jasbasissets")
        self.f10detocc = F10occ(
            fort10=self.fort10,
            start_keyword=self.f10detocc_start_keyword,
            end_keyword=self.f10detocc_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10detocc")
        self.f10jasocc = F10occ(
            fort10=self.fort10,
            start_keyword=self.f10jasocc_start_keyword,
            end_keyword=self.f10jasocc_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10jasocc")
        self.f10detmatrix = F10detmat(
            fort10=self.fort10,
            detmat=self.f10header.det_mat_nonzero,
            complex_flag=self.f10header.complex_flag,
            start_keyword=self.f10detmat_start_keyword,
            end_keyword=self.f10jasmat_start_keyword,
            in_place=self.in_place,
        )
        # logger.debug("ff10detmatrix")
        self.f10detmat_sym = F10matsymmetry(
            fort10=self.fort10,
            num_component=self.f10header.iessw,
            start_keyword=self.f10detmatsym_start_keyword,
            end_keyword=self.f10detmatsym_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10detmat_sym")
        self.f10jasmatrix = F10jasmat(
            fort10=self.fort10,
            jasmat=self.f10header.jas_mat_nonzero,
            start_keyword=self.f10jasmat_start_keyword,
            end_keyword=self.f10jasmat_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10jasmatrix")
        self.f10jasmat_sym = F10matsymmetry(
            fort10=self.fort10,
            num_component=self.f10header.iesfree,
            start_keyword=self.f10jasmatsym_start_keyword,
            end_keyword=self.f10jasmatsym_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10jasmat_sym")
        self.f10detbasis_sym = F10basissymmetry(
            fort10=self.fort10,
            num_component=self.f10header.eq_det_atomic_par,
            start_keyword=self.f10detbasis_sym_start_keyword,
            end_keyword=self.f10detbasis_sym_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10detbasis_sym")
        self.f10jasbasis_sym = F10basissymmetry(
            fort10=self.fort10,
            num_component=self.f10header.eq_3_body_atomic_par,
            start_keyword=self.f10jasbasis_sym_start_keyword,
            end_keyword=self.f10jasbasis_sym_end_keyword,
            in_place=self.in_place,
        )
        # logger.debug("f10jasbasis_sym")
        # logger.debug("Init End")

    # properties!!
    @property
    def pp_flag(self):
        # logger.debug("pseudo_check")
        atomic_numbers = self.f10structure.atomic_numbers
        valence_electrons = self.f10structure.valence_electrons
        diff_list = np.array(atomic_numbers) - np.array(valence_electrons)
        if np.sum(diff_list) != 0.0:  # PP at least one element:
            return True
        else:  # all-electron
            return False

    @property
    def pbc_flag(self):
        return self.f10structure.pbc_flag

    @property
    def det_contraction_flag(self):
        logger.debug("contraction_check")
        shell_index = self.f10detbasissets.shell_index
        if len(shell_index) != len(list(set(shell_index))):
            return True
        else:
            return False

    @property
    def complex_flag(self):
        return self.f10header.complex_flag

    @property
    def ansatz_type(self):
        ansatz_list = ["sd", "agps", "agpu", "agpn", "pf", "pfn"]
        if self.f10header.nel > 0:
            # agp
            if self.f10detbasissets.has_mo:
                if self.f10structure.pbc_flag:
                    if self.f10detbasissets.num_mo == self.f10header.neldn * 2 + (
                        self.f10header.nelup - self.f10header.neldn
                    ):
                        ansatz = "sd"
                    else:
                        ansatz = "agpn"
                else:
                    if self.f10detbasissets.num_mo == self.f10header.neldn * 1 + (
                        self.f10header.nelup - self.f10header.neldn
                    ):
                        ansatz = "sd"
                    else:
                        ansatz = "agpn"
            else:
                if self.f10jastwobody.num_jas_param < 0:
                    ansatz = "agpu"
                else:
                    ansatz = "agps"
        else:
            # pf
            if self.f10detbasissets.has_mo:
                ansatz = "pfn"
            else:
                ansatz = "pf"

        if ansatz not in ansatz_list:
            raise ValueError

        return ansatz


class F10header:
    def __init__(self, fort10: str, in_place: bool = True):
        self.start_keyword = "Nelup"
        self.end_keyword = "Ion coordinates"
        self.fort10 = fort10
        self.in_place = in_place
        self.read_flag = False

        # Value()
        self.__Nelup = Value()
        self.__Nel = Value()
        self.__Ion = Value()
        self.__Shell_Det = Value()
        self.__Shell_Jas = Value()
        self.__Jas_2body = Value()
        self.__Det = Value()
        self.__three_body_atomic_par = Value()
        self.__Det_mat_nonzero = Value()
        self.__Jas_mat_nonzero = Value()
        self.__Eq_Det_atomic_par = Value()
        self.__Eq_3_body_atomic_par = Value()
        self.__iesfree = Value()
        self.__iessw = Value()
        self.__ieskinr = Value()
        self.__io_flag = Value()

        # auxiliary values
        self.complex_flag = False
        self.crystal_jastrow_flag = True
        self.tilted_flag = False
        self.pfaff_flag = True

    def read(self):
        if not self.read_flag:
            with open(self.fort10, "r"):
                lineno = self.start_lineno
                line = pygetline(filename=self.fort10, lineno=lineno).split()
                self.__Nelup = Value(
                    value=int(line[0]),
                    lineno=lineno,
                    index=0,
                    file=self.fort10,
                )
                self.__Nel = Value(
                    value=int(line[1]),
                    lineno=lineno,
                    index=1,
                    file=self.fort10,
                )
                self.__Ion = Value(
                    value=int(line[2]),
                    lineno=lineno,
                    index=2,
                    file=self.fort10,
                )
                lineno += 2
                line = pygetline(filename=self.fort10, lineno=lineno).split()
                self.__Shell_Det = Value(
                    value=int(line[0]),
                    lineno=lineno,
                    index=0,
                    file=self.fort10,
                )
                self.__Shell_Jas = Value(
                    value=int(line[1]),
                    lineno=lineno,
                    index=1,
                    file=self.fort10,
                )
                lineno += 2
                line = pygetline(filename=self.fort10, lineno=lineno).split()
                self.__Jas_2body = Value(
                    value=int(line[0]),
                    lineno=lineno,
                    index=0,
                    file=self.fort10,
                )
                self.__Det = Value(
                    value=int(line[1]),
                    lineno=lineno,
                    index=1,
                    file=self.fort10,
                )
                self.__three_body_atomic_par = Value(
                    value=int(line[2]),
                    lineno=lineno,
                    index=2,
                    file=self.fort10,
                )
                lineno += 2
                line = pygetline(filename=self.fort10, lineno=lineno).split()
                self.__Det_mat_nonzero = Value(
                    value=int(line[0]),
                    lineno=lineno,
                    index=0,
                    file=self.fort10,
                )
                self.__Jas_mat_nonzero = Value(
                    value=int(line[1]),
                    lineno=lineno,
                    index=1,
                    file=self.fort10,
                )
                lineno += 2
                line = pygetline(filename=self.fort10, lineno=lineno).split()
                self.__Eq_Det_atomic_par = Value(
                    value=int(line[0]),
                    lineno=lineno,
                    index=0,
                    file=self.fort10,
                )
                self.__Eq_3_body_atomic_par = Value(
                    value=int(line[1]),
                    lineno=lineno,
                    index=1,
                    file=self.fort10,
                )
                lineno += 2
                line = pygetline(filename=self.fort10, lineno=lineno).split()
                self.__iesfree = Value(
                    value=int(line[0]),
                    lineno=lineno,
                    index=0,
                    file=self.fort10,
                )
                self.__iessw = Value(
                    value=int(line[1]),
                    lineno=lineno,
                    index=1,
                    file=self.fort10,
                )
                self.__ieskinr = Value(
                    value=int(line[2]),
                    lineno=lineno,
                    index=2,
                    file=self.fort10,
                )
                self.__io_flag = Value(
                    value=int(line[3]),
                    lineno=lineno,
                    index=3,
                    file=self.fort10,
                )

            self.read_flag = True

            # logger.debug("check----1")
            # logger.debug(self.nelup)
            # logger.debug(self.neldn)
            # logger.debug(self.nel)
            # logger.debug(self.natom)
            # logger.debug(self.shell_det)
            # logger.debug(self.shell_jas)
            # logger.debug(self.jas_2body)
            # logger.debug(self.det)
            # logger.debug(self.three_body_atomic_par)
            # logger.debug(self.det_mat_nonzero)
            # logger.debug(self.jas_mat_nonzero)
            # logger.debug(self.eq_det_atomic_par)
            # logger.debug(self.eq_3_body_atomic_par)
            # logger.debug(self.iesfree)
            # logger.debug(self.iessw)
            # logger.debug(self.ieskinr)
            # logger.debug("check-1")
            # logger.debug(self.io_flag)

            # set auxiliary values
            # logger.debug("check0")
            if self.nel < 0:
                self.pfaff_flag = True
            else:
                self.pfaff_flag = False
            # logger.debug("check1")
            if self.shell_det < 0:
                self.complex_flag = True
            else:
                self.complex_flag = False
            # logger.debug("check2")
            if self.shell_jas < 0:
                self.crystal_jastrow_flag = True
            else:
                self.crystal_jastrow_flag = False
            # logger.debug("check3")

            # logger.debug("End init of IO_fort10")

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1

    @property
    def nelup(self):
        self.read()
        return np.abs(self.__Nelup.v)

    @nelup.setter
    def nelup(self, value):
        self.read()
        self.__Nelup.replace(value=value, in_place=self.in_place)

    @property
    def nel(self):
        self.read()
        return np.abs(self.__Nel.v)

    @nel.setter
    def nel(self, value):
        self.read()
        self.__Nel.replace(value=value, in_place=self.in_place)

    @property
    def natom(self):
        self.read()
        return self.__Ion.v

    @natom.setter
    def natom(self, value):
        self.read()
        self.__Ion.replace(value=value, in_place=self.in_place)

    @property
    def shell_det(self):
        self.read()
        return self.__Shell_Det.v

    @shell_det.setter
    def shell_det(self, value):
        self.read()
        self.__Shell_Det.replace(value=value, in_place=self.in_place)

    @property
    def shell_jas(self):
        self.read()
        return self.__Shell_Jas.v

    @shell_jas.setter
    def shell_jas(self, value):
        self.read()
        self.__Shell_Jas.replace(value=value, in_place=self.in_place)

    @property
    def jas_2body(self):
        self.read()
        return self.__Jas_2body.v

    @jas_2body.setter
    def jas_2body(self, value):
        self.read()
        self.__Jas_2body.replace(value=value, in_place=self.in_place)

    @property
    def det(self):
        self.read()
        return self.__Det.v

    @det.setter
    def det(self, value):
        self.read()
        self.__Det.replace(value=value, in_place=self.in_place)

    @property
    def three_body_atomic_par(self):
        self.read()
        return self.__three_body_atomic_par.v

    @three_body_atomic_par.setter
    def three_body_atomic_par(self, value):
        self.read()
        self.__three_body_atomic_par.replace(value=value, in_place=self.in_place)

    @property
    def det_mat_nonzero(self):
        self.read()
        return self.__Det_mat_nonzero.v

    @det_mat_nonzero.setter
    def det_mat_nonzero(self, value):
        self.read()
        self.__Det_mat_nonzero.replace(value=value, in_place=self.in_place)

    @property
    def jas_mat_nonzero(self):
        self.read()
        return self.__Jas_mat_nonzero.v

    @jas_mat_nonzero.setter
    def jas_mat_nonzero(self, value):
        self.read()
        self.__Jas_mat_nonzero.replace(value=value, in_place=self.in_place)

    @property
    def eq_det_atomic_par(self):
        self.read()
        return self.__Eq_Det_atomic_par.v

    @eq_det_atomic_par.setter
    def eq_det_atomic_par(self, value):
        self.read()
        self.__Eq_Det_atomic_par.replace(value=value, in_place=self.in_place)

    @property
    def eq_3_body_atomic_par(self):
        self.read()
        return self.__Eq_3_body_atomic_par.v

    @eq_3_body_atomic_par.setter
    def eq_3_body_atomic_par(self, value):
        self.read()
        self.__Eq_3_body_atomic_par.replace(value=value, in_place=self.in_place)

    @property
    def iesfree(self):
        self.read()
        return self.__iesfree.v

    @iesfree.setter
    def iesfree(self, value):
        self.read()
        self.__iesfree.replace(value=value, in_place=self.in_place)

    @property
    def iessw(self):
        self.read()
        return self.__iessw.v

    @iessw.setter
    def iessw(self, value):
        self.read()
        self.__iessw.replace(value=value, in_place=self.in_place)

    @property
    def ieskinr(self):
        self.read()
        return self.__ieskinr.v

    @ieskinr.setter
    def ieskinr(self, value):
        self.read()
        self.__ieskinr.replace(value=value, in_place=self.in_place)

    @property
    def io_flag(self):
        self.read()
        return self.__io_flag.v

    @io_flag.setter
    def io_flag(self, value):
        self.read()
        self.__io_flag.replace(value=value, in_place=self.in_place)

    # auxilary properties
    @property
    def neldn(self):
        return np.abs(self.__Nel.v) - np.abs(self.__Nelup.v)


class F10structure:
    def __init__(
        self,
        fort10: str,
        Nel: int,
        Ion: int,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.nel = Nel
        self.ion = Ion
        self.in_place = in_place
        self.read_flag = False
        self.__pbc = False
        self.__tilted_flag = False

        self.__atomic_numbers = []
        self.__valence_electrons = []
        self.__positions = [[]]

    def read(self):
        if not self.read_flag:
            # read PBC flag
            first_line = pygetline(filename=self.fort10, lineno=0)
            if "PBC_C" in first_line:
                self.__pbc = True
                self.__tilted_flag = False
            elif "PBC_T" in first_line:
                self.__pbc = True
                self.__tilted_flag = True
            else:
                self.__pbc = False
                self.__tilted_flag = False

            # read cell parameters only for PBC case
            if self.__pbc:
                if not self.__tilted_flag:
                    lineno = 1
                    line = pygetline(filename=self.fort10, lineno=lineno).split()
                    self.__r_s = Value(
                        value=float(line[0]),
                        lineno=lineno,
                        index=0,
                        file=self.fort10,
                    )  # electron density parameter
                    self.__Ly_Lx = Value(
                        value=float(line[1]),
                        lineno=lineno,
                        index=1,
                        file=self.fort10,
                    )  # Ly/Lx
                    self.__Lz_Lx = Value(
                        value=float(line[2]),
                        lineno=lineno,
                        index=2,
                        file=self.fort10,
                    )  # Lz/Lx
                    self.__phase_up_1 = Value(
                        value=float(line[3]),
                        lineno=lineno,
                        index=3,
                        file=self.fort10,
                    )  # phase 1 up
                    self.__phase_up_2 = Value(
                        value=float(line[4]),
                        lineno=lineno,
                        index=4,
                        file=self.fort10,
                    )  # phase 2 up
                    self.__phase_up_3 = Value(
                        value=float(line[5]),
                        lineno=lineno,
                        index=5,
                        file=self.fort10,
                    )  # phase 3 up
                    self.__phase_dn_1 = Value(
                        value=float(line[6]),
                        lineno=lineno,
                        index=6,
                        file=self.fort10,
                    )  # phase 1 dn
                    self.__phase_dn_2 = Value(
                        value=float(line[7]),
                        lineno=lineno,
                        index=7,
                        file=self.fort10,
                    )  # phase 2 dn
                    self.__phase_dn_3 = Value(
                        value=float(line[8]),
                        lineno=lineno,
                        index=8,
                        file=self.fort10,
                    )  # phase 3 dn

                    # convert
                    # these units are bohr
                    x_length = self.__r_s.v * pow(
                        (4 * math.pi * self.nel)
                        / (3 * self.__Ly_Lx.v * self.__Lz_Lx.v),
                        1.0 / 3.0,
                    )
                    y_length = x_length * self.__Ly_Lx.v
                    z_length = x_length * self.__Lz_Lx.v

                    self.__vec_a = np.array([x_length, 0.0, 0.0], dtype=float)
                    self.__vec_b = np.array([0.0, y_length, 0.0], dtype=float)
                    self.__vec_c = np.array([0.0, 0.0, z_length], dtype=float)

                else:
                    lineno = 1
                    line = pygetline(filename=self.fort10, lineno=lineno).split()
                    self.__vec_a_1 = Value(
                        value=float(line[0]),
                        lineno=lineno,
                        index=0,
                        file=self.fort10,
                    )
                    self.__vec_a_2 = Value(
                        value=float(line[1]),
                        lineno=lineno,
                        index=1,
                        file=self.fort10,
                    )
                    self.__vec_a_3 = Value(
                        value=float(line[2]),
                        lineno=lineno,
                        index=2,
                        file=self.fort10,
                    )
                    self.__vec_b_1 = Value(
                        value=float(line[3]),
                        lineno=lineno,
                        index=3,
                        file=self.fort10,
                    )
                    self.__vec_b_2 = Value(
                        value=float(line[4]),
                        lineno=lineno,
                        index=4,
                        file=self.fort10,
                    )
                    self.__vec_b_3 = Value(
                        value=float(line[5]),
                        lineno=lineno,
                        index=5,
                        file=self.fort10,
                    )
                    self.__vec_c_1 = Value(
                        value=float(line[6]),
                        lineno=lineno,
                        index=6,
                        file=self.fort10,
                    )
                    self.__vec_c_2 = Value(
                        value=float(line[7]),
                        lineno=lineno,
                        index=7,
                        file=self.fort10,
                    )
                    self.__vec_c_3 = Value(
                        value=float(line[8]),
                        lineno=lineno,
                        index=8,
                        file=self.fort10,
                    )
                    self.__phase_up_1 = Value(
                        value=float(line[9]),
                        lineno=lineno,
                        index=9,
                        file=self.fort10,
                    )  # phase 1 up
                    self.__phase_up_2 = Value(
                        value=float(line[10]),
                        lineno=lineno,
                        index=10,
                        file=self.fort10,
                    )  # phase 2 up
                    self.__phase_up_3 = Value(
                        value=float(line[11]),
                        lineno=lineno,
                        index=11,
                        file=self.fort10,
                    )  # phase 3 up
                    self.__phase_dn_1 = Value(
                        value=float(line[12]),
                        lineno=lineno,
                        index=12,
                        file=self.fort10,
                    )  # phase 1 dn
                    self.__phase_dn_2 = Value(
                        value=float(line[13]),
                        lineno=lineno,
                        index=13,
                        file=self.fort10,
                    )  # phase 2 dn
                    self.__phase_dn_3 = Value(
                        value=float(line[14]),
                        lineno=lineno,
                        index=14,
                        file=self.fort10,
                    )  # phase 3 dn

                    self.__vec_a = np.array(
                        [self.__vec_a_1.v, self.__vec_a_2.v, self.__vec_a_3.v]
                    )
                    self.__vec_b = np.array(
                        [self.__vec_b_1.v, self.__vec_b_2.v, self.__vec_b_3.v]
                    )
                    self.__vec_c = np.array(
                        [self.__vec_c_1.v, self.__vec_c_2.v, self.__vec_c_3.v]
                    )

            else:
                self.__vec_a = np.array([0.0, 0.0, 0.0], dtype=float)
                self.__vec_b = np.array([0.0, 0.0, 0.0], dtype=float)
                self.__vec_c = np.array([0.0, 0.0, 0.0], dtype=float)

            # read structures
            self.__atomic_numbers = []
            self.__valence_electrons = []
            self.__positions = []
            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be false (for loop), otherwise so slow

            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            for i in range(self.ion):
                p[5 * i + 0].type(float)
                p[5 * i + 1].type(float)
                p[5 * i + 2].type(float)
                p[5 * i + 3].type(float)
                p[5 * i + 4].type(float)
                self.__valence_electrons.append(p[5 * i + 0])
                self.__atomic_numbers.append(p[5 * i + 1])
                self.__positions.append([p[5 * i + 2], p[5 * i + 3], p[5 * i + 4]])

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def structure(self):
        self.read()
        __cell = Cell(vec_a=self.vec_a, vec_b=self.vec_b, vec_c=self.vec_c)
        __structure = Structure(
            cell=__cell,
            atomic_numbers=self.atomic_numbers,
            element_symbols=[return_element_symbol(z) for z in self.atomic_numbers],
            positions=self.positions,
        )
        return __structure

    @property
    def ortho_flag(self):
        self.read()
        return not self.__tilted_flag

    @property
    def pbc_flag(self):
        self.read()
        return self.__pbc

    @property
    def vec_a(self):
        self.read()
        return self.__vec_a

    @property
    def vec_b(self):
        self.read()
        return self.__vec_b

    @property
    def vec_c(self):
        self.read()
        return self.__vec_c

    @property
    def norm_vec_a(self):
        self.read()
        return np.linalg.norm(self.__vec_a)

    @property
    def norm_vec_b(self):
        self.read()
        return np.linalg.norm(self.__vec_b)

    @property
    def norm_vec_c(self):
        self.read()
        return np.linalg.norm(self.__vec_c)

    @property
    def atomic_numbers(self):
        # logger.debug("atomic numbers")
        self.read()
        return [i.v for i in self.__atomic_numbers]

    @property
    def valence_electrons(self):
        # logger.debug("valence electrons")
        self.read()
        return [i.v for i in self.__valence_electrons]

    @property
    def positions(self):
        self.read()
        return np.array(
            [
                [positions[0].v, positions[1].v, positions[2].v]
                for positions in self.__positions
            ]
        )

    @property
    def phase_up(self):
        self.read()
        return [self.__phase_up_1.v, self.__phase_up_2.v, self.__phase_up_3.v]

    @property
    def phase_dn(self):
        self.read()
        return [self.__phase_dn_1.v, self.__phase_dn_2.v, self.__phase_dn_3.v]

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10forceconstraint:
    def __init__(
        self,
        fort10: str,
        ieskinr: int,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.ieskinr = ieskinr
        self.in_place = in_place
        self.read_flag = False

        self.__constraint_num = []
        self.__constraint_index = []
        self.__atom_label = []
        self.__direction = []

    def read(self):
        if not self.read_flag:
            # read constraints
            self.__constraint_num = []
            self.__constraint_index = []
            self.__atom_label = []
            self.__direction = []

            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be false (for for loop); otherwise very slow
            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            c_index = 0
            for i in range(self.ieskinr):
                p[c_index].type(int)
                self.__constraint_num.append(p[c_index])
                num_constraints = p[c_index].v

                for num in range(np.abs(num_constraints)):
                    p[c_index + 2 * num + 1].type(int)
                    p[c_index + 2 * num + 2].type(int)
                    self.__constraint_index.append(i)
                    self.__atom_label.append(p[c_index + 2 * num + 1])
                    self.__direction.append(p[c_index + 2 * num + 2])

                c_index += 2 * num_constraints + 1

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def constraint_num(self):
        self.read()
        return [i.v for i in self.__constraint_num]

    @property
    def constraint_index(self):
        self.read()
        return [i for i in self.__constraint_index]

    @property
    def atom_label(self):
        self.read()
        return [i.v for i in self.__atom_label]

    @property
    def direction(self):
        self.read()
        return [i.v for i in self.__direction]

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10jastwobody:
    def __init__(
        self,
        fort10: str,
        jas_2body: int,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.read_flag = False
        self.jastrow_type = jas_2body
        self.in_place = in_place

        # Values()
        self.__num_jas_param = 0
        self.__twobody_list = []
        self.__onebody_list = []

    def read(self):
        if not self.read_flag:
            # read det basis sets, molecular orbitals, and hybrid orbitals
            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be false
            # a manual clearcache is better to avoid unexpected
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            p[0].type(int)
            self.__num_jas_param = p[0]
            num_jas_param = np.abs(p[0].v)
            num_twobody, flag_onebody = return_num_twobody_and_flag_onebody(
                jastrow_type=self.jastrow_type
            )
            # logger.debug(self.jastrow_type)
            # logger.debug(num_twobody)
            # logger.debug(flag_onebody)
            if num_twobody == 0:
                # why here the explicit "if" is needed?
                # because self.__num_param is not always 0!
                # -> jastype=0 with symagp=.false.
                self.__twobody_list = []
                self.__onebody_list = []

            elif flag_onebody:  # yes
                num_twobody_jas = num_twobody
                num_onebody_jas = num_jas_param - num_twobody
                for i in range(num_twobody_jas):
                    p[1 + i].type(float)
                    self.__twobody_list.append(p[1 + i])
                for i in range(num_onebody_jas):
                    p[1 + num_twobody_jas + i].type(float)
                    self.__onebody_list.append(p[1 + num_twobody_jas + i])

            else:
                num_twobody_jas = num_twobody
                for i in range(num_twobody_jas):
                    p[1 + i].type(float)
                    self.__twobody_list.append(p[1 + i])
                self.__onebody_list = []

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def num_jas_param(self):
        self.read()
        return self.__num_jas_param.v

    @property
    def twobody_list(self):
        self.read()
        return [i.v for i in self.__twobody_list]

    @property
    def onebody_list(self):
        self.read()
        return [i.v for i in self.__onebody_list]

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10detbasissets:
    def __init__(
        self,
        fort10: str,
        shell_det: int,
        complex_flag: bool,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.read_flag = False
        self.shell_det = shell_det
        self.complex_flag = complex_flag
        self.in_place = in_place

        # Values()
        self.__atom_label = []  # starting from 1!
        self.__param_num = []
        self.__shell_multiplicity = []
        self.__shell_ang_mom_turbo_notation = []
        self.__exponent = []
        self.__coefficient = []
        self.__coefficient_imag = []
        # auxiliary values
        self.__shell_index = []

        # Values()
        self.__mo_atom_label = []  # starting from 1!
        self.__mo_param_num = []
        self.__mo_shell_multiplicity = []
        self.__mo_shell_ang_mom_turbo_notation = []
        self.__mo_prim_index = []
        # dimension (prim_ao.num, mo.num) prim_ao.num != ao.num
        # in the turborvb implementation.
        self.__mo_coefficient = []
        # dimension (prim_ao.num, mo.num) prim_ao.num != ao.num
        # in the turborvb implementation.
        self.__mo_coefficient_imag = []
        # dimension (prim_ao.num, mo.num) prim_ao.num != ao.num
        # in the turborvb implementation.

        # Values()
        self.__hyb_atom_label = []  # starting from 1!
        self.__hyb_param_num = []
        self.__hyb_shell_multiplicity = []
        self.__hyb_shell_ang_mom_turbo_notation = []
        self.__hyb_prim_index = []
        # dimension (prim_ao.num, hyb.num) prim_ao.num != ao.num
        # in the turborvb implementation.
        self.__hyb_coefficient = []
        # dimension (prim_ao.num, hyb.num) prim_ao.num != ao.num
        # in the turborvb implementation.
        self.__hyb_coefficient_imag = []
        # dimension (prim_ao.num, hyb.num) prim_ao.num != ao.num
        # in the turborvb implementation.

    def read(self):
        if not self.read_flag:
            # logger.debug("Reading!!")
            # read det basis sets, molecular orbitals, and hybrid orbitals
            # logger.debug("Check0")
            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be False otherwise very slow.
            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            # logger.debug("Check1")
            for line in lines:
                # logger.debug(line)
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            # logger.debug("Check2")
            for num in range(abs(self.shell_det)):
                p[0].type(int)
                multiplicity = p[0]
                p[1].type(int)
                param_num = p[1]
                p[2].type(int)
                shell_ang_mom_turbo_notation = p[2]
                p[3].type(int)
                atom_label = p[3]

                orb_type_chr = str(return_orb_type_chr(shell_ang_mom_turbo_notation.v))

                if orb_type_chr == "hyb":
                    self.__hyb_shell_multiplicity.append(multiplicity)
                    self.__hyb_param_num.append(param_num)
                    self.__hyb_shell_ang_mom_turbo_notation.append(
                        shell_ang_mom_turbo_notation
                    )
                    self.__hyb_atom_label.append(atom_label)

                    hyb_prim_index = []
                    hyb_coefficient = []
                    hyb_coefficient_imag = []
                    num = int(param_num.v / 2)

                    if not self.complex_flag:  # real
                        for n in p[4 : 4 + num]:
                            n.type(int)
                            hyb_prim_index.append(n)
                        for n in p[4 + num : 4 + param_num.v]:
                            n.type(float)
                            hyb_coefficient.append(n)

                    else:  # complex
                        for n in p[4 : 4 + num]:
                            n.type(int)
                            hyb_prim_index.append(n)
                        for i, n in enumerate(
                            p[4 + num : 4 + int(param_num.v * 3.0 / 2.0)]
                        ):
                            if i % 2 == 0:
                                n.type(float)
                                hyb_coefficient.append(n)
                            else:
                                n.type(float)
                                hyb_coefficient_imag.append(n)

                    self.__hyb_prim_index.append(hyb_prim_index)
                    self.__hyb_coefficient.append(hyb_coefficient)
                    self.__hyb_coefficient_imag.append(hyb_coefficient_imag)

                elif orb_type_chr == "mol":
                    self.__mo_shell_multiplicity.append(multiplicity)
                    self.__mo_param_num.append(param_num)
                    self.__mo_shell_ang_mom_turbo_notation.append(
                        shell_ang_mom_turbo_notation
                    )
                    self.__mo_atom_label.append(atom_label)

                    mo_prim_index = []
                    mo_coefficient = []
                    mo_coefficient_imag = []
                    ind_num = int(param_num.v / 2)

                    if not self.complex_flag:  # real case
                        for n in p[4 : 4 + ind_num]:
                            n.type(int)
                            mo_prim_index.append(n)
                        for n in p[4 + ind_num : 4 + param_num.v]:
                            n.type(float)
                            mo_coefficient.append(n)

                    else:  # complex case!!
                        for n in p[4 : 4 + ind_num]:
                            n.type(int)
                            mo_prim_index.append(n)
                        for i, n in enumerate(
                            p[4 + ind_num : 4 + int(param_num.v * 3.0 / 2.0)]
                        ):
                            if i % 2 == 0:
                                n.type(float)
                                mo_coefficient.append(n)
                            else:
                                n.type(float)
                                mo_coefficient_imag.append(n)

                    self.__mo_prim_index.append(mo_prim_index)
                    self.__mo_coefficient.append(mo_coefficient)
                    self.__mo_coefficient_imag.append(mo_coefficient_imag)

                elif orb_type_chr in {"s", "p", "d", "f", "g", "h", "i"}:
                    shell_index = len(self.__shell_multiplicity)
                    self.__shell_multiplicity.append(multiplicity)
                    self.__param_num.append(param_num)
                    self.__shell_ang_mom_turbo_notation.append(
                        shell_ang_mom_turbo_notation
                    )
                    self.__atom_label.append(atom_label)

                    # logger.debug(f"multi= {multiplicity.v}")
                    # logger.debug(f"param_num = {param_num.v}")
                    # logger.debug(
                    #    f"orb_type_num = {shell_ang_mom_turbo_notation.v}"
                    # )
                    # logger.debug(f"atom_label = {atom_label.v}")
                    contraction_flag = return_contraction_flag(
                        shell_ang_mom_turbo_notation.v
                    )

                    if contraction_flag:
                        num = int(param_num.v / 2)

                        if not self.complex_flag:  # real case
                            for n in p[4 : 4 + num]:
                                n.type(float)
                                self.__exponent.append(n)
                                self.__shell_index.append(shell_index)
                            for n in p[4 + num : 4 + param_num.v]:
                                n.type(float)
                                self.__coefficient.append(n)
                        else:  # complex case
                            for n in p[4 : 4 + num]:
                                n.type(float)
                                self.__exponent.append(n)
                                self.__shell_index.append(shell_index)
                            for i, n in enumerate(
                                p[4 + num : 4 + int(param_num.v * 3.0 / 2.0)]
                            ):
                                if i % 2 == 0:
                                    n.type(float)
                                    self.__coefficient.append(n)
                                else:
                                    n.type(float)
                                    self.__coefficient_imag.append(n)
                    else:
                        p[4].type(float)
                        self.__exponent.append(p[4])
                        self.__shell_index.append(shell_index)
                        self.__coefficient.append(Value())

                if not self.complex_flag:  # real case
                    p = p[4 + param_num.v :]
                else:  # complex case
                    p = p[4 + int(param_num.v * 3.0 / 2.0) :]

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def det_basis_sets(self):
        self.read()
        return Det_Basis_sets(
            # det basis
            nucleus_index=[a.v - 1 for a in self.__atom_label],
            shell_ang_mom=[int((a.v - 1) / 2) for a in self.__shell_multiplicity],
            shell_ang_mom_turbo_notation=[
                a.v for a in self.__shell_ang_mom_turbo_notation
            ],
            shell_factor=[1.0] * len(self.__atom_label),
            shell_index=[a for a in self.__shell_index],
            exponent=[a.v for a in self.__exponent],
            coefficient=[a.v if a.v is not None else 1.0 for a in self.__coefficient],
            coefficient_imag=[a.v for a in self.__coefficient_imag],
            prim_factor=[1.0] * len(self.__exponent),
            # hybrid
            hyb_nucleus_index=[a.v - 1 for a in self.__hyb_atom_label],
            hyb_param_num=[a.v for a in self.__hyb_param_num],
            hyb_shell_ang_mom=[0 for a in self.__hyb_shell_multiplicity],
            hyb_shell_ang_mom_turbo_notation=[
                a.v for a in self.__hyb_shell_ang_mom_turbo_notation
            ],
            hyb_prim_index=[
                [a.v - 1 for a in hyb_prim_list]
                for hyb_prim_list in self.__hyb_prim_index
            ],
            hyb_coefficient=[
                [a.v for a in hyb_coeff_list]
                for hyb_coeff_list in self.__hyb_coefficient
            ],
            hyb_coefficient_imag=[
                [a.v for a in hyb_coeff_imag_list]
                for hyb_coeff_imag_list in self.__hyb_coefficient_imag
            ],
        )

    @property
    def shell_index(self):
        self.read()
        return self.__shell_index
        # self.read(); return [[i.v for i in j] for j in self.__shell_index]

    @property
    def has_mo(self):
        self.read()
        if len(self.__mo_atom_label) > 0:
            return True
        else:
            return False

    @property
    def num_mo(self):
        self.read()
        return len(self.__mo_atom_label)

    @property
    def mo_coefficient(self):
        self.read()
        return [[i.v for i in j] for j in self.__mo_coefficient]

    @property
    def mo_coefficient_imag(self):
        self.read()
        return [[i.v for i in j] for j in self.__mo_coefficient_imag]

    @property
    def hyb_coefficient(self):
        self.read()
        return [[i.v for i in j] for j in self.__hyb_coefficient]

    @property
    def hyb_coefficient_imag(self):
        self.read()
        return [[i.v for i in j] for j in self.__hyb_coefficient_imag]

    @mo_coefficient.setter
    def mo_coefficient(self, new_mo_coefficient):
        self.read()
        # logger.debug(len(new_mo_coefficient))
        # logger.debug(len(self.__mo_coefficient))
        # logger.debug(new_mo_coefficient[0][-1])
        assert len(new_mo_coefficient) == len(self.__mo_coefficient)
        total_num_sed = len(new_mo_coefficient) * len(new_mo_coefficient[0])
        logger.debug(f"Total num sed = {total_num_sed}")
        if self.in_place:
            """with gnu-sed (very slow after improvement!!)
            logger.debug("replace by gnu-sed")
            self.replace_mo_coeff_with_sed(
                self.__mo_coefficient, new_mo_coefficient
            )
            """
            # """ with readlines
            logger.debug("replace by sed-python (real part)")
            self.replace_mo_coeff_pure_python(self.__mo_coefficient, new_mo_coefficient)
            # """
        # todo, self.__mo_coefficient itself should be replaced
        # with new_mo_coefficient!! at present, fort.10 is not updated.

    @mo_coefficient_imag.setter
    def mo_coefficient_imag(self, new_mo_coefficient_imag):
        self.read()
        # logger.debug(len(new_mo_coefficient))
        # logger.debug(len(self.__mo_coefficient_imag))
        # logger.debug(new_mo_coefficient[0][-1])
        assert len(new_mo_coefficient_imag) == len(self.__mo_coefficient_imag)
        total_num_sed = len(new_mo_coefficient_imag) * len(new_mo_coefficient_imag[0])
        logger.debug(f"Total num sed = {total_num_sed}")
        if self.in_place:
            """with gnu-sed (very slow after improvement!!)
            logger.debug("replace by gnu-sed")
            self.replace_mo_coeff_with_sed(
                self.__mo_coefficient_imag, new_mo_coefficient_imag
            )
            """
            # """ with readlines
            logger.debug("replace by sed-python (img. part)")
            self.replace_mo_coeff_pure_python(
                self.__mo_coefficient_imag, new_mo_coefficient_imag
            )
            # """
        # todo, self.__mo_coefficient_imag itself should be replaced
        # with new_mo_coefficient!! at present, fort.10 is not updated.

    def replace_mo_coeff_with_sed(self, old_mo_coefficient, new_mo_coefficient):
        start_sed = time.time()
        file_list = []
        lineno_list = []
        value_list = []
        index_list = []
        for mo_coeff, new_mo_coeff in zip(old_mo_coefficient, new_mo_coefficient):
            assert len(mo_coeff) == len(new_mo_coeff)
            lineno = -1
            files = []
            index_same_lineno = []
            new_coeff_same_lineno = []
            # logger.debug(mo_coeff)
            # logger.debug(new_mo_coeff)
            for coeff, new_coeff in zip(mo_coeff, new_mo_coeff):
                l = coeff.l
                # v = coeff.v
                i = coeff.i
                f = coeff.f

                if lineno == -1:
                    lineno = l
                    files.append(f)
                    index_same_lineno.append(i)
                    new_coeff_same_lineno.append(new_coeff)

                else:
                    if lineno == l:
                        files.append(f)
                        index_same_lineno.append(i)
                        new_coeff_same_lineno.append(new_coeff)

                    else:
                        assert len(set(files)) == 1
                        file = files[0]

                        file_list.append(file)
                        lineno_list.append(lineno)
                        value_list.append(new_coeff_same_lineno)
                        index_list.append(index_same_lineno)

                        files = []
                        index_same_lineno = []
                        new_coeff_same_lineno = []

                        lineno = l
                        files.append(f)
                        index_same_lineno.append(i)
                        new_coeff_same_lineno.append(new_coeff)

            if len(files) != 0:
                assert len(set(files)) == 1
                file = files[0]
                file_list.append(file)
                lineno_list.append(lineno)
                value_list.append(new_coeff_same_lineno)
                index_list.append(index_same_lineno)

        assert len(set(file_list)) == 1
        file = file_list[0]
        pysed_replace_lines(
            file=file,
            lineno_list=lineno_list,
            value_list=value_list,
            index_list=index_list,
        )
        end_sed = time.time()
        logger.info("elapsed time for sed:{:f}".format(end_sed - start_sed) + "[sec]")

    def replace_mo_coeff_pure_python(self, old_mo_coefficient, new_mo_coefficient):
        line_no_list = []
        index_list = []
        w_mo_coeff_list = []
        for old_mo_coeff, new_mo_coeff in zip(old_mo_coefficient, new_mo_coefficient):
            if len(old_mo_coeff) != len(new_mo_coeff):
                logger.error(
                    f"len(old_mo_coeff):{len(old_mo_coeff)} != len(new_mo_coeff):{len(new_mo_coeff)}"
                )
                raise ValueError
            line_no_list += [coeff.l for coeff in old_mo_coeff]
            index_list += [coeff.i for coeff in old_mo_coeff]
            w_mo_coeff_list += new_mo_coeff

        available_memory = psutil.virtual_memory().available  # Byte
        size_of_fort10 = os.path.getsize(self.fort10)  # Byte
        logger.debug(f"available memory = {available_memory >> 20} MiB")
        logger.debug(f"fort10 size = {size_of_fort10 >> 20} MiB")

        if available_memory * 0.90 > size_of_fort10:
            logger.debug("available memory > fort10 size")
            # version with readlines (more memory, but faster)
            with open(self.fort10, "r") as f:
                lines = f.readlines()

            # """ straightforward, but established way
            tqdm_out = TqdmToLogger(logger, level=logging.INFO)

            # new way!! a big loop (line_no_list) is iterated only once!!
            line_no_index_list_dict = {line_no: [] for line_no in line_no_list}
            for i, line_no in enumerate(line_no_list):
                line_no_index_list_dict[line_no].append(i)

            for line_no in tqdm(
                list(set(line_no_list)),
                file=tqdm_out,
                miniters=int(len(list(set(line_no_list))) / 10),
                maxinterval=1.0e5,
            ):
                # old way, very slow!! because
                # a big loop (line_no_list) is iterated twice!!
                """
                line_no_index_list = [
                    i for i, x in enumerate(line_no_list) if x == line_no
                ]
                """
                # new way
                line_no_index_list = line_no_index_list_dict[line_no]

                r_index_list = [index_list[i] for i in line_no_index_list]

                r_new_mo_coeff_list = [w_mo_coeff_list[i] for i in line_no_index_list]

                if len(r_index_list) != len(set(r_index_list)):
                    logger.error("Duplicated r_index!! It should not happen.")
                    raise ValueError
                line = lines[line_no].split()

                # old way (explicit for loop)
                """
                for r_index, r_new_mo_coeff in zip(
                    r_index_list, r_new_mo_coeff_list
                ):
                    line[r_index] = r_new_mo_coeff
                """
                # new way (implicit for loop)
                line = [
                    r_new_mo_coeff_list[r_index_list.index(i)]
                    if i in r_index_list
                    else l
                    for i, l in enumerate(line)
                ]
                lines[line_no] = " ".join(list(map(str, line))) + "\n"

            # """

            with open(self.fort10, "w") as f:
                f.writelines(lines)

        else:
            logger.debug("available memory < fort10 size")

            # version with readline (less memory, but super slow!!)
            sed_fort10 = self.fort10 + "_sed_tmp"
            with open(self.fort10, "r") as f:
                with open(sed_fort10, "w") as fw:
                    line = f.readline()
                    cnt = 0
                    while line:
                        if cnt % 1000 == 0:
                            logger.debug(f"cnt={cnt}")
                        line_no_index_list = [
                            i for i, x in enumerate(line_no_list) if x == cnt
                        ]
                        r_index_list = [index_list[i] for i in line_no_index_list]
                        r_new_mo_coeff_list = [
                            w_mo_coeff_list[i] for i in line_no_index_list
                        ]
                        line = line.split()
                        for r_index, r_new_mo_coeff in zip(
                            r_index_list, r_new_mo_coeff_list
                        ):
                            line[r_index] = r_new_mo_coeff
                        line = list(map(str, line))
                        line = " ".join(line) + "\n"
                        fw.write(line)
                        line = f.readline()
                        cnt += 1
            os.remove(self.fort10)
            shutil.copy(sed_fort10, self.fort10)
            os.rename(sed_fort10, self.fort10)

    @property
    def coefficient(self):
        self.read()
        return [i.v for i in self.__coefficient]

    @property
    def coefficient_imag(self):
        self.read()
        return [i.v for i in self.__coefficient_imag]

    @property
    def exponent(self):
        self.read()
        return [i.v for i in self.__exponent]

    @property
    def shell_multiplicity(self):
        self.read()
        return [i.v for i in self.__shell_multiplicity]

    @property
    def hyb_atom_label(self):
        self.read()
        return [i.v for i in self.__hyb_atom_label]

    @property
    def hyb_param_num(self):
        self.read()
        return [i.v for i in self.__hyb_param_num]

    @property
    def hyb_shell_multiplicity(self):
        self.read()
        return [i.v for i in self.__hyb_shell_multiplicity]

    @property
    def hyb_shell_ang_mom_turbo_notation(self):
        self.read()
        return [i.v for i in self.__hyb_shell_ang_mom_turbo_notation]

    @property
    def hyb_prim_index(self):
        self.read()
        return [i.v for i in self.__hyb_prim_index]

    """ to be deleted
    @property
    def hyb_coefficient(self):
        self.read()
        return [i.v for i in self.__hyb_coefficient]

    @property
    def hyb_coefficient_imag(self):
        self.read()
        return [i.v for i in self.__hyb_coefficient_imag]
    """

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10jasbasissets:
    def __init__(
        self,
        fort10: str,
        shell_jas: int,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.read_flag = False
        self.shell_jas = shell_jas
        self.in_place = in_place

        # Values()
        self.__atom_label = []  # starting from 1!
        self.__param_num = []
        self.__shell_multiplicity = []
        self.__shell_ang_mom_turbo_notation = []
        self.__exponent = []
        self.__coefficient = []

        # auxiliary values
        self.__shell_index = []

    def read(self):
        if not self.read_flag:
            # read jas basis sets, molecular orbitals, and hybrid orbitals
            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be False, otherwise very slow
            # a manual clearcache is better to avoid unexpected behavior
            linecache.clearcache()
            # print(lines)
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            for num in range(abs(self.shell_jas)):
                p[0].type(int)
                multiplicity = p[0]
                p[1].type(int)
                param_num = p[1]
                p[2].type(int)
                shell_ang_mom_turbo_notation = p[2]
                p[3].type(int)
                atom_label = p[3]

                shell_index = len(self.__shell_multiplicity)
                self.__shell_multiplicity.append(multiplicity)
                self.__param_num.append(param_num)
                self.__shell_ang_mom_turbo_notation.append(shell_ang_mom_turbo_notation)
                self.__atom_label.append(atom_label)

                # logger.debug(f"multi= {multiplicity.v}")
                # logger.debug(f"param_num = {param_num.v}")
                # logger.debug(
                #    f"orb_type_num = {shell_ang_mom_turbo_notation.v}"
                # )
                # logger.debug(f"atom_label = {atom_label.v}")

                orb_type_chr = str(return_orb_type_chr(shell_ang_mom_turbo_notation.v))

                if orb_type_chr in {"s", "p", "d", "f", "g", "h", "i"}:
                    contraction_flag = return_contraction_flag(
                        shell_ang_mom_turbo_notation.v
                    )

                    if contraction_flag:
                        num = int(param_num.v / 2)
                        for n in p[4 : 4 + num]:
                            n.type(float)
                            self.__exponent.append(n)
                            self.__shell_index.append(shell_index)
                        for n in p[4 + num : 4 + param_num.v]:
                            n.type(float)
                            self.__coefficient.append(n)

                    else:
                        p[4].type(float)
                        self.__exponent.append(p[4])
                        self.__shell_index.append(shell_index)
                        self.__coefficient.append(Value())

                elif orb_type_chr in {"jas_const"}:
                    self.__shell_index.append(shell_index)
                    self.__exponent.append(Value())
                    self.__coefficient.append(Value())

                p = p[4 + param_num.v :]

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def jas_basis_sets(self):
        self.read()
        return Jas_Basis_sets(
            nucleus_index=[a.v - 1 for a in self.__atom_label[:-1]],
            shell_ang_mom=[int((a.v - 1) / 2) for a in self.__shell_multiplicity[:-1]],
            shell_ang_mom_turbo_notation=[
                a.v for a in self.__shell_ang_mom_turbo_notation[:-1]
            ],
            shell_factor=[1.0] * len(self.__atom_label[:-1]),
            shell_index=[a for a in self.__shell_index[:-1]],
            exponent=[a.v for a in self.__exponent[:-1]],
            coefficient=[
                a.v if a.v is not None else 1.0 for a in self.__coefficient[:-1]
            ],
            prim_factor=[1.0] * len(self.__atom_label[:-1]),
        )

    @property
    def coefficient(self):
        self.read()
        return [i.v for i in self.__coefficient]

    @property
    def exponent(self):
        self.read()
        return [i.v for i in self.__exponent]

    @property
    def shell_multiplicity(self):
        self.read()
        return [i.v for i in self.__shell_multiplicity]

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10occ:
    def __init__(
        self,
        fort10: str,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        # logger.debug("F10occ starts")
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.in_place = in_place
        self.read_flag = False
        # logger.debug("F10occ ends")

        # Values()
        self.__occupation = []

    def read(self):
        if not self.read_flag:
            # read det basis sets, molecular orbitals, and hybrid orbitals
            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be false, otherwise very slow.
            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            # print(lines)
            # assert len(lines) == abs(np.sum(self.shell_multiplicity))
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            for num in range(len(lines)):
                p[num].type(int)
                occ = p[num]
                self.__occupation.append(occ)

                logger.debug(f"multi= {occ.v}")

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def occ(self):
        self.read()
        return [i.v for i in self.__occupation]

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10detmat:
    def __init__(
        self,
        fort10: str,
        detmat: int,
        complex_flag: bool,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.in_place = in_place
        self.read_flag = False
        self.detmat = detmat
        self.complex_flag = complex_flag

        # Values()
        self.__row = []
        self.__col = []
        self.__coeff_real = []
        self.__coeff_imag = []  # to be implemented.

    def read(self):
        if not self.read_flag:
            # read det basis sets, molecular orbitals, and hybrid orbitals
            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be False otherwise very slow
            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            c_index = 0
            for num in range(abs(self.detmat)):
                p[c_index + 0].type(int)
                row = p[c_index + 0]
                p[c_index + 1].type(int)
                col = p[c_index + 1]
                p[c_index + 2].type(float)
                coeff_real = p[c_index + 2]
                if self.complex_flag:
                    p[c_index + 3].type(float)
                    coeff_imag = p[c_index + 3]

                if not self.complex_flag:
                    self.__row.append(row)
                    self.__col.append(col)
                    self.__coeff_real.append(coeff_real)
                    c_index += 3
                else:  # complex
                    self.__row.append(row)
                    self.__col.append(col)
                    self.__coeff_real.append(coeff_real)
                    self.__coeff_imag.append(coeff_imag)
                    c_index += 4

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def row(self):
        self.read()
        return [i.v for i in self.__row]

    @row.setter
    def row(self, value_list):
        assert len(self.__row) == len(value_list)
        self.read()
        for row, value in zip(self.__row, value_list):
            if row.v != value:
                row.replace(value=value, in_place=self.in_place)

    @property
    def col(self):
        self.read()
        return [i.v for i in self.__col]

    @col.setter
    def col(self, value_list):
        assert len(self.__col) == len(value_list)
        self.read()
        for col, value in zip(self.__col, value_list):
            if col.v != value:
                col.replace(value=value, in_place=self.in_place)

    @property
    def coeff_real(self):
        self.read()
        return [i.v for i in self.__coeff_real]

    @coeff_real.setter
    def coeff_real(self, value_list):
        assert len(self.__coeff_real) == len(value_list)
        self.read()
        for coeff_real, value in zip(self.__coeff_real, value_list):
            if coeff_real.v != value:
                coeff_real.replace(value=value, in_place=self.in_place)

    @property
    def coeff_imag(self):
        self.read()
        return [i.v for i in self.__coeff_imag]

    @coeff_imag.setter
    def coeff_imag(self, value_list):
        assert len(self.__coeff_imag) == len(value_list)
        self.read()
        for coeff_imag, value in zip(self.__coeff_imag, value_list):
            if coeff_imag.v != value:
                coeff_imag.replace(value=value, in_place=self.in_place)

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10jasmat:
    def __init__(
        self,
        fort10: str,
        jasmat: int,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.in_place = in_place
        self.read_flag = False
        self.jasmat = jasmat

        # Values()
        self.__row = []
        self.__col = []
        self.__coeff = []

    def read(self):
        if not self.read_flag:
            # read det basis sets, molecular orbitals, and hybrid orbitals
            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be false, otherwise very slow.
            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            c_index = 0
            for num in range(abs(self.jasmat)):
                p[c_index + 0].type(int)
                row = p[c_index + 0]
                p[c_index + 1].type(int)
                col = p[c_index + 1]
                p[c_index + 2].type(float)
                coeff_real = p[c_index + 2]

                self.__row.append(row)
                self.__col.append(col)
                self.__coeff.append(coeff_real)

                c_index += 3

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def row(self):
        self.read()
        return [i.v for i in self.__row]

    @row.setter
    def row(self, value_list):
        assert len(self.__row) == len(value_list)
        self.read()
        for row, value in zip(self.__row, value_list):
            if row.v != value:
                row.replace(value=value, in_place=self.in_place)

    @property
    def col(self):
        self.read()
        return [i.v for i in self.__col]

    @col.setter
    def col(self, value_list):
        assert len(self.__col) == len(value_list)
        self.read()
        for col, value in zip(self.__col, value_list):
            if col.v != value:
                col.replace(value=value, in_place=self.in_place)

    @property
    def coeff(self):
        self.read()
        return [i.v for i in self.__coeff]

    @coeff.setter
    def coeff(self, value_list):
        assert len(self.__coeff) == len(value_list)
        self.read()
        for coeff, value in zip(self.__coeff, value_list):
            if coeff.v != value:
                coeff.replace(value=value, in_place=self.in_place)

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10matsymmetry:
    def __init__(
        self,
        fort10: str,
        num_component: int,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.num_component = num_component
        self.in_place = in_place
        self.read_flag = False

        self.__constraint_num = []
        self.__constraint_index = []
        self.__row = []
        self.__column = []

    def read(self):
        if not self.read_flag:
            # read constraints
            self.__constraint_num = []
            self.__constraint_index = []
            self.__row = []
            self.__column = []

            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clearcache should be false otherwise very slow.
            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            c_index = 0
            for i in range(self.num_component):
                p[c_index].type(int)
                self.__constraint_num.append(p[c_index])
                num_constraints = np.abs(p[c_index].v)
                logger.debug(num_constraints)
                for num in range(num_constraints):
                    p[c_index + 2 * num + 1].type(int)
                    p[c_index + 2 * num + 2].type(int)
                    logger.debug(p[c_index + 2 * num + 1].v)
                    logger.debug(p[c_index + 2 * num + 2].v)
                    self.__constraint_index.append(i)
                    self.__row.append(p[c_index + 2 * num + 1])
                    self.__column.append(p[c_index + 2 * num + 2])

                c_index += 2 * num_constraints + 1

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def constraint_num(self):
        self.read()
        return [i.v for i in self.__constraint_num]

    @constraint_num.setter
    def constraint_num(self, value_list):
        assert len(self.__constraint_num) == len(value_list)
        self.read()
        for constraint, value in zip(self.__constraint_num, value_list):
            if constraint.v != value:
                constraint.replace(value=value, in_place=self.in_place)

    @property
    def constraint_index(self):
        self.read()
        return [i for i in self.__constraint_index]

    @property
    def row(self):
        self.read()
        return [i.v for i in self.__row]

    @row.setter
    def row(self, value_list):
        assert len(self.__row) == len(value_list)
        self.read()
        for row, value in zip(self.__row, value_list):
            if row.v != value:
                row.replace(value=value, in_place=self.in_place)

    @property
    def column(self):
        self.read()
        return [i.v for i in self.__column]

    @column.setter
    def column(self, value_list):
        assert len(self.__column) == len(value_list)
        self.read()
        for column, value in zip(self.__column, value_list):
            if column.v != value:
                column.replace(value=value, in_place=self.in_place)

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


class F10basissymmetry:
    def __init__(
        self,
        fort10: str,
        num_component: int,
        start_keyword: str,
        end_keyword: str,
        in_place: bool = True,
    ):
        self.start_keyword = start_keyword
        self.end_keyword = end_keyword
        self.fort10 = fort10
        self.num_component = num_component
        self.in_place = in_place
        self.read_flag = False

        self.__constraint_num = []
        self.__constraint_index = []
        self.__basis_index = []

    def read(self):
        if not self.read_flag:
            # read constraints
            self.__constraint_num = []
            self.__constraint_index = []
            self.__basis_index = []

            lines = [
                pygetline(filename=self.fort10, lineno=lineno, clearcache=False)
                for lineno in range(self.start_lineno, self.end_lineno + 1)
            ]  # clear cache should be false otherwise very slow.
            # a manual clearcache is better to avoid
            # unexpected behaviour of pygetline
            linecache.clearcache()
            lineno = self.start_lineno
            p = []
            for line in lines:
                for i, value in enumerate(line.split()):
                    p.append(
                        Value(
                            value=value,
                            lineno=lineno,
                            index=i,
                            file=self.fort10,
                        )
                    )
                lineno += 1

            c_index = 0
            for i in range(self.num_component):
                p[c_index].type(int)
                self.__constraint_num.append(p[c_index])
                num_constraints = np.abs(p[c_index].v)
                for num in range(num_constraints):
                    p[c_index + num + 1].type(int)
                    self.__constraint_index.append(i)
                    self.__basis_index.append(p[c_index + num + 1])

                c_index += num_constraints + 1

            self.read_flag = True

    def write(self):
        raise NotImplementedError

    @property
    def constraint_num(self):
        self.read()
        return [i.v for i in self.__constraint_num]

    @property
    def constraint_index(self):
        self.read()
        return [i for i in self.__constraint_index]

    @property
    def basis_index(self):
        self.read()
        return [i.v for i in self.__basis_index]

    @property
    def start_lineno(self):
        return pygrep_lineno(self.fort10, self.start_keyword) + 1

    @property
    def end_lineno(self):
        return pygrep_lineno(self.fort10, self.end_keyword) - 1


if __name__ == "__main__":
    from logging import getLogger

    logger = getLogger("pyturbo")
    logger.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter("%(name)s - %(levelname)s - %(lineno)d - %(message)s")
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples
    # pyturbo modules
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from utils.env import pyturbo_root

    os.chdir(os.path.join(pyturbo_root, "tests", "fort10"))
    fort10 = IO_fort10("fort.10_SiO2_complex_mo")
    logger.debug(fort10.f10detbasissets.mo_coefficient)
    logger.debug(fort10.f10detbasissets.mo_coefficient)
    logger.debug(fort10.f10detbasissets.mo_coefficient_imag)
    logger.debug(fort10.f10detbasissets.coefficient)
    logger.debug(fort10.pp_flag)
    coeff = fort10.f10detbasissets.mo_coefficient
    coeff[0][0] = 1000
    fort10.f10detbasissets.mo_coefficient = coeff
    fort10 = IO_fort10("fort.10_hydrogen")
    fort10 = IO_fort10("fort.10_benzene")
    fort10 = IO_fort10("fort.10_SiO2")
