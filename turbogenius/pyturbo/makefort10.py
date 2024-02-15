#!python -u
# -*- coding: utf-8 -*-

"""

pyturbo: makefort10 related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# python modules
import os
import re
import numpy as np
from decimal import Decimal
from typing import Optional

# set logger
from logging import getLogger, StreamHandler, Formatter

# pyturbo modules
from turbogenius.pyturbo.utils.env import pyturbo_data_dir
from turbogenius.pyturbo.utils.env import turbo_makefort10_run_command
from turbogenius.pyturbo.utils.utility import file_check, file_check_flag
from turbogenius.pyturbo.utils.utility import (
    return_num_twobody_and_flag_onebody,
)
from turbogenius.pyturbo.structure import Structure
from turbogenius.pyturbo.namelist import Namelist
from turbogenius.pyturbo.fortranIO import FortranIO
from turbogenius.pyturbo.pseudopotentials import Pseudopotentials
from turbogenius.pyturbo.basis_set import Jas_Basis_sets, Det_Basis_sets
from turbogenius.pyturbo.utils.execute import run

logger = getLogger("pyturbo").getChild(__name__)


class Makefort10(FortranIO):
    def __init__(
        self,
        structure: Optional[Structure] = None,
        det_basis_sets: Optional[Det_Basis_sets] = None,
        jas_basis_sets: Optional[Jas_Basis_sets] = None,
        pseudo_potentials: Optional[Pseudopotentials] = None,
        namelist: Optional[Namelist] = None,
    ):
        if structure is None:
            structure = Structure()
        if det_basis_sets is None:
            det_basis_sets = Det_Basis_sets()
        if jas_basis_sets is None:
            jas_basis_sets = Jas_Basis_sets()
        if pseudo_potentials is None:
            pseudo_potentials = Pseudopotentials()
        if namelist is None:
            namelist = Namelist()

        """
        input values
        """
        self.structure = structure
        self.det_basis_sets = det_basis_sets
        self.jas_basis_sets = jas_basis_sets
        self.pseudo_potentials = pseudo_potentials
        self.namelist = namelist

    def __str__(self):

        output = [
            "TurboRVB makefort10 python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        assert self.structure.natom == self.det_basis_sets.nuclei_num

    def generate_input(self, input_name: str, basis_sets_unique_element: bool = True):
        # pseudo potential generation (pseudo.dat)
        self.pseudo_potentials.write_pseudopotential_turborvb_file()

        # input_name
        self.namelist.write(input_name)
        logger.info(f"{input_name} has been generated.")

        ###############################
        # ATOMIC_POSITIONS
        ###############################
        output = []
        output.append("ATOMIC_POSITIONS \n")
        atomic_numbers_done = []
        atomic_numbers_done_counter = []
        dummy_atomic_number_shift = 0.01
        atomic_numbers_shifted_list = []

        for num in range(self.structure.natom):
            atomic_number = self.structure.atomic_numbers[num]

            if self.structure.has_celldm:
                logger.warning("Cartesian coord.")
                x, y, z = self.structure.positions[num]
            else:
                logger.warning("Fractional coord.")
                x, y, z = self.structure.positions_frac[num]

            # we need this operation even if basis_sets_unique_element is True
            # because special treaty is needed for PP case.
            if atomic_number in atomic_numbers_done:
                atomic_numbers_done_counter[
                    atomic_numbers_done.index(atomic_number)
                ] += 1
                shift = (
                    atomic_numbers_done_counter[
                        atomic_numbers_done.index(atomic_number)
                    ]
                    * dummy_atomic_number_shift
                )
            else:
                atomic_numbers_done_counter.append(1)
                shift = 1 * dummy_atomic_number_shift
                atomic_numbers_done.append(atomic_number)

            if num in self.pseudo_potentials.nucleus_index:
                # PP case
                z_core = self.pseudo_potentials.z_core[num]
                valence_electron = atomic_number - z_core
                if (
                    basis_sets_unique_element
                ):  # NOT only hydrogen and helium cases!! e.g., for Li, Be, etc...
                    # if z_core == 0 and (atomic_number == 1 or atomic_number == 2): # hydrogen and helium cases!!
                    if z_core == 0:  # PPs that do not remove any electron
                        shift = dummy_atomic_number_shift
                    else:
                        shift = 0.0
                else:
                    shift = shift

            else:
                # all-electron
                valence_electron = atomic_number
                if basis_sets_unique_element:
                    shift = 0.0  # reset shift!!
                else:
                    shift = shift

            atomic_numbers_shifted_list.append(shift)

            output.append(
                " {:.8f}  {:.8f}  {:.14f}  {:.14f}  {:.14f}\n".format(
                    valence_electron, atomic_number + shift, x, y, z
                )
            )

        output.append("/\n\n")

        with open(input_name, mode="a") as f:
            f.writelines(output)

        ###############################
        # Basis sets
        ###############################
        output = []
        atomic_numbers_done = []

        # check hybrid orbials are 0 or consistent with the number of nuclei
        if len(self.det_basis_sets.number_of_additional_hybrid_orbitals) != 0:
            det_add_hybridflag = True
            assert (
                len(self.det_basis_sets.number_of_additional_hybrid_orbitals)
                == self.structure.natom
            )
        else:
            det_add_hybridflag = False
        if len(self.jas_basis_sets.number_of_additional_hybrid_orbitals) != 0:
            jas_add_hybridflag = True
            assert (
                len(self.jas_basis_sets.number_of_additional_hybrid_orbitals)
                == self.structure.natom
            )
        else:
            jas_add_hybridflag = False

        for nucleus in range(self.structure.natom):
            atomic_number = self.structure.atomic_numbers[nucleus]
            if (
                basis_sets_unique_element is True
                and atomic_number in atomic_numbers_done
            ):
                continue
            shift = atomic_numbers_shifted_list[nucleus]
            fake_atomic_number = atomic_number + shift
            if fake_atomic_number.is_integer():
                label = "ATOM_{:d}".format(int(fake_atomic_number))
            else:
                # label = "ATOM_{:f}".format(fake_atomic_number) # :f format does not work.
                label = "ATOM_{:f}".format(
                    Decimal(str(round(fake_atomic_number, 8))).normalize()
                )
            nshelldet = len(
                [
                    i
                    for i, x in enumerate(self.det_basis_sets.nucleus_index)
                    if x == nucleus
                ]
            ) + len(
                [
                    i
                    for i, x in enumerate(self.det_basis_sets.hyb_nucleus_index)
                    if x == nucleus
                ]
            )
            nshelljas = len(
                [
                    i
                    for i, x in enumerate(self.jas_basis_sets.nucleus_index)
                    if x == nucleus
                ]
            )
            if det_add_hybridflag:
                ndet_hyb = self.det_basis_sets.number_of_additional_hybrid_orbitals[
                    nucleus
                ]
            else:
                ndet_hyb = 0
            if jas_add_hybridflag:
                njas_hyb = self.jas_basis_sets.number_of_additional_hybrid_orbitals[
                    nucleus
                ]
            else:
                njas_hyb = 0
            output.append(f"{label}\n")
            output.append("&shells\n")
            output.append(f" nshelldet={nshelldet}\n")
            output.append(f" nshelljas={nshelljas}\n")
            if ndet_hyb != 0:
                output.append(f" ndet_hyb={ndet_hyb}\n")
            if njas_hyb != 0:
                output.append(f" njas_hyb={njas_hyb}\n")
            output.append("/\n")

            for k, basis_sets in enumerate((self.det_basis_sets, self.jas_basis_sets)):
                # basis set
                shell_index = [
                    i for i, x in enumerate(basis_sets.nucleus_index) if x == nucleus
                ]
                for shell in shell_index:
                    shell_ang_mom = basis_sets.shell_ang_mom[shell]
                    shell_ang_mom_turbo = basis_sets.shell_ang_mom_turbo_notation[shell]
                    prim_index = [
                        i for i, x in enumerate(basis_sets.shell_index) if x == shell
                    ]

                    if len(prim_index) > 1:
                        logger.debug("Contracted shell!")
                        # contracted shell
                        exponent_list = [
                            basis_sets.exponent[prim] for prim in prim_index
                        ]
                        coefficient_list = [
                            basis_sets.coefficient[prim] for prim in prim_index
                        ]
                        if basis_sets.complex_flag:
                            coefficient_imag_list = [
                                basis_sets.coefficient_imag[prim] for prim in prim_index
                            ]
                        logger.debug(exponent_list)
                        logger.debug(coefficient_list)

                        num_param = len(exponent_list) + len(coefficient_list)

                        output.append(
                            "  {:d}   {:d}   {:d}\n".format(
                                2 * shell_ang_mom + 1,
                                num_param,
                                shell_ang_mom_turbo,
                            )
                        )
                        if basis_sets_unique_element:
                            nucleus_index_label = 1
                        else:
                            nucleus_index_label = nucleus + 1
                        output.append("  {:d}".format(nucleus_index_label))
                        for exponent in exponent_list:
                            output.append(" {:12f}".format(exponent))
                        if basis_sets.complex_flag:  # complex case
                            for coefficient, coefficient_imag in zip(
                                coefficient_list, coefficient_imag_list
                            ):
                                output.append(
                                    " {:12f} {:12f}".format(
                                        coefficient, coefficient_imag
                                    )
                                )
                        else:  # real case
                            for coefficient in coefficient_list:
                                output.append(" {:12f}".format(coefficient))
                        output.append("\n")
                    else:
                        logger.debug("Uncontracted shell!")
                        prim_index = prim_index[0]
                        # uncontracted shell
                        exponent = basis_sets.exponent[prim_index]
                        logger.debug(exponent)
                        output.append(
                            "  {:d}   {:d}   {:d}\n".format(
                                2 * shell_ang_mom + 1, 1, shell_ang_mom_turbo
                            )
                        )
                        if basis_sets_unique_element:
                            nucleus_index_label = 1
                        else:
                            nucleus_index_label = nucleus + 1
                        output.append(
                            "  {:d}   {:.12f}\n".format(nucleus_index_label, exponent)
                        )

                # hybrid orbitals (k==0):
                if k == 0:  # hybrid orbital:
                    hyb_index_list = [
                        i
                        for i, x in enumerate(basis_sets.hyb_nucleus_index)
                        if x == nucleus
                    ]
                    for hyb_index in hyb_index_list:
                        hyb_shell_ang_mom = basis_sets.hyb_shell_ang_mom[hyb_index]
                        hyb_shell_ang_mom_turbo = (
                            basis_sets.hyb_shell_ang_mom_turbo_notation[hyb_index]
                        )
                        hyb_param_num = basis_sets.hyb_param_num[hyb_index]

                        output.append(
                            "  {:d}   {:d}   {:d}\n".format(
                                2 * hyb_shell_ang_mom + 1,
                                hyb_param_num,
                                hyb_shell_ang_mom_turbo,
                            )
                        )

                        if basis_sets_unique_element:
                            nucleus_index_label = 1
                        else:
                            nucleus_index_label = nucleus + 1
                        output.append("  {:d}".format(nucleus_index_label))

                        min_hyb_prim_index = np.min(
                            basis_sets.hyb_prim_index[hyb_index]
                        )
                        for hyb_prim_index in basis_sets.hyb_prim_index[hyb_index]:
                            output.append(
                                " {:d}".format(hyb_prim_index - min_hyb_prim_index + 1)
                            )
                        output.append("\n")
                        if basis_sets.complex_flag:  # complex case
                            for coefficient, coefficient_imag in zip(
                                basis_sets.hyb_coefficient[hyb_index],
                                basis_sets.hyb_coefficient_imag[hyb_index],
                            ):
                                output.append(
                                    " {:12f} {:12f}".format(
                                        coefficient, coefficient_imag
                                    )
                                )
                            output.append("\n")
                        else:
                            for coefficient in basis_sets.hyb_coefficient[hyb_index]:
                                output.append(" {:12f}".format(coefficient))
                            output.append("\n")

                if k == 0:  # det part:
                    output.append("# Parameters atomic Jastrow wf \n")
                else:  # jas part
                    output.append("\n")

            atomic_numbers_done.append(atomic_number)

        with open(input_name, "a") as f:
            f.writelines(output)

    def run(
        self,
        input_name: str = "makefort10.input",
        output_name: str = "out_make",
    ):
        run(
            turbo_makefort10_run_command,
            input_name=input_name,
            output_name=output_name,
        )

    def check_results(self, output_names: Optional[list] = None):
        if output_names is None:
            output_names = ["out_make"]
        flags = []
        for output_name in output_names:
            file_check(output_name)
            with open(output_name, "r") as f:
                lines = f.readlines()
            if any([re.match(r".*Effective.*rs.*", line) for line in lines]):
                if file_check_flag("fort.10_new"):
                    flags.append(True)
            else:
                flags.append(False)
        return flags

    @staticmethod
    def read_default_namelist(
        structure: Optional[Structure] = None,
        jastrow_type: int = -6,
        neldiff: int = 0,
        nel: int = -1,
        complex_flag: bool = False,
    ):
        if structure is None:
            structure = Structure()
        makefort10_default_file = os.path.join(
            pyturbo_data_dir, "makefort10", "makefort10.input"
        )
        namelist = Namelist.parse_namelist_from_file(makefort10_default_file)

        namelist.set_parameter("natoms", structure.natom, "&system")
        namelist.set_parameter("ntyp", structure.ntyp, "&system")

        # PBC case
        if structure.pbc_flag:
            # use celldm() notation for makefort10.input
            if structure.has_celldm:
                namelist.set_parameter("posunits", "bohr", "&system")
                # non-orthorhombic cell
                if structure.tilted_flag:
                    namelist.set_parameter(
                        "celldm(1)", structure.cell.celldm_1, "&system"
                    )
                    namelist.set_parameter(
                        "celldm(2)", structure.cell.celldm_2, "&system"
                    )
                    namelist.set_parameter(
                        "celldm(3)", structure.cell.celldm_3, "&system"
                    )
                    namelist.set_parameter(
                        "celldm(4)", structure.cell.celldm_4, "&system"
                    )
                    namelist.set_parameter(
                        "celldm(5)", structure.cell.celldm_5, "&system"
                    )
                    namelist.set_parameter(
                        "celldm(6)", structure.cell.celldm_6, "&system"
                    )
                    namelist.set_parameter("yes_tilted", ".true.", "&system")
                # orthorhomic cell
                else:
                    namelist.set_parameter(
                        "celldm(1)", structure.cell.celldm_1, "&system"
                    )
                    namelist.set_parameter(
                        "celldm(2)", structure.cell.celldm_2, "&system"
                    )
                    namelist.set_parameter(
                        "celldm(3)", structure.cell.celldm_3, "&system"
                    )
                    namelist.comment_out("celldm(4)")
                    namelist.comment_out("celldm(5)")
                    namelist.comment_out("celldm(6)")
                    namelist.set_parameter("yes_tilted", ".false.", "&system")

            # use at() notation for makefort10.input
            else:
                namelist.set_parameter("posunits", "crystal", "&system")
                norm_ref = structure.cell.norm_vec_a  # bohr
                at_1_1 = structure.cell.vec_a[0] / norm_ref
                at_2_1 = structure.cell.vec_a[1] / norm_ref
                at_3_1 = structure.cell.vec_a[2] / norm_ref

                at_1_2 = structure.cell.vec_b[0] / norm_ref
                at_2_2 = structure.cell.vec_b[1] / norm_ref
                at_3_2 = structure.cell.vec_b[2] / norm_ref

                at_1_3 = structure.cell.vec_c[0] / norm_ref
                at_2_3 = structure.cell.vec_c[1] / norm_ref
                at_3_3 = structure.cell.vec_c[2] / norm_ref

                # Bohr
                namelist.set_parameter("L_read", norm_ref, "&system")
                namelist.set_parameter("at(1,1)", at_1_1, "&system")
                namelist.set_parameter("at(2,1)", at_2_1, "&system")
                namelist.set_parameter("at(3,1)", at_3_1, "&system")
                namelist.set_parameter("at(1,2)", at_1_2, "&system")
                namelist.set_parameter("at(2,2)", at_2_2, "&system")
                namelist.set_parameter("at(3,2)", at_3_2, "&system")
                namelist.set_parameter("at(1,3)", at_1_3, "&system")
                namelist.set_parameter("at(2,3)", at_2_3, "&system")
                namelist.set_parameter("at(3,3)", at_3_3, "&system")

                namelist.comment_out("celldm(1)")
                namelist.comment_out("celldm(2)")
                namelist.comment_out("celldm(3)")
                namelist.comment_out("celldm(4)")
                namelist.comment_out("celldm(5)")
                namelist.comment_out("celldm(6)")

                # non-orthorhombic cell
                if structure.tilted_flag:
                    namelist.set_parameter("yes_tilted", ".true.", "&system")
                # orthorhomic cell
                else:
                    namelist.set_parameter("yes_tilted", ".false.", "&system")

            namelist.set_parameter("pbcfort10", ".true.", "&system")
            namelist.set_parameter("yes_crystal", ".true.", "&electrons")
            namelist.set_parameter("yes_crystalj", ".true.", "&electrons")

            if complex_flag:
                namelist.set_parameter("complexfort10", ".true.", "&system")
            else:
                namelist.set_parameter("complexfort10", ".false.", "&system")

        # Open system
        else:
            namelist.set_parameter("posunits", "bohr", "&system")
            if complex_flag:
                raise ValueError
            namelist.set_parameter("pbcfort10", ".false.", "&system")
            namelist.set_parameter("yes_crystal", ".false.", "&electrons")
            namelist.set_parameter("yes_crystalj", ".false.", "&electrons")
            namelist.set_parameter("complexfort10", ".false.", "&system")

        # onebody and twobody Jastrow parts
        namelist.set_parameter("twobody", jastrow_type, "&electrons")
        num_twobody, flag_onebody = return_num_twobody_and_flag_onebody(
            jastrow_type=jastrow_type
        )

        # onebody part
        if flag_onebody:
            for ntyp in range(structure.ntyp):
                namelist.set_parameter(f"onebodypar({ntyp + 1})", 1.0, "&electrons")
            namelist.comment_out("onebodypar")
        else:
            namelist.comment_out("onebodypar")

        # twobody part
        for i in range(num_twobody):
            namelist.set_parameter(f"twobodypar({i+1})", 1.0, "&electrons")
        namelist.comment_out("twobodypar")

        # spin (i.e., neldiff)
        namelist.set_parameter("neldiff", neldiff, "&electrons")

        # electrons (i.e., number of electrons)
        if nel != -1:
            namelist.set_parameter("nel", nel, "&electrons")

        return namelist

    @staticmethod
    def read_namelist_from_file(file):
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist

    @classmethod
    def parse_from_default_namelist(cls):
        raise NotImplementedError

    @classmethod
    def parse_from_file(cls, file):
        raise NotImplementedError


if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter("%(name)s - %(levelname)s - %(lineno)d - %(message)s")
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples
