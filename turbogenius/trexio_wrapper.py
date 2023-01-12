#!/usr/bin/env python
# coding: utf-8

"""

TREXIO related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# import python modules
import os

# logger
from logging import getLogger, StreamHandler, Formatter

# import trexio
import trexio

logger = getLogger("Turbo-Genius").getChild(__name__)


class Trexio_wrapper_r:
    """

    This class is a wrapper for the TREXIO program

    Attributes:
         trexio_file (str): name of TREXIO file

    """

    def __init__(self, trexio_file):
        # prefix and file names
        logger.info(f"TREXIO file = {trexio_file}")

        self.hdf5_filename = trexio_file
        file_r = trexio.File(
            os.path.join(self.hdf5_filename),
            mode="r",
            back_end=trexio.TREXIO_HDF5,
        )

        # check if the system is PBC or not.
        self.periodic = trexio.read_pbc_periodic(file_r)

        if self.periodic:
            logger.info("Crystal (Periodic boundary condition)")
            self.cell_a = trexio.read_cell_a(file_r)
            self.cell_b = trexio.read_cell_b(file_r)
            self.cell_c = trexio.read_cell_c(file_r)
            self.k_point = trexio.read_pbc_k_point(file_r)
        else:
            logger.info("Molecule (Open boundary condition)")

        # read electron num
        self.num_ele_up = trexio.read_electron_up_num(file_r)
        self.num_ele_dn = trexio.read_electron_dn_num(file_r)
        self.num_ele_total = self.num_ele_up + self.num_ele_dn

        # read structure info.
        self.nucleus_num_r = trexio.read_nucleus_num(file_r)
        self.labels_r = trexio.read_nucleus_label(file_r)
        self.charges_r = trexio.read_nucleus_charge(file_r)
        self.coords_r = trexio.read_nucleus_coord(file_r)

        # Reading basis sets info
        self.basis_type = trexio.read_basis_type(file_r)
        self.basis_shell_num = trexio.read_basis_shell_num(file_r)
        self.basis_shell_index = trexio.read_basis_shell_index(file_r)
        self.basis_prim_num = trexio.read_basis_prim_num(file_r)
        self.basis_nucleus_index = trexio.read_basis_nucleus_index(file_r)
        self.basis_shell_ang_mom = trexio.read_basis_shell_ang_mom(file_r)
        self.basis_shell_factor = trexio.read_basis_shell_factor(file_r)
        self.basis_shell_index = trexio.read_basis_shell_index(file_r)
        self.basis_exponent = trexio.read_basis_exponent(file_r)
        self.basis_coefficient = trexio.read_basis_coefficient(file_r)
        self.basis_prim_factor = trexio.read_basis_prim_factor(file_r)

        # Pseudo potentials info
        self.ecp_max_ang_mom_plus_1 = trexio.read_ecp_max_ang_mom_plus_1(
            file_r
        )
        self.ecp_z_core = trexio.read_ecp_z_core(file_r)
        self.ecp_num = trexio.read_ecp_num(file_r)
        self.ecp_ang_mom = trexio.read_ecp_ang_mom(file_r)
        self.ecp_nucleus_index = trexio.read_ecp_nucleus_index(file_r)
        self.ecp_exponent = trexio.read_ecp_exponent(file_r)
        self.ecp_coefficient = trexio.read_ecp_coefficient(file_r)
        self.ecp_power = trexio.read_ecp_power(file_r)

        # ao info
        self.ao_cartesian = trexio.read_ao_cartesian(file_r)
        self.ao_num = trexio.read_ao_num(file_r)
        self.ao_shell = trexio.read_ao_shell(file_r)
        self.ao_normalization = trexio.read_ao_normalization(file_r)

        # mo info
        self.mo_type = trexio.read_mo_type(file_r)
        self.mo_num = trexio.read_mo_num(file_r)
        self.mo_occupation = trexio.read_mo_occupation(file_r)
        self.mo_coefficient = trexio.read_mo_coefficient(file_r)
        if trexio.has_mo_coefficient_im(file_r):
            logger.info("The WF is complex")
            self.mo_coefficient_imag = trexio.read_mo_coefficient_im(file_r)
            self.complex_flag = True
        else:
            logger.info("The WF is real")
            self.complex_flag = False

        file_r.close()


if __name__ == "__main__":
    import sys

    logger = getLogger("Turbo-Genius")
    logger.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples
    sys.path.append(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")
    )
    from utils_workflows.env import turbo_genius_root

    trexio_test_dir = os.path.join(
        turbo_genius_root, "tests", "trexio_to_turborvb"
    )

    os.chdir(trexio_test_dir)

    # diamond
    Trexio_wrapper_r(trexio_file="diamond_single_k.hdf5")
