#!python
# -*- coding: utf-8 -*-

"""

convertfort10mol genius related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# python modules
import os
import numpy as np

# Logger
from logging import getLogger, StreamHandler, Formatter

# pyturbo modules
from turbogenius.pyturbo.convertfort10mol import Convertfort10mol
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.utils_workflows.env import turbo_genius_root
from turbogenius.geniusIO import GeniusIO

logger = getLogger("Turbo-Genius").getChild(__name__)


class Convertfort10mol_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo convertfortmol class

    Attributes:
         fort10 (str): input fort.10 WF file
         add_random_mo (bool): flag to add random MOs
         additional_mo (int): The number of added MOs
         grid_size (float): grid size for x,y,z (Bohr)
    """

    def __init__(
        self,
        fort10="fort.10_in",
        add_random_mo=True,
        grid_size=0.10,
        additional_mo=0,
    ):

        self.fort10 = fort10
        self.grid_size = grid_size

        self.convertfort10mol = Convertfort10mol.parse_from_default_namelist(
            in_fort10=fort10
        )
        self.io_fort10 = IO_fort10(self.fort10)
        nmol = self.io_fort10.f10header.neldn + additional_mo

        self.convertfort10mol.set_parameter(
            parameter="nmol", value=nmol, namelist="&molec_info"
        )

        if self.io_fort10.f10structure.pbc_flag:
            self.convertfort10mol.comment_out(parameter="ax")
            self.convertfort10mol.comment_out(parameter="ay")
            self.convertfort10mol.comment_out(parameter="az")

        if not add_random_mo:
            self.convertfort10mol.set_parameter(
                parameter="epsdgm", value=1.0e-14, namelist="&control"
            )
            # set L_box
            if self.io_fort10.f10structure.pbc_flag:
                # for crystals, Lx, Ly, and Lz are cells
                Lx = self.io_fort10.f10structure.norm_vec_a
                Ly = self.io_fort10.f10structure.norm_vec_b
                Lz = self.io_fort10.f10structure.norm_vec_c
                logger.info("Lbox is the norms of the lattice vectors")
                logger.info(f"Lx={Lx}, Ly={Ly}, Lz={Lz}")
                ax = self.grid_size
                ay = self.grid_size
                az = self.grid_size
                self.convertfort10mol.set_parameter(
                    parameter="nx", value=ax, namelist="&mesh_info"
                )
                self.convertfort10mol.set_parameter(
                    parameter="ny", value=ay, namelist="&mesh_info"
                )
                self.convertfort10mol.set_parameter(
                    parameter="nz", value=az, namelist="&mesh_info"
                )
                self.convertfort10mol.comment_out(parameter="ax")
                self.convertfort10mol.comment_out(parameter="ay")
                self.convertfort10mol.comment_out(parameter="az")
            else:
                # +- 7.5 bohr from the edges.
                pos = self.io_fort10.f10structure.positions
                Lx = np.max(pos[:, 0]) - np.min(pos[:, 0]) + 15.0
                Ly = np.max(pos[:, 1]) - np.min(pos[:, 1]) + 15.0
                Lz = np.max(pos[:, 2]) - np.min(pos[:, 2]) + 15.0
                logger.info(
                    "Lbox is set to +- 7.5 bohr from the edges of the molecules."
                )
                logger.info(f"Lx={Lx}, Ly={Ly}, Lz={Lz}")
                ax = self.grid_size
                ay = self.grid_size
                az = self.grid_size
                nx = int(Lx / ax)
                ny = int(Ly / ay)
                nz = int(Lz / az)
                logger.info(f"nx={nx}, ny={ny}, nz={nz}")
                self.convertfort10mol.set_parameter(
                    parameter="ax", value=ax, namelist="&mesh_info"
                )
                self.convertfort10mol.set_parameter(
                    parameter="ay", value=ay, namelist="&mesh_info"
                )
                self.convertfort10mol.set_parameter(
                    parameter="az", value=az, namelist="&mesh_info"
                )
                self.convertfort10mol.set_parameter(
                    parameter="nx", value=nx, namelist="&mesh_info"
                )
                self.convertfort10mol.set_parameter(
                    parameter="ny", value=ny, namelist="&mesh_info"
                )
                self.convertfort10mol.set_parameter(
                    parameter="nz", value=nz, namelist="&mesh_info"
                )

        else:
            logger.info("Random MO coefficients will be added to fort.10_in")

        self.energy = None

    def run_all(
        self,
        input_name: str = "convertfort10mol.input",
        output_name: str = "out_mol",
    ) -> None:
        """
        Generate input files and run the command.

        Args:
            input_name (str): input file name
            output_name (str): output file name

        """
        self.generate_input(input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)

    def generate_input(
        self, input_name: str = "convertfort10mol.input"
    ) -> None:
        """
        Generate input file.

        Args:
            input_name (str): input file name

        """
        self.convertfort10mol.generate_input(input_name=input_name)

    def run(
        self,
        input_name: str = "convertfort10mol.input",
        output_name: str = "out_mol",
    ):
        """
        Run the command.

        Args:
            input_name (str): input file name
            output_name (str): output file name
        """
        self.convertfort10mol.run(
            input_name=input_name, output_name=output_name
        )
        flags = self.convertfort10mol.check_results(output_names=[output_name])
        assert all(flags)

    def check_results(self, output_names: list = ["out_mol"]) -> bool:
        """
        Check the result.

        Args:
            output_names (list): a list of output file names
        Return:
            bool: True if all the runs were successful, False if an error is detected in the files.
        """
        return self.convertfort10mol.check_results(output_names=output_names)


if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "convertfort10mol"))

    # moved to examples
    convertfort10mol_genius = Convertfort10mol_genius()
    convertfort10mol_genius.generate_input()
