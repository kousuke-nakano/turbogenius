#!python
# -*- coding: utf-8 -*-

"""

convertfort10 genius related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# python modules
import os
import numpy as np
from typing import Optional

# Logger
from logging import getLogger, StreamHandler, Formatter

# pyturbo modules
from turbogenius.pyturbo.convertfort10 import Convertfort10
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.utils_workflows.env import turbo_genius_root
from turbogenius.geniusIO import GeniusIO

logger = getLogger("Turbo-Genius").getChild(__name__)


class Convertfort10_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo convertfort10 class

    Attributes:
         in_fort10 (str): fort.10 WF file (input)
         out_fort10 (str): fort.10 WF file (template)
         grid_size (float): grid size for xyz (bohr)
    """

    def __init__(
        self,
        in_fort10: str = "fort.10_in",
        out_fort10: str = "fort.10_out",
        grid_size: float = 0.10,
    ):

        self.in_fort10 = in_fort10
        self.out_fort10 = out_fort10
        self.grid_size = grid_size

        self.convertfort10 = Convertfort10.parse_from_default_namelist(
            in_fort10=in_fort10, out_fort10=out_fort10
        )
        self.io_fort10 = IO_fort10(self.in_fort10)

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
            nx = int(Lx / ax)
            ny = int(Ly / ay)
            nz = int(Lz / az)
            if nx%2 != 0:
                nx+=1
            if ny%2 != 0:
                ny+=1
            if nz%2 != 0:
                nz+=1
            self.convertfort10.set_parameter(parameter="nx", value=nx, namelist="&mesh_info")
            self.convertfort10.set_parameter(parameter="ny", value=ny, namelist="&mesh_info")
            self.convertfort10.set_parameter(parameter="nz", value=nz, namelist="&mesh_info")
            self.convertfort10.comment_out(parameter="ax")
            self.convertfort10.comment_out(parameter="ay")
            self.convertfort10.comment_out(parameter="az")
        else:
            pos = self.io_fort10.f10structure.positions
            Lx = np.max(pos[:, 0]) - np.min(pos[:, 0]) + 6.0
            Ly = np.max(pos[:, 1]) - np.min(pos[:, 1]) + 6.0
            Lz = np.max(pos[:, 2]) - np.min(pos[:, 2]) + 6.0
            logger.info(
                "Lbox is set to +- 3.0 bohr from the edges of the molecules."
            )
            logger.info(f"Lx={Lx}, Ly={Ly}, Lz={Lz}")
            ax = self.grid_size
            ay = self.grid_size
            az = self.grid_size
            nx = int(Lx / ax)
            ny = int(Ly / ay)
            nz = int(Lz / az)
            if nx%2 != 0:
                nx+=1
            if ny%2 != 0:
                ny+=1
            if nz%2 != 0:
                nz+=1
            logger.info(f"nx={nx}, ny={ny}, nz={nz}")
            self.convertfort10.set_parameter(
                parameter="ax", value=ax, namelist="&mesh_info"
            )
            self.convertfort10.set_parameter(
                parameter="ay", value=ay, namelist="&mesh_info"
            )
            self.convertfort10.set_parameter(
                parameter="az", value=az, namelist="&mesh_info"
            )
            self.convertfort10.set_parameter(
                parameter="nx", value=nx, namelist="&mesh_info"
            )
            self.convertfort10.set_parameter(
                parameter="ny", value=ny, namelist="&mesh_info"
            )
            self.convertfort10.set_parameter(
                parameter="nz", value=nz, namelist="&mesh_info"
            )

    def run_all(
        self,
        input_name: str = "convertfort10.input",
        output_name: str = "out_conv",
    ) -> None:
        """
        Generate input files and run the command.

        Args:
            input_name (str): input file name
            output_name (str): output file name

        """
        self.generate_input(input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)

    def generate_input(self, input_name: str = "convertfort10.input") -> None:
        """
        Generate input file.

        Args:
            input_name (str): input file name

        """
        self.convertfort10.generate_input(input_name=input_name)

    def run(self, input_name="convertfort10.input", output_name="out_conv"):
        """
        Run the command.

        Args:
            input_name (str): input file name
            output_name (str): output file name
        """
        self.convertfort10.run(input_name=input_name, output_name=output_name)
        flags = self.convertfort10.check_results(output_names=[output_name])
        assert all(flags)

    def check_results(self, output_names: Optional[list] = None) -> bool:
        """
        Check the result.

        Args:
            output_names (list): a list of output file names
        Returns:
            bool: True if all the runs were successful, False if an error is detected in the files.
        """
        if output_names is None:
            output_names = ["out_conv"]
        return self.convertfort10.check_results(output_names=output_names)


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

    os.chdir(os.path.join(turbo_genius_root, "tests", "convertfort10"))

    # moved to examples
    convertfort10_genius = Convertfort10_genius()
    convertfort10_genius.generate_input()
    convertfort10_genius.run()
