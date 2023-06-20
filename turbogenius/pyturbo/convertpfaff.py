#!python
# -*- coding: utf-8 -*-
"""

pyturbo: convertpfaff related classes and methods

"""

# python modules
from typing import Optional

# pyturbo modules
from turbogenius.pyturbo.fortranIO import FortranIO
from turbogenius.pyturbo.utils.env import turbo_convertfortpfaff_run_command
from turbogenius.pyturbo.utils.utility import file_check
from turbogenius.pyturbo.utils.execute import run
from turbogenius.pyturbo.io_fort10 import IO_fort10

from logging import getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)


class Convertpfaff(FortranIO):
    """

    This class is a wrapper of turborvb convertpfaff.x

    Attributes:
         in_fort10 (str): input fort.10
         out_fort10 (str): output fort.10

    """

    def __init__(self, in_fort10: str = "fort.10_in", out_fort10: str = "fort.10_out"):

        """
        input values
        """
        file_check(in_fort10)
        file_check(out_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag:
            file_check("pseudo.dat")
        self.in_fort10 = in_fort10
        self.out_fort10 = out_fort10

    def __str__(self) -> str:
        return str(f"{self.__class__.__name__} class")

    def sanity_check(self) -> None:
        """
        Sanity check (to be implemented.)

        """
        pass

    def generate_input(self) -> None:
        """
        Generate input file.

        Args:

        """
        logger.info("No input has been generated. \n")

    def run(
        self,
        rotate_flag: bool = False,
        rotate_angle: int = 0,
        scale_mean_field: int = 1000,
        output_name: str = "out_pfaff",
    ):
        """
        Run the command.

        Args:
            rotate_flag (bool): rotate the angle of spins
            rotate_angle (int): rotation angle (* pi radian)
            scale_mean_filed (int): to be filled
            output_name (str): output file name
        """
        if not rotate_flag:
            fort10 = IO_fort10("fort.10_in")
            if fort10.f10header.nelup != fort10.f10header.neldn:  # spin is fintite
                logger.info(
                    f"norate for an unpaired case. put scale_mean_field={scale_mean_field}"
                )
                cmd = f"echo {scale_mean_field} | {turbo_convertfortpfaff_run_command} norotate"
            else:  # closed-shell
                cmd = f"{turbo_convertfortpfaff_run_command} norotate"

        if rotate_flag:
            cmd = f"echo {rotate_angle} | {turbo_convertfortpfaff_run_command} rotate"

        logger.info(f"cmd = {cmd}")
        run(cmd, input_name=None, output_name=output_name)

    def check_results(self, output_names: Optional[list] = None) -> list:
        """
        Check the result.

        Args:
            output_names (list): a list of output file names
        Returns:
            bool: True if all the runs were successful, False if an error is detected in the files.
        """
        if output_names is None:
            output_names = ["out_pfaff"]
        return [True]

    @staticmethod
    def read_default_namelist() -> None:
        """
        Read default namelist values from turbogenius database

        Args:
        """
        pass

    @staticmethod
    def read_namelist_from_file() -> None:
        """
        Read namelist values from a specified file

        Args:
        """
        pass

    @classmethod
    def parse_from_default_namelist(
        cls, in_fort10: str = "fort.10_in", out_fort10: str = "fort.10_out"
    ):
        """
        Read default namelist values from turbogenius database

        Args:
            in_fort10 (str): input fort.10
            out_fort10 (str): template fort.10
        Returns:
            cls: cls with default namelist values taken from the database
        """
        return cls(in_fort10=in_fort10, out_fort10=out_fort10)

    @classmethod
    def parse_from_file(
        cls,
        file: Optional[str] = None,
        in_fort10: str = "fort.10_in",
        out_fort10: str = "fort.10_out",
    ):
        """
        Read namelist values from a specified file

        Args:
            file (str): filename
            in_fort10 (str): input fort.10
            out_fort10 (str): template fort.10
        Returns:
            cls: cls with namelist values read from the specified file
        """
        return cls(in_fort10=in_fort10, out_fort10=out_fort10)


if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter("%(name)s - %(levelname)s - %(lineno)d - %(message)s")
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples.
