#!python
# -*- coding: utf-8 -*-

"""

pyturbo: convertfort10 related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# python modules
import os
import re
from typing import Optional

# pyturbo modules
from turbogenius.pyturbo.namelist import Namelist
from turbogenius.pyturbo.fortranIO import FortranIO
from turbogenius.pyturbo.utils.env import turbo_convertfort10_run_command
from turbogenius.pyturbo.utils.env import pyturbo_data_dir
from turbogenius.pyturbo.utils.utility import file_check
from turbogenius.pyturbo.utils.execute import run
from turbogenius.pyturbo.io_fort10 import IO_fort10


from logging import getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)


class Convertfort10(FortranIO):
    """

    This class is a wrapper of turborvb convertfort10.x

    Attributes:
         in_fort10 (str): input fort.10
         out_fort10 (str): template fort.10
         namelist (Namelist): fortran namelist for convertfort10.x

    """

    def __init__(
        self,
        in_fort10: str = "fort.10_in",
        out_fort10: str = "fort.10_out",
        namelist: Optional[Namelist] = None,
    ):
        if namelist is None:
            namelist = Namelist()

        """
        input values
        """
        file_check(in_fort10)
        file_check(out_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag:
            file_check("pseudo.dat")
        self.in_fort10 = in_fort10
        self.out_fort10 = out_fort10
        self.namelist = namelist

    def __str__(self):

        output = [
            "TurboRVB convertfort10 python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self) -> None:
        """
        Sanity check

        """
        pass

    def generate_input(self, input_name: str = "convertfort10.input") -> None:
        """
        Generate input file.

        Args:
            input_name (str): input file name
        """
        self.namelist.write(input_name)
        logger.info(f"{input_name} has been generated.")

    def run(
        self,
        input_name: str = "convertfort10.input",
        output_name: str = "out_conv",
    ):
        """
        Run the command.

        Args:
            input_name (str): input file name
            output_name (str): output file name
        """
        run(
            turbo_convertfort10_run_command,
            input_name=input_name,
            output_name=output_name,
        )

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
        flags = []
        for output_name in output_names:
            file_check(output_name)
            with open(output_name, "r") as f:
                lines = f.readlines()
            if any([re.match(r".*Overlap.*square.*", line) for line in lines]):
                flags.append(True)
            else:
                flags.append(False)
        return flags

    @staticmethod
    def read_default_namelist(
        in_fort10: str = "fort.10_in", out_fort10: str = "fort.10_out"
    ):  # -> namelist
        """
        Read default namelist values from turbogenius database

        Args:
            in_fort10 (str): input fort.10
            out_fort10 (str): template fort.10
        Returns:
            Namelist: default namelist values taken from the database
        """
        convertfort10_default_file = os.path.join(
            pyturbo_data_dir, "convertfort10", "convertfort10.input"
        )
        namelist = Namelist.parse_namelist_from_file(
            convertfort10_default_file
        )
        # To be implemented
        # set ax,ay,az,nx,ny,nz, depending of fort.10s
        return namelist

    @staticmethod
    def read_namelist_from_file(file: str):  # -> namelist
        """
        Read namelist values from a specified file

        Args:
            file (str): filename
        Returns:
            Namelist: namelist values read from the specified file
        """
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist

    @classmethod
    def parse_from_default_namelist(
        cls, in_fort10: str = "fort.10_in", out_fort10: str = "fort.10_out"
    ):  # -> cls
        """
        Read default namelist values from turbogenius database

        Args:
            in_fort10 (str): input fort.10
            out_fort10 (str): template fort.10
        Returns:
            cls: cls with default namelist values taken from the database
        """
        namelist = cls.read_default_namelist(
            in_fort10=in_fort10, out_fort10=out_fort10
        )
        return cls(
            in_fort10=in_fort10, out_fort10=out_fort10, namelist=namelist
        )

    @classmethod
    def parse_from_file(
        cls,
        file: str,
        in_fort10: str = "fort.10_in",
        out_fort10: str = "fort.10_out",
    ):  # -> cls
        """
        Read namelist values from a specified file

        Args:
            file (str): filename
            in_fort10 (str): input fort.10
            out_fort10 (str): template fort.10
        Returns:
            cls: cls with namelist values read from the specified file
        """
        namelist = Namelist.parse_namelist_from_file(file)
        return cls(
            in_fort10=in_fort10, out_fort10=out_fort10, namelist=namelist
        )


if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples.
