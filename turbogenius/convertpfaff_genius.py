#!python
# -*- coding: utf-8 -*-

"""

convertpfaff genius related classes and methods


"""

# python modules
import os
from typing import Optional

# Logger
from logging import getLogger, StreamHandler, Formatter

# turbogenius modules
from turbogenius.pyturbo.convertpfaff import Convertpfaff
from turbogenius.utils_workflows.env import turbo_genius_root
from turbogenius.geniusIO import GeniusIO

logger = getLogger("Turbo-Genius").getChild(__name__)


class Convertpfaff_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo Convertpfaff class

    Attributes:
         in_fort10 (str): fort.10 WF file (input WF file)
         out_fort10 (str): fort.10 WF file (template WF file)
    """

    def __init__(
        self,
        in_fort10: str = "fort.10_in",
        out_fort10: str = "fort.10_out",
    ):

        self.in_fort10 = in_fort10
        self.out_fort10 = out_fort10

        self.convertpfaff = Convertpfaff.parse_from_default_namelist(
            in_fort10=in_fort10, out_fort10=out_fort10
        )

    def run_all(
        self,
        rotate_flag: bool = False,
        rotate_angle: float = 0,
        scale_mean_field: int = 1000,
        output_name: str = "out_pfaff",
    ) -> None:
        """
        Generate input files and run the command.

        Args:
            rotate_flag (bool): rotation flag, True or False
            rotate_angle (float): rotation angle
            scale_mean_field (int): scaling mean field
            output_name (str): output file name

        """
        self.generate_input()
        self.run(
            rotate_flag=rotate_flag,
            rotate_angle=rotate_angle,
            scale_mean_field=scale_mean_field,
            output_name=output_name,
        )

    def run(
        self,
        rotate_flag: float = False,
        rotate_angle: float = 0,
        scale_mean_field: int = 1000,
        output_name: str = "out_pfaff",
    ) -> None:
        """
        Run the command.

        Args:
            rotate_flag (bool): rotation flag, True or False
            rotate_angle (float): rotation angle
            scale_mean_field (int): scaling mean field
            output_name (str): output file name
        """
        self.convertpfaff.run(
            rotate_flag=rotate_flag,
            rotate_angle=rotate_angle,
            scale_mean_field=scale_mean_field,
            output_name=output_name,
        )
        flags = self.convertpfaff.check_results(output_names=[output_name])
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
            output_names = ["out_pfaff"]
        return self.convertpfaff.check_results(output_names=output_names)


if __name__ == "__main__":
    from turbogenius.convertfort10_genius import Convertfort10_genius

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
