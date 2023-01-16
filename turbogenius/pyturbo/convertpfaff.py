#!python
# -*- coding: utf-8 -*-
"""

pyturbo: convertpfaff related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# pyturbo modules
from turbogenius.pyturbo.fortranIO import FortranIO
from turbogenius.pyturbo.utils.env import turbo_convertfortpfaff_run_command
from turbogenius.pyturbo.utils.utility import file_check
from turbogenius.pyturbo.utils.execute import run
from turbogenius.pyturbo.io_fort10 import IO_fort10

from logging import getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)


class Convertpfaff(FortranIO):
    def __init__(self, in_fort10="fort.10_in", out_fort10="fort.10_out"):

        """
        input values
        """
        file_check(in_fort10)
        file_check(out_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag:
            file_check("pseudo.dat")
        self.in_fort10 = in_fort10
        self.out_fort10 = out_fort10

    def __str__(self):

        output = [
            "TurboRVB convertpfaff python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        pass

    def generate_input(self):
        logger.info(f"No input has been generated. \n")

    def run(
        self,
        rotate_flag=False,
        rotate_angle=0,
        scale_mean_field=1000,
        output_name="out_pfaff",
    ):

        if not rotate_flag:
            fort10 = IO_fort10("fort.10_in")
            if (
                fort10.f10header.nelup != fort10.f10header.neldn
            ):  # spin is fintite
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

    def check_results(self, output_names=["out_pfaff"]):
        return [True]

    @staticmethod
    def read_default_namelist():
        pass

    @staticmethod
    def read_namelist_from_file():
        pass

    @classmethod
    def parse_from_default_namelist(
        cls, in_fort10="fort.10_in", out_fort10="fort.10_out"
    ):
        return cls(in_fort10=in_fort10, out_fort10=out_fort10)

    @classmethod
    def parse_from_file(cls, in_fort10="fort.10_in", out_fort10="fort.10_out"):
        return cls(in_fort10=in_fort10, out_fort10=out_fort10)


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
