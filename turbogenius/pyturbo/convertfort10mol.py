#!python
# -*- coding: utf-8 -*-

"""

pyturbo: convertfort10mol related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""
# python modules
import os
import re

# pyturbo modules
from turbogenius.pyturbo.namelist import Namelist
from turbogenius.pyturbo.fortranIO import FortranIO
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.pyturbo.utils.env import turbo_convertfort10mol_run_command
from turbogenius.pyturbo.utils.env import pyturbo_data_dir
from turbogenius.pyturbo.utils.utility import file_check
from turbogenius.pyturbo.utils.execute import run

from logging import getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)


class Convertfort10mol(FortranIO):
    def __init__(
        self,
        in_fort10="fort.10_in",
        namelist=Namelist(),
    ):

        """
        input values
        """
        self.in_fort10 = in_fort10
        self.namelist = namelist

        file_check(in_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag:
            file_check("pseudo.dat")

    def __str__(self):

        output = [
            "TurboRVB convertfort10mol python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        pass

    def generate_input(self, input_name="convertfort10mol.input"):
        self.namelist.write(input_name)
        logger.info(f"{input_name} has been generated.")

    def run(self, input_name="convertfort10mol.input", output_name="out_mol"):
        run(
            turbo_convertfort10mol_run_command,
            input_name=input_name,
            output_name=output_name,
        )

    def check_results(self, output_names=["out_mol"]):
        flags = []
        for output_name in output_names:
            file_check(output_name)
            with open(output_name, "r") as f:
                lines = f.readlines()
            if any(
                [
                    re.match(r".*Time.*change.*fort\.10.*", line)
                    for line in lines
                ]
            ):
                flags.append(True)
            else:
                flags.append(False)
        return flags

    @staticmethod
    def read_default_namelist(in_fort10="fort.10_in"):
        convertfort10mol_default_file = os.path.join(
            pyturbo_data_dir, "convertfort10mol", "convertfort10mol.input"
        )
        namelist = Namelist.parse_namelist_from_file(
            convertfort10mol_default_file
        )
        io_fort10 = IO_fort10(in_fort10)
        nmol = int(io_fort10.f10header.nel / 2)
        namelist.set_parameter("nmol", nmol, "&molec_info")
        return namelist

    @staticmethod
    def read_namelist_from_file(file):
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist

    @classmethod
    def parse_from_default_namelist(cls, in_fort10="fort.10_in"):
        namelist = cls.read_default_namelist(in_fort10=in_fort10)
        return cls(in_fort10=in_fort10, namelist=namelist)

    @classmethod
    def parse_from_file(cls, file, in_fort10="fort.10_in"):
        namelist = Namelist.parse_namelist_from_file(file)
        return cls(in_fort10=in_fort10, namelist=namelist)


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
