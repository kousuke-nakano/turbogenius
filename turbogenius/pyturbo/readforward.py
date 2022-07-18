#!python
# -*- coding: utf-8 -*-

"""

pyturbo: readforward related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

#python modules
import os, sys
import re

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from namelist import Namelist
from fortranIO import FortranIO
from utils.env import turbo_readforward_run_command
from utils.env import pyturbo_data_dir, pyturbo_root
from utils.utility import file_check
from io_fort10 import IO_fort10
from utils.execute import run

from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo').getChild(__name__)

class Readforward(FortranIO):

    def __init__(self,
                 in_fort10='fort.10',
                 namelist=Namelist(),
                 ):

        """
        input values
        """
        file_check(in_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag: file_check("pseudo.dat")
        self.in_fort10=in_fort10
        self.namelist=namelist

    def __str__(self):

        output = [
            "TurboRVB readforward python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        pass

    def generate_input(self, input_name="readforward.input"):
        self.namelist.write(input_name)
        logger.info(f"{input_name} has been generated.")

    def run(self, input_name="datasvmc.input", output_name="out_readforward"):
        run(turbo_readforward_run_command, input_name=input_name, output_name=output_name)

    def check_results(self, output_names=["out_readforward"]):
        flags = []
        for output_name in output_names:
            file_check(output_name)
            with open(output_name, "r") as f:
                lines = f.readlines()
            if any([re.match(r".*total.*#.*bin.*considered.*", line) for line in lines]):
                flags.append(True)
            else:
                flags.append(False)
        return flags

    @staticmethod
    def read_default_namelist(in_fort10="fort.10"):
        readforward_default_file = os.path.join(pyturbo_data_dir, "readforward", "readforward.input")
        namelist = Namelist.parse_namelist_from_file(readforward_default_file)
        return namelist
    @staticmethod
    def read_namelist_from_file(file):
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist
    @classmethod
    def parse_from_default_namelist(cls, in_fort10="fort.10"):
        namelist = cls.read_default_namelist(in_fort10=in_fort10)
        return cls(
            in_fort10=in_fort10,
            namelist=namelist
        )
    @classmethod
    def parse_from_file(cls, file, in_fort10="fort.10"):
        namelist = Namelist.parse_namelist_from_file(file)
        return cls(
            in_fort10=in_fort10,
            namelist=namelist
        )

if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples.

    readforward_test_dir=os.path.join(pyturbo_root, "tests", "readforward")
    os.chdir(readforward_test_dir)

    readforward=Readforward().parse_from_default_namelist()
    readforward.generate_input()