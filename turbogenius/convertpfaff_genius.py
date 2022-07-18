#!python
# -*- coding: utf-8 -*-

"""

convertpfaff genius related classes and methods


"""

#python modules
import os, sys
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import click

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pyturbo.convertpfaff import Convertpfaff
from pyturbo.io_fort10 import IO_fort10

#turbo-genius modules
from utils_workflows.env import turbo_genius_root
from geniusIO import GeniusIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)

from turbo_genius_cli import cli, decorate_grpost, header
@cli.command(short_help = "readforward_genius")
@decorate_grpost
@click.option("-rotate", "rotate_flag",
              help= 'rotate_flag',
              is_flag=True,
              default=False)
@click.option("-angle", "rotate_angle",
              help= 'Specify rotate_angle',
              default = 0.0,
              type = float)
@click.option("-scale", "scale_mean_field",
              help= 'Specify scale_mean_field',
              default = 1000,
              type = int)
@header
def convertpfaff(
            g:bool,r:bool,post:bool,
            operation:bool,
            log_level:str,
            rotate_flag:bool,
            rotate_angle:float,
            scale_mean_field:int
):
    pkl_name="convertpfaff_genius_cli.pkl"
    root_dir=os.getcwd()
    pkl_file=os.path.join(root_dir, pkl_name)

    if g:
        convertpfaff_genius=Convertpfaff_genius(
            in_fort10='fort.10_in',
            out_fort10='fort.10_out'
        )
        convertpfaff_genius.generate_input()
        with open(pkl_file, "wb") as f:
            pickle.dump(convertpfaff_genius, f)
    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                convertpfaff_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        convertpfaff_genius.run(rotate_flag=rotate_flag, rotate_angle=rotate_angle, scale_mean_field=scale_mean_field)

    if post:
        if post:
            os.chdir(root_dir)
            try:
                with open(pkl_file, "rb") as f:
                    convertpfaff_genius = pickle.load(f)
            except FileNotFoundError:
                logger.error("Did you generate your input file using turbogenius?")
                raise FileNotFoundError
            flags = convertpfaff_genius.check_results()
            if all(flags):
                logger.info("Job was successful.")
            else:
                logger.info("Job was failure. See the output file.")
                return
            convertpfaff_genius.check_results()

class Convertpfaff_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo Convertpfaff class

    Attributes:
         in_fort10 (str): fort.10 WF file (input WF file)
         out_fort10 (str): fort.10 WF file (template WF file)
    """
    def __init__(self,
                 in_fort10:str="fort.10_in",
                 out_fort10:str="fort.10_out",
                 ):

        self.in_fort10 = in_fort10
        self.out_fort10 = out_fort10

        self.convertpfaff=Convertpfaff.parse_from_default_namelist(in_fort10=in_fort10, out_fort10=out_fort10)

    def run_all(self, rotate_flag:bool=False, rotate_angle:float=0, scale_mean_field:int=1000, output_name:str="out_pfaff")->None:
        """
            Generate input files and run the command.

            Args:
                rotate_flag (bool): rotation flag, True or False
                rotate_angle (float): rotation angle
                scale_mean_field (int): scaling mean field
                output_name (str): output file name

        """
        self.generate_input()
        self.run(rotate_flag=rotate_flag, rotate_angle=rotate_angle, scale_mean_field=scale_mean_field, output_name=output_name)

    def run(self, rotate_flag:float=False, rotate_angle:float=0, scale_mean_field:int=1000, output_name:str="out_pfaff")->None:
        """
            Run the command.

            Args:
                rotate_flag (bool): rotation flag, True or False
                rotate_angle (float): rotation angle
                scale_mean_field (int): scaling mean field
                output_name (str): output file name
        """
        self.convertpfaff.run(rotate_flag=rotate_flag, rotate_angle=rotate_angle, scale_mean_field=scale_mean_field, output_name=output_name)
        flags=self.convertpfaff.check_results(output_names=[output_name])
        assert all(flags)

    def check_results(self, output_names:list=["out_pfaff"])->bool:
        """
            Check the result.

            Args:
                output_names (list): a list of output file names
            Returns:
                bool: True if all the runs were successful, False if an error is detected in the files.
        """
        return self.convertpfaff.check_results(output_names=output_names)

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "convertfort10"))

    # moved to examples
    convertfort10_genius=Convertfort10_genius()
    convertfort10_genius.generate_input()
    convertfort10_genius.run()
    