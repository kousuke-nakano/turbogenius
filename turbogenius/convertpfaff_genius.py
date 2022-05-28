#!python
# -*- coding: utf-8 -*-

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
            g,r,post,
            operation,
            log_level,
            rotate_flag,
            rotate_angle,
            scale_mean_field
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

    def __init__(self,
                 in_fort10="fort.10_in",
                 out_fort10="fort.10_out",
                 ):

        self.in_fort10 = in_fort10
        self.out_fort10 = out_fort10

        self.convertpfaff=Convertpfaff.parse_from_default_namelist(in_fort10=in_fort10, out_fort10=out_fort10)

    def run_all(self, rotate_flag=False, rotate_angle=0, scale_mean_field=1000, output_name="out_pfaff"):
        self.generate_input()
        self.run(rotate_flag=rotate_flag, rotate_angle=rotate_angle, scale_mean_field=scale_mean_field, output_name=output_name)

    def generate_input(self):
        pass

    def run(self, rotate_flag=False, rotate_angle=0, scale_mean_field=1000, output_name="out_pfaff"):
        self.convertpfaff.run(rotate_flag=rotate_flag, rotate_angle=rotate_angle, scale_mean_field=scale_mean_field, output_name=output_name)
        flags=self.convertpfaff.check_results(output_names=[output_name])
        assert all(flags)

    def check_results(self, output_names=["out_pfaff"]):
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
    