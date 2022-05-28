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
from pyturbo.io_fort10 import IO_fort10
from pyturbo.readforward import Readforward

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
@click.option("-corr", "corr_sampling",
              help= 'correlated sampling',
              is_flag=True,
              default=True)
@click.option("-bin", "bin_block",
              help= 'Specify bin_block',
              default = 1,
              type = int)
@click.option("-warmup", "warmupblocks",
              help= 'Specify warmupblocks',
              default = 1,
              type = int)
@header
def readforward(
            g,r,post,
            operation,
            log_level,
            bin_block,
            warmupblocks,
            corr_sampling,
):
    pkl_name="readforward_genius_cli.pkl"
    root_dir=os.getcwd()
    pkl_file=os.path.join(root_dir, pkl_name)

    if g:
        readforward_genius=Readforward_genius(
            in_fort10='fort.10',
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            corr_sampling=corr_sampling
        )
        readforward_genius.generate_input()
        with open(pkl_file, "wb") as f:
            pickle.dump(readforward_genius, f)
    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                readforward_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        readforward_genius.run()

    if post:
        if post:
            os.chdir(root_dir)
            try:
                with open(pkl_file, "rb") as f:
                    readforward_genius = pickle.load(f)
            except FileNotFoundError:
                logger.error("Did you generate your input file using turbogenius?")
                raise FileNotFoundError
            flags = readforward_genius.check_results()
            if all(flags):
                logger.info("Job was successful.")
            else:
                logger.info("Job was failure. See the output file.")
                return
            readforward_genius.check_results()

class Readforward_genius(GeniusIO):

    def __init__(self,
                 in_fort10="fort.10_in",
                 corr_fort10="fort.10_corr",
                 bin_block=10,
                 warmupblocks=2,
                 corr_sampling=True
                 ):

        self.in_fort10 = in_fort10
        self.corr_fort10 = corr_fort10

        self.readforward=Readforward.parse_from_default_namelist(in_fort10=in_fort10)
        self.io_fort10 = IO_fort10(self.in_fort10)

        self.readforward.set_parameter(parameter="bin_length", value=bin_block, namelist="&corrfun")
        self.readforward.set_parameter(parameter="initial_bin", value=warmupblocks, namelist="&corrfun")
        if corr_sampling:
            self.readforward.set_parameter(parameter="correlated_samp", value='.true.', namelist="&corrfun")
        else:
            pass
            #self.convertfort10.comment_out(parameter="nz")

    def run_all(self, input_name="datasvmc.input", output_name="out_readforward"):
        self.generate_input(input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)

    def generate_input(self, input_name="datasvmc.input"):
        self.readforward.generate_input(input_name=input_name)

    def run(self, input_name="datasvmc.input", output_name="out_readforward"):
        self.readforward.run(input_name=input_name, output_name=output_name)
        flags=self.readforward.check_results(output_names=[output_name])
        assert all(flags)

    def check_results(self, output_names=["out_readforward"]):
        return self.readforward.check_results(output_names=output_names)

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
    