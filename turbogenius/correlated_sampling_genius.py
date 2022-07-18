#!python
# -*- coding: utf-8 -*-

"""

Correlated sampling related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

#python modules
import os, sys
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import click
import pickle

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pyturbo.io_fort10 import IO_fort10
from pyturbo.vmc import VMC
from pyturbo.readforward import Readforward

#turbo-genius modules
from utils_workflows.env import turbo_genius_root
from geniusIO import GeniusIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)

from turbo_genius_cli import cli, decorate_grpost, header
@cli.command(short_help = "correlated_sampling_genius")
@decorate_grpost
@click.option("-steps", "vmcsteps",
              help= 'Specify vmcsteps',
              default = 1000,
              type = int)
@click.option("-bin", "bin_block",
              help= 'Specify bin_block',
              default = 2,
              type = int)
@click.option("-warmup", "warmupblocks",
              help= 'Specify warmupblocks',
              default = 1,
              type = int)
@click.option("-nw", "num_walkers",
              help= 'Specify num_walkers',
              default = -1,
              type = int)
@click.option("-maxtime", "maxtime",
              help= 'Specify maxtime',
              default = 3600,
              type = int)
@click.option("-twist", "twist_average",
              help= 'flag for twist_average',
              is_flag=True,
              type=bool
              )
@header
def correlated_sampling(
            g:bool,r:bool,post:bool,
            operation:bool,
            log_level:str,
            vmcsteps:int,
            bin_block:int,
            warmupblocks:int,
            num_walkers:int,
            maxtime:int,
            twist_average:bool
)->bool:
    pkl_name="correlated_sampling_genius_cli.pkl"
    root_dir=os.getcwd()
    pkl_file=os.path.join(root_dir, pkl_name)

    if g:
        readforward_genius=Correlated_sampling_genius(
            in_fort10='fort.10',
            corr_fort10="fort.10_corr",
            vmcsteps=vmcsteps,
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            num_walkers=num_walkers,
            maxtime=maxtime,
            twist_average=twist_average,
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

class Correlated_sampling_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo LRDMC class

    Attributes:
         in_fort10 (str): fort.10 WF file
         corr_fort10 (str): fort.10 WF file
         vmcsteps (int): total number of MCMC steps.
         bin_block (int): binning length
         warmupblocks (int): the number of disregarded blocks
         num_walkers (int): The number of walkers, -1 (default) = the number of MPI processes
         maxtime (int): Maxtime (sec.)
         twist_average (bool): Twist average flag, True or False
         kpoints (list): k Monkhorst-Pack grids, [kx,ky,kz,nx,ny,nz], kx,y,z-> grids, nx,y,z-> shift=0, noshift=1.
    """
    def __init__(self,
                 in_fort10:str="fort.10_in",
                 corr_fort10:str="fort.10_corr",
                 vmcsteps:int = 100,
                 bin_block:int = 10,
                 warmupblocks:int = 2,
                 num_walkers:int=-1,  # default -1 -> num of MPI process.
                 maxtime:int=172800,
                 twist_average:bool=False,
                 kpoints:list=[1,1,1,0,0,0],
                 ):

        self.in_fort10 = in_fort10
        self.corr_fort10 = corr_fort10

        self.vmcsteps = vmcsteps
        self.bin_block = bin_block
        self.warmupblocks = warmupblocks
        self.num_walkers = num_walkers
        self.maxtime = maxtime
        self.twist_average=twist_average
        self.kpoints=kpoints

        # VMC
        self.vmc=VMC.parse_from_default_namelist(in_fort10=in_fort10, twist_average=twist_average)
        if vmcsteps < 40 * bin_block + bin_block * warmupblocks:
            logger.warning(f"vmcsteps = {vmcsteps} is too small! < 40 * bin_block + bin_block * warmupblocks = {40 * bin_block + bin_block * warmupblocks}")
            logger.warning(f"vmcsteps = {vmcsteps} is set to 40 * bin_block + bin_block * warmupblocks = {40 * bin_block + bin_block * warmupblocks}")
            vmcsteps = 40 * bin_block + bin_block * warmupblocks
        self.vmc.set_parameter(parameter="ngen", value=vmcsteps, namelist="&simulation")
        self.vmc.set_parameter(parameter="maxtime", value=maxtime, namelist="&simulation")
        if num_walkers != -1: self.vmc.set_parameter(parameter="nw", value=num_walkers, namelist="&simulation")

        self.vmc.set_parameter(parameter="iread", value=3, namelist="&readio")

        # kpoints
        if self.twist_average: # not 0 (= not False)!!
            if self.twist_average == 1: # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.vmc.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                self.vmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.vmc.set_parameter(parameter="kp_type", value=1, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nkx", value=nkx, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nky", value=nky, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nkz", value=nkz, namelist="&kpoints")
                self.vmc.set_parameter(parameter="kx", value=kx, namelist="&kpoints")
                self.vmc.set_parameter(parameter="ky", value=ky, namelist="&kpoints")
                self.vmc.set_parameter(parameter="kz", value=kz, namelist="&kpoints")
                self.vmc.set_parameter(parameter="skip_equivalence", value='.true.', namelist="&kpoints")
                self.vmc.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
            elif self.twist_average == 2: # k-points are set from the user
                assert len(self.kpoints) == 2
                kpoints_up, kpoints_dn = self.kpoints
                assert len(kpoints_up) == len(kpoints_dn)
                for kup, kdn in zip(kpoints_up, kpoints_dn):
                    assert len(kup) == 4 # kx, ky, kz, wkp for up
                    assert len(kdn) == 4 # kx, ky, kz, wkp for dn
                self.vmc.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                self.vmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.vmc.set_parameter(parameter="kp_type", value=2, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nk1", value=len(kpoints_up)+len(kpoints_dn), namelist="&kpoints")
                self.vmc.manual_kpoints=self.kpoints

            else:
                logger.error(f"twist_average = {self.twist_average} is not implemented.")
                raise NotImeplementedError

        # readforward
        self.readforward=Readforward.parse_from_default_namelist(in_fort10=in_fort10)
        self.readforward.set_parameter(parameter="bin_length", value=bin_block, namelist="&corrfun")
        self.readforward.set_parameter(parameter="initial_bin", value=warmupblocks, namelist="&corrfun")
        self.readforward.set_parameter(parameter="correlated_samp", value='.true.', namelist="&corrfun")

    def run_all(self, input_name:str="datasvmc.input", vmc_output_name:str="out_vmc", readforward_output_name:str="out_readforward")->None:
        """
            Generate input files and run the command.

            Args:
                input_name (str): input file name
                vmc_output_name (str): vmc output file name
                readforward_output_name (str): readforward output file name

        """
        self.vmc.generate_input(input_name=input_name)
        self.vmc.run(input_name=input_name, output_name=vmc_output_name)
        self.readforward.generate_input(input_name=input_name)
        self.readforward.run(input_name=input_name, output_name=readforward_output_name)

    def generate_input(self, input_name:str="datasvmc.input")->None:
        """
            Generate input file.

            Args:
                input_name (str): input file name

        """
        self.vmc.generate_input(input_name=input_name)
        self.readforward.generate_input(input_name="readforward.input")

    def run(self, input_name:str="datasvmc.input", vmc_output_name:str="out_vmc", readforward_output_name:str="out_readforward")->None:
        """
            Run the command.

            Args:
                input_name (str): input file name
                vmc_output_name (str): vmc output file name
                readforward_output_name (str): readforward output file name
        """
        self.vmc.run(input_name=input_name, output_name=vmc_output_name)
        flags=self.vmc.check_results(output_names=[vmc_output_name])
        assert all(flags)
        self.readforward.run(input_name=input_name, output_name=readforward_output_name)
        flags=self.readforward.check_results(output_names=[readforward_output_name])
        assert all(flags)

    def check_results(self, vmc_output_names:list=["out_vmc"], readforward_output_names:list=["out_readforward"])->bool:
        """
            Check the result.

            Args:
                vmc_output_names (list): a list of output file names
                readforward_output_names (list): a list of output file names
            Returns:
                bool: True if all the runs were successful, False if an error is detected in the files.
        """
        return self.readforward.check_results(output_names=readforward_output_names) + self.vmc.check_results(output_names=vmc_output_names)

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
    