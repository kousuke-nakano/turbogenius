#!python
# -*- coding: utf-8 -*-

"""

lrdmc genius related classes and methods

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
from pyturbo.lrdmc import LRDMC
from pyturbo.io_fort10 import IO_fort10

#turbo-genius modules
from utils_workflows.env import turbo_genius_root
from utils_workflows.utility import get_nonlocalmoves_setting
from geniusIO import GeniusIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)


from turbo_genius_cli import cli, decorate_grpost, header
@cli.command(short_help = "lrdmc_genius")
@decorate_grpost
@click.option("-steps", "lrdmcsteps",
              help= 'Specify lrdmcsteps',
              default = 1000,
              type = int)
@click.option("-bin", "bin_block",
              help= 'Specify bin_block',
              default = 1,
              type = int)
@click.option("-corr", "correcting_factor",
              help= 'Specify correcting_factor',
              default = 1,
              type = int)
@click.option("-alat", "alat",
              help= 'Specify alat',
              default = -0.20,
              type = float)
@click.option("-etry", "etry",
              help= 'Specify etry',
              default = 0.0,
              type = float)
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
@click.option("-force", "force_calc_flag",
              help= 'flag for force_calc_flag',
              is_flag=True,
              type=bool
              )
@click.option("-nonlocal", "nonlocalmoves",
              help= 'Specify nonlocalmoves, tmove, dla, dlatm',
              default = "tmove",
              type = str)
@header
def lrdmc(
            g:bool,r:bool,post:bool,
            operation:bool,
            log_level:str,
            lrdmcsteps:int,
            bin_block:int,
            warmupblocks:int,
            correcting_factor:int,
            alat:float,
            etry:float,
            num_walkers:int,  # default -1 -> num of MPI process.
            maxtime:int,
            twist_average:bool,
            force_calc_flag:bool,
            nonlocalmoves:str # tmove, dla, dlatm
):
    pkl_name="lrdmc_genius_cli.pkl"
    root_dir=os.getcwd()
    pkl_file=os.path.join(root_dir, pkl_name)

    if g:
        lrdmc_genius=LRDMC_genius(
            lrdmcsteps=lrdmcsteps,
            alat=alat,
            etry=etry,
            num_walkers=num_walkers,  # default -1 -> num of MPI process.
            maxtime=maxtime,
            twist_average=twist_average,
            force_calc_flag=force_calc_flag,
            nonlocalmoves=nonlocalmoves  # tmove, dla, dlatm
        )
        lrdmc_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(lrdmc_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                lrdmc_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        lrdmc_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                lrdmc_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags=lrdmc_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        lrdmc_genius.check_results()
        lrdmc_genius.compute_energy_and_forces(bin_block=bin_block, warmupblocks=warmupblocks, correcting_factor=correcting_factor)

class LRDMC_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo LRDMC class

    Attributes:
         fort10 (str): fort.10 WF file
         lrdmcsteps (int): total number of MCMC steps.
         alat (float): Lattice space (Bohr)
         etry (float): Trial Energy (Ha)
         num_walkers (int): The number of walkers, -1 (default) = the number of MPI processes
         maxtime (int): Maxtime (sec.)
         twist_average (bool): Twist average flag, True or False
         kpoints (list): k Monkhorst-Pack grids, [kx,ky,kz,nx,ny,nz], kx,y,z-> grids, nx,y,z-> shift=0, noshift=1.
         force_calc_flag (bool): if True, compute energy and force, if False, compute only energy
         nonlocalmoves (str): Treatment of locality approximation, choose from "tmove", "dla", "dlatm"
    """
    def __init__(self,
                 fort10:str="fort.10",
                 lrdmcsteps:int = 100,
                 alat:float = -0.20,
                 etry:float = 0.0,
                 num_walkers:int = -1,  # default -1 -> num of MPI process.
                 maxtime:int = 172800,
                 twist_average:bool = False,
                 kpoints:list = [1, 1, 1, 0, 0, 0],
                 force_calc_flag:bool = False,
                 nonlocalmoves:str = "tmove" # tmove, dla, dlatm
                 ):

        self.force_calc_flag=force_calc_flag
        self.twist_average = twist_average
        self.kpoints=kpoints

        self.estimated_time_for_1_generation=None

        self.energy=None
        self.energy_error=None
        self.forces=None #np.array([[]]) # 3 * natom matrix
        self.forces_error=None #np.array([[]])  # 3 * natom matrix

        self.lrdmc=LRDMC.parse_from_default_namelist(in_fort10=fort10, twist_average=twist_average)
        self.lrdmc.set_parameter(parameter="ngen", value=lrdmcsteps, namelist="&simulation")
        self.lrdmc.set_parameter(parameter="maxtime", value=maxtime, namelist="&simulation")
        if num_walkers != -1: self.lrdmc.set_parameter(parameter="nw", value=num_walkers, namelist="&simulation")

        self.lrdmc.set_parameter(parameter="etry", value=etry, namelist="&dmclrdmc")
        self.lrdmc.set_parameter(parameter="alat", value=alat, namelist="&dmclrdmc")

        typereg, npow = get_nonlocalmoves_setting(nonlocalmoves=nonlocalmoves)
        self.lrdmc.set_parameter(parameter="typereg", value=typereg, namelist="&dmclrdmc")
        self.lrdmc.set_parameter(parameter="npow", value=npow, namelist="&dmclrdmc")

        # kpoints
        if self.twist_average: # not 0 (= not False)!!
            if self.twist_average == 1: # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.lrdmc.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                #self.lrdmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.lrdmc.set_parameter(parameter="kp_type", value=1, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="nk1", value=nkx, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="nk2", value=nky, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="nk3", value=nkz, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="k1", value=kx, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="k2", value=ky, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="k3", value=kz, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="skip_equivalence", value='.true.', namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
            elif self.twist_average == 2: # k-points are set from the user
                assert len(self.kpoints) == 2
                kpoints_up, kpoints_dn = self.kpoints
                assert len(kpoints_up) == len(kpoints_dn)
                for kup, kdn in zip(kpoints_up, kpoints_dn):
                    assert len(kup) == 4 # kx, ky, kz, wkp for up
                    assert len(kdn) == 4 # kx, ky, kz, wkp for dn
                self.lrdmc.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                #self.lrdmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.lrdmc.set_parameter(parameter="kp_type", value=2, namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="nk1", value=len(kpoints_up), namelist="&kpoints")
                self.lrdmc.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
                self.lrdmc.manual_kpoints=self.kpoints

            else:
                logger.error(f"twist_average = {self.twist_average} is not implemented.")
                raise NotImeplementedError

    def run_all(self, bin_block:int=10, warmupblocks:int=2, correcting_factor:int=2, cont:bool=False, input_name:str="datasfn.input", output_name:str="out_fn")->None:
        """
            Generate input files and run the command.

            Args:
                bin_block (int): binning length
                warmupblocks (int): the number of disregarded blocks,
                correcting_factor (int): correcting factors
                cont (bool): if True, continuation run (i.e., iopt=0), if False, starting from scratch (i.e., iopt=1).
                input_name (str): input file name
                output_name (str): output file name

        """
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        self.compute_energy_and_forces(bin_block=bin_block, warmupblocks=warmupblocks, correcting_factor=correcting_factor)

    def generate_input(self, cont:bool=False, input_name:str="datasfn.input")->None:
        """
            Generate input file.

            Args:
                cont (bool): if True, continuation run (i.e., iopt=0), if False, starting from scratch (i.e., iopt=1).
                input_name (str): input file name

        """
        if cont: self.lrdmc.set_parameter("iopt", 0, "&simulation")
        self.lrdmc.generate_input(input_name=input_name)

    def run(self, input_name:str="datasfn.input", output_name:str="out_fn")->None:
        """
            Run the command.

            Args:
                input_name (str): input file name
                output_name (str): output file name
        """
        self.lrdmc.run(input_name=input_name, output_name=output_name)
        flags=self.lrdmc.check_results(output_names=[output_name])
        assert all(flags)

    def store_result(self, bin_block:int=10, warmupblocks:int=2, correcting_factor:int=2, output_names:list=["out_fn"], rerun:bool=False)->None:
        """
            Store results. This procedure stores estimated_time_for_1_generation, energy, and energy_error.
            This method is needed for storing data and access to them later.

            Args:
                bin_block (int): binning length
                warmupblocks (int): the number of disregarded blocks
                correcting_factor (int): correcting factors
                output_names (list): a list of output file names
                rerun (bool): if true, compute energy and force again even if there are energy and force files.
        """
        self.estimated_time_for_1_generation = self.get_estimated_time_for_1_generation(output_names=output_names)
        self.energy, self.energy_error = self.lrdmc.get_energy(init=warmupblocks, correct=correcting_factor, bin=bin_block, rerun=rerun)

    def compute_energy_and_forces(self, bin_block:int=10, warmupblocks:int=2, correcting_factor:int=2, rerun:bool=False)->None:
        """
            Compute energy and forces

            Args:
                bin_block (int): binning length
                warmupblocks (int): the number of disregarded blocks
                correcting_factor (int): correcting factors
                rerun (bool): if true, compute energy and force again even if there are energy and force files.
        """
        self.energy, self.energy_error = self.lrdmc.get_energy(init=warmupblocks, correct=correcting_factor, bin=bin_block, rerun=rerun)

    def get_estimated_time_for_1_generation(self, output_names:list=["out_fn"])->float:
        """
            This procedure stores estimated_time_for_1_generation.

            Args:
                output_names (list): a list of output file names

            Return:
                float: estimated_time_for_1_generation.
        """
        return self.lrdmc.get_estimated_time_for_1_generation(output_names=output_names)

    def check_results(self, output_names:list=["out_fn"])->bool:
        """
            Check the result.

            Args:
                output_names (list): a list of output file names
            Return:
                bool: True if all the runs were successful, False if an error is detected in the files.
        """
        return self.lrdmc.check_results(output_names=output_names)

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("INFO")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "lrdmc"))

    # moved to examples