#!python
# -*- coding: utf-8 -*-

"""

vmc genius related classes and methods

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
from pyturbo.vmc import VMC
from pyturbo.io_fort10 import IO_fort10

#turbo-genius modules
from utils_workflows.env import turbo_genius_root
from tools_genius import copy_jastrow
from geniusIO import GeniusIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)

from turbo_genius_cli import cli, decorate_grpost, header
@cli.command(short_help = "vmc_genius")
@decorate_grpost
@click.option("-steps", "vmcsteps",
              help= 'Specify vmcsteps',
              default = 1000,
              type = int)
@click.option("-bin", "bin_block",
              help= 'Specify bin_block',
              default = 1,
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
@click.option("-kpts", "kpoints",
              help= 'kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]',
              nargs=6,
              default=[0,0,0,0,0,0],
              type=int
              )
@click.option("-force", "force_calc_flag",
              help= 'flag for force_calc_flag',
              is_flag=True,
              type=bool
              )
@header
def vmc(
            g:bool,r:bool,post:bool,
            operation:bool,
            log_level:str,
            vmcsteps:int,
            bin_block:int,
            warmupblocks:int,
            num_walkers:int,
            maxtime:int,
            twist_average:bool,
            kpoints:list,
            force_calc_flag:bool
)->None:
    pkl_name="vmc_genius_cli.pkl"
    root_dir=os.getcwd()
    pkl_file=os.path.join(root_dir, pkl_name)

    if g:
        vmc_genius=VMC_genius(
            vmcsteps=vmcsteps,
            num_walkers=num_walkers,
            maxtime=maxtime,
            twist_average=twist_average,
            kpoints=kpoints,
            force_calc_flag=force_calc_flag
        )
        vmc_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(vmc_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                vmc_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        vmc_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                vmc_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags=vmc_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        vmc_genius.check_results()
        vmc_genius.compute_energy_and_forces(bin_block=bin_block, warmupblocks=warmupblocks)

class VMC_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo VMC class

    Attributes:
         fort10 (str): fort.10 WF file
         vmcsteps (int): total number of MCMC steps.
         num_walkers (int): The number of walkers, -1 (default) = the number of MPI processes
         maxtime (int): Maxtime (sec.)
         twist_average (bool): Twist average flag, True or False
         kpoints (list): k Monkhorst-Pack grids, [kx,ky,kz,nx,ny,nz], kx,y,z-> grids, nx,y,z-> shift=0, noshift=1.
         force_calc_flag (bool): if True, compute energy and force, if False, compute only energy
    """
    def __init__(self,
                 fort10:str="fort.10",
                 vmcsteps:int = 100,
                 num_walkers:int=-1,  # default -1 -> num of MPI process.
                 maxtime:int=172800,
                 twist_average:bool=False,
                 kpoints:list=[1,1,1,0,0,0],
                 force_calc_flag:bool=True
                 ):

        self.force_calc_flag = force_calc_flag
        self.vmcsteps = vmcsteps
        self.num_walkers = num_walkers
        self.maxtime = maxtime
        self.twist_average=twist_average
        self.kpoints=kpoints

        self.estimated_time_for_1_generation = None

        self.vmc=VMC.parse_from_default_namelist(in_fort10=fort10, twist_average=twist_average)

        self.vmc.set_parameter(parameter="ngen", value=vmcsteps, namelist="&simulation")
        self.vmc.set_parameter(parameter="maxtime", value=maxtime, namelist="&simulation")
        if num_walkers != -1: self.vmc.set_parameter(parameter="nw", value=num_walkers, namelist="&simulation")

        self.energy=None
        self.energy_error=None
        self.forces=None  #np.array([[]]) # 3 * natom matrix
        self.forces_error=None #np.array([[]])  # 3 * natom matrix

        #pseudo integration
        self.vmc.set_parameter(parameter="npsamax", value=4, namelist="&pseudo")

        # should be dicussed with Sandro.
        # Do you want to compute forces? if not, it is better to switch epscut=0.0 off.
        if not self.force_calc_flag:
           self.vmc.set_parameter(parameter="epscut", value=0.0, namelist="&vmc")
        else:
            self.vmc.comment_out(parameter="epscut")
            self.vmc.set_parameter(parameter="ieskin", value=1, namelist="&parameters")

        # kpoints
        if self.twist_average: # not 0 (= not False)!!
            if self.twist_average == 1: # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.vmc.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                #self.vmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.vmc.set_parameter(parameter="kp_type", value=1, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nk1", value=nkx, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nk2", value=nky, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nk3", value=nkz, namelist="&kpoints")
                self.vmc.set_parameter(parameter="k1", value=kx, namelist="&kpoints")
                self.vmc.set_parameter(parameter="k2", value=ky, namelist="&kpoints")
                self.vmc.set_parameter(parameter="k3", value=kz, namelist="&kpoints")
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
                #self.vmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.vmc.set_parameter(parameter="kp_type", value=2, namelist="&kpoints")
                self.vmc.set_parameter(parameter="nk1", value=len(kpoints_up), namelist="&kpoints")
                self.vmc.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
                self.vmc.manual_kpoints=self.kpoints

            else:
                logger.error(f"twist_average = {self.twist_average} is not implemented.")
                raise NotImeplementedError


    def run_all(self, cont:bool=False, input_name:str="datasvmc.input", output_name:str="out_vmc")->None:
        """
            Generate input files and run the command.

            Args:
                input_name (str): input file name
                output_name (str): output file name

        """
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        self.compute_energy_and_forces()

    def generate_input(self, cont:bool=False, input_name:str="datasvmc.input")->None:
        """
            Generate input file.

            Args:
                cont (bool): if True, continuation run (i.e., iopt=0), if False, starting from scratch (i.e., iopt=1).
                input_name (str): input file name

        """
        if cont: self.vmc.set_parameter("iopt", 0, "&simulation")
        self.vmc.generate_input(input_name=input_name)

    def run(self, input_name="datasvmc.input", output_name="out_vmc"):
        """
            Run the command.

            Args:
                input_name (str): input file name
                output_name (str): output file name
        """
        self.vmc.run(input_name=input_name, output_name=output_name)
        flags=self.vmc.check_results(output_names=[output_name])
        assert all(flags)

    def store_result(self, bin_block:int=10, warmupblocks:int=5, output_names:list=["out_vmc"], rerun:bool=False)->bool:
        """
            Store results. This procedure stores estimated_time_for_1_generation, energy, and energy_error.
            This method is needed for storing data and access to them later.

            Args:
                bin_block (int): binning length
                warmupblocks (int): the number of disregarded blocks
                output_names (list): a list of output file names
                rerun (bool): if true, compute energy and force again even if there are energy and force files.
        """
        self.estimated_time_for_1_generation = self.get_estimated_time_for_1_generation(output_names=output_names)
        self.energy, self.energy_error = self.vmc.get_energy(init=warmupblocks, bin=bin_block, rerun=rerun)

    def compute_energy_and_forces(self, bin_block:int=10, warmupblocks:int=5, rerun:bool=False)->None:
        """
            Compute energy and forces

            Args:
                bin_block (int): binning length
                warmupblocks (int): the number of disregarded blocks
                rerun (bool): if true, compute energy and force again even if there are energy and force files.
        """
        self.energy, self.energy_error = self.vmc.get_energy(init=warmupblocks, bin=bin_block, rerun=rerun)
        if self.force_calc_flag:
            self.forces, self.forces_error = self.vmc.get_forces(init=warmupblocks, bin=bin_block, rerun=False)

    def get_estimated_time_for_1_generation(self, output_names:list=["out_vmc"])->float:
        """
            This procedure stores estimated_time_for_1_generation.

            Args:
                output_names (list): a list of output file names

            Return:
                float: estimated_time_for_1_generation.
        """
        return self.vmc.get_estimated_time_for_1_generation(output_names=output_names)

    def check_results(self, output_names:list=["out_vmc"])->bool:
        """
            Check the result.

            Args:
                output_names (list): a list of output file names
            Return:
                bool: True if all the runs were successful, False if an error is detected in the files.
        """
        return self.vmc.check_results(output_names=output_names)

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "vmc"))

    # moved to examples