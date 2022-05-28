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
            g,r,post,
            operation,
            log_level,
            vmcsteps,
            bin_block,
            warmupblocks,
            num_walkers,
            maxtime,
            twist_average,
            kpoints,
            force_calc_flag
):
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

    def __init__(self,
                 fort10="fort.10",
                 vmcsteps = 100,
                 num_walkers=-1,  # default -1 -> num of MPI process.
                 maxtime=172800,
                 twist_average=False,
                 kpoints=[1,1,1,0,0,0],
                 force_calc_flag=True
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


    def run_all(self, cont=False, input_name="datasvmc.input", output_name="out_vmc"):
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        self.compute_energy_and_forces()

    def generate_input(self, cont=False, input_name="datasvmc.input"):
        if cont: self.vmc.set_parameter("iopt", 0, "$systems")
        self.vmc.generate_input(input_name=input_name)

    def run(self, input_name="datasvmc.input", output_name="out_vmc"):
        self.vmc.run(input_name=input_name, output_name=output_name)
        flags=self.vmc.check_results(output_names=[output_name])
        assert all(flags)

    def store_result(self, bin_block=10, warmupblocks=5, output_names=["out_vmc"]):
        self.estimated_time_for_1_generation = self.get_estimated_time_for_1_generation(output_names=output_names)
        self.energy, self.energy_error = self.vmc.get_energy(init=warmupblocks, bin=bin_block)

    def compute_energy_and_forces(self, bin_block=10, warmupblocks=5):
        self.energy, self.energy_error = self.vmc.get_energy(init=warmupblocks, bin=bin_block)
        if self.force_calc_flag:
            self.forces, self.forces_error = self.vmc.get_forces(init=warmupblocks, bin=bin_block)

    def get_estimated_time_for_1_generation(self, output_names=["out_vmc"]):
        return self.vmc.get_estimated_time_for_1_generation(output_names=output_names)

    def check_results(self, output_names=["out_vmc"]):
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