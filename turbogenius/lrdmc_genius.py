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
from pyturbo.lrdmc import LRDMC
from pyturbo.lrdmc_extrapolation import LRDMC_extrapolation
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
            g,r,post,
            operation,
            log_level,
            lrdmcsteps,
            bin_block,
            warmupblocks,
            correcting_factor,
            alat,
            etry,
            num_walkers,  # default -1 -> num of MPI process.
            maxtime,
            twist_average,
            force_calc_flag,
            nonlocalmoves # tmove, dla, dlatm
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

    def __init__(self,
                 fort10="fort.10",
                 lrdmcsteps = 100,
                 alat=-0.20,
                 etry=0.0,
                 num_walkers=-1,  # default -1 -> num of MPI process.
                 maxtime=172800,
                 twist_average=False,
                 force_calc_flag=False,
                 nonlocalmoves="tmove" # tmove, dla, dlatm
                 ):

        self.force_calc_flag=force_calc_flag
        self.twist_average = twist_average

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

        if self.twist_average:  # not 0 (= not False)!!
            logger.error("twist_average=true is not implemented.")
            raise NotImplemtedError

    def run_all(self, bin_block=10, warmupblocks=2, correcting_factor=2, cont=False, input_name="datasfn.input", output_name="out_fn"):
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        self.compute_energy_and_forces(bin_block=bin_block, warmupblocks=warmupblocks, correcting_factor=correcting_factor)

    def generate_input(self, cont=False, input_name="datasfn.input"):
        if cont: self.lrdmc.set_parameter("iopt", 0, "$systems")
        self.lrdmc.generate_input(input_name=input_name)

    def run(self, input_name="datasfn.input", output_name="out_fn"):
        self.lrdmc.run(input_name=input_name, output_name=output_name)
        flags=self.lrdmc.check_results(output_names=[output_name])
        assert all(flags)

    def store_result(self, bin_block=10, warmupblocks=2, correcting_factor=2, output_names=["out_fn"]):
        self.estimated_time_for_1_generation = self.get_estimated_time_for_1_generation(output_names=output_names)
        self.energy, self.energy_error = self.lrdmc.get_energy(init=warmupblocks, correct=correcting_factor, bin=bin_block)

    def compute_energy_and_forces(self, bin_block=10, warmupblocks=2, correcting_factor=2):
        self.energy, self.energy_error = self.lrdmc.get_energy(init=warmupblocks, correct=correcting_factor, bin=bin_block)

    def get_estimated_time_for_1_generation(self, output_names=["out_fn"]):
        return self.lrdmc.get_estimated_time_for_1_generation(output_names=output_names)

    def check_results(self, output_names=["out_fn"]):
        return self.lrdmc.check_results(output_names=output_names)

class LRDMC_ext_genius(GeniusIO):

    def __init__(self,
                 fort10="fort.10",
                 lrdmcsteps = 100,
                 maxtime=172800,
                 alat_list=[],
                 etry=0.0,
                 num_walkers=-1,  # default -1 -> num of MPI process.
                 polynominal_order=2,
                 twist_average=False,
                 force_calc_flag=False,
                 nonlocalmoves="tmove"  # tmove, dla, dlatm
                 ):

        self.force_calc_flag=force_calc_flag

        self.energy=None
        self.energy_error=None
        self.forces=None #np.array([[]]) # 3 * natom matrix
        self.forces_error=None #np.array([[]])  # 3 * natom matrix

        self.polynominal_order = polynominal_order

        self.lrdmc_ext=LRDMC_extrapolation.parse_from_default_namelist(in_fort10=fort10,
                                                                       alat_list=alat_list,
                                                                       etry=etry,
                                                                       twist_average=twist_average)
        self.lrdmc_ext.set_parameter(parameter="ngen", value=lrdmcsteps, namelist="&simulation")
        self.lrdmc_ext.set_parameter(parameter="maxtime", value=maxtime, namelist="&simulation")
        if num_walkers != -1: self.lrdmc_ext.set_parameter(parameter="nw", value=num_walkers, namelist="&simulation")

        typereg, npow = get_nonlocalmoves_setting(nonlocalmoves=nonlocalmoves)
        self.lrdmc_ext.set_parameter(parameter="typereg", value=typereg, namelist="&dmclrdmc")
        self.lrdmc_ext.set_parameter(parameter="npow", value=npow, namelist="&dmclrdmc")

    def run_all(self, bin_block=10, warmupblocks=2, correcting_factor=2, cont=False, input_name="datasfn.input", output_name="out_fn"):
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        self.do_extrapolation(bin_block=bin_block, warmupblocks=warmupblocks, correcting_factor=correcting_factor, input_name=input_name, output_name=output_name)

    def generate_input(self, cont=False, input_name="datasfn.input"):
        if cont: self.lrdmc_ext.set_parameter("iopt", 0, "$systems")
        self.lrdmc_ext.generate_input(input_name=input_name)

    def run(self, input_name="datasfn.input", output_name="out_fn"):
        self.lrdmc_ext.run(input_name=input_name, output_name=output_name)
        flags=self.lrdmc_ext.check_results(output_name=output_name)
        assert all(flags)

    def do_extrapolation(self, bin_block=10, warmupblocks=2, correcting_factor=2, input_name="datasfn.input", output_name="out_fn"):
        flags=self.lrdmc_ext.check_results(output_name=output_name)
        assert all(flags)
        e,err=self.lrdmc_ext.get_extrapolated_energy(degree_poly=self.polynominal_order,
                                               input_name=input_name, graph_plot=True,
                                               init=warmupblocks, correct=correcting_factor,
                                               bin=bin_block, num_proc=1)
        self.energy=e; self.energy_error=err
        if self.force_calc_flag: raise NotImplementedError

    def check_results(self, output_names=["out_fn"]):
        self.lrdmc_ext.check_results(output_names=output_names)

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