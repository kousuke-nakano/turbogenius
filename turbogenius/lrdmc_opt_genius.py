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
from pyturbo.lrdmcopt import LRDMCopt
from pyturbo.io_fort10 import IO_fort10

#turbo-genius modules
from utils_workflows.env import turbo_genius_root
from utils_workflows.utility import get_optimizer_flags, get_nonlocalmoves_setting
from geniusIO import GeniusIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)

from turbo_genius_cli import cli, decorate_grpost, header
@cli.command(short_help = "lrdmc_genius")
@decorate_grpost
@click.option("-lrdmcoptsteps", "lrdmcoptsteps",
              help= 'Specify lrdmcoptsteps',
              default = 1000,
              type = int)
@click.option("-optwarmup", "optwarmupsteps",
              help= 'Specify optwarmupsteps',
              default = 1,
              type = int)
@click.option("-steps", "steps",
              help= 'Specify steps per one iteration',
              default = 1,
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
              default = 1,
              type = int)
@click.option("-optimizer", "optimizer",
              help= 'Specify optimizer, sr or lr',
              default = "lr",
              type = str)
@click.option("-learn", "learning_rate",
              help= 'Specify learning_rate',
              default = 0.35,
              type = float)
@click.option("-reg", "regularization",
              help= 'Specify regularization',
              default = 0.001,
              type = float)
@click.option("-alat", "alat",
              help= 'Specify alat',
              default = -0.20,
              type = float)
@click.option("-etry", "etry",
              help= 'Specify etry',
              default = 0.0,
              type = float)
@click.option("-nonlocal", "nonlocalmoves",
              help= 'Specify nonlocalmoves, tmove, dla, dlatm',
              default = "tmove",
              type = str)
@click.option("-opt_onebody", "opt_onebody",
              help= 'flag for opt_onebody',
              is_flag=True,
              type = bool
              )
@click.option("-opt_twobody", "opt_twobody",
              help= 'flag for opt_twobody',
              is_flag=True,
              type=bool
              )
@click.option("-opt_det_mat", "opt_det_mat",
              help= 'flag for opt_det_mat',
              is_flag=False,
              type=bool
              )
@click.option("-opt_jas_mat", "opt_jas_mat",
              help= 'flag for opt_jas_mat',
              is_flag=True,
              type=bool
              )
@click.option("-opt_det_basis_exp", "opt_det_basis_exp",
              help= 'flag for opt_det_basis_exp',
              is_flag=False,
              type=bool
              )
@click.option("-opt_jas_basis_exp", "opt_jas_basis_exp",
              help= 'flag for opt_jas_basis_exp',
              is_flag=False,
              type=bool
              )
@click.option("-opt_det_basis_coeff", "opt_det_basis_coeff",
              help= 'flag for opt_det_basis_coeff',
              is_flag=False,
              type=bool
              )
@click.option("-opt_jas_basis_coeff", "opt_jas_basis_coeff",
              help= 'flag for opt_jas_basis_coeff',
              is_flag=False,
              type=bool
              )
@click.option("-twist", "twist_average",
              help= 'flag for twist_average',
              is_flag=False,
              type=bool
              )
@click.option("-kpts", "kpoints",
              help= 'kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]',
              nargs=6,
              default=[0,0,0,0,0,0],
              type=int
              )
@header
def lrdmcopt(
            g,r,post,
            operation,
            log_level,
            lrdmcoptsteps,
            optwarmupsteps,
            steps,
            bin_block,
            warmupblocks,
            num_walkers,  # default -1 -> num of MPI process.
            maxtime,
            optimizer,
            learning_rate,
            regularization,
            alat,
            etry,
            nonlocalmoves,  # tmove, dla, dlatm
            opt_onebody,
            opt_twobody,
            opt_det_mat,
            opt_jas_mat,
            opt_det_basis_exp,
            opt_jas_basis_exp,
            opt_det_basis_coeff,
            opt_jas_basis_coeff,
            twist_average,
            kpoints,
):
    if g:
        lrdmcopt_genius=LRDMCopt_genius(
            lrdmcoptsteps=lrdmcoptsteps,
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            steps=steps,
            num_walkers=num_walkers,  # default -1 -> num of MPI process.
            maxtime=maxtime,
            optimizer=optimizer,
            learning_rate=learning_rate,
            regularization=regularization,
            alat=alat,
            etry=etry,
            nonlocalmoves=nonlocalmoves,
            opt_onebody=opt_onebody,
            opt_twobody=opt_twobody,
            opt_det_mat=opt_det_mat,
            opt_jas_mat=opt_jas_mat,
            opt_det_basis_exp=opt_det_basis_exp,
            opt_jas_basis_exp=opt_jas_basis_exp,
            opt_det_basis_coeff=opt_det_basis_coeff,
            opt_jas_basis_coeff=opt_jas_basis_coeff,
            twist_average=twist_average,
            kpoints=kpoints,
        )
        lrdmcopt_genius.generate_input()

    if r:
        lrdmcopt_genius.run()

    if post:
        lrdmcopt.check_results()
        # avarege? optwarmupsteps=optwarmupsteps,

class LRDMCopt_genius(GeniusIO):

    def __init__(self,
                 fort10="fort.10",
                 lrdmcoptsteps = 100,
                 steps = 10,
                 bin_block = 1,
                 warmupblocks = 1,
                 num_walkers=-1,  # default -1 -> num of MPI process.
                 maxtime=172800,
                 optimizer="sr",
                 learning_rate=0.02,
                 regularization=0.001,
                 alat=-0.20,
                 etry=0.0,
                 nonlocalmoves="tmove",  # tmove, dla, dlatm
                 opt_onebody=True,
                 opt_twobody=True,
                 opt_det_mat=False,
                 opt_jas_mat=True,
                 opt_det_basis_exp=False,
                 opt_jas_basis_exp=False,
                 opt_det_basis_coeff=False,
                 opt_jas_basis_coeff=False,
                 twist_average=False,
                 kpoints=[1, 1, 1, 0, 0, 0],
                 ):

        self.fort10=fort10
        self.twist_average = twist_average

        self.energy = None
        self.energy_error = None
        self.estimated_time_for_1_generation = None

        optimizer_number, iesdonebodyoff, iesdtwobodyoff, twobodyoff, iesd, iesfree, iessw, iesup, iesm = \
            get_optimizer_flags(
                    optimizer=optimizer,
                    opt_onebody=opt_onebody,
                    opt_twobody=opt_twobody,
                    opt_det_mat=opt_det_mat,
                    opt_jas_mat=opt_jas_mat,
                    opt_det_basis_exp=opt_det_basis_exp,
                    opt_jas_basis_exp=opt_jas_basis_exp,
                    opt_det_basis_coeff=opt_det_basis_coeff,
                    opt_jas_basis_coeff=opt_jas_basis_coeff,
                    qmc_type="lrdmc"
            )

        self.lrdmcopt=LRDMCopt.parse_from_default_namelist(in_fort10=fort10, twist_average=twist_average)

        self.lrdmcopt.set_parameter(parameter="itestr4", value=optimizer_number, namelist="&simulation")
        self.lrdmcopt.set_parameter(parameter="ngen", value=lrdmcoptsteps * steps, namelist="&simulation")
        self.lrdmcopt.set_parameter(parameter="maxtime", value=maxtime, namelist="&simulation")
        if num_walkers != -1: self.lrdmcopt.set_parameter(parameter="nw", value=num_walkers, namelist="&simulation")
        self.lrdmcopt.set_parameter(parameter="nweight", value=steps, namelist="&optimization")
        self.lrdmcopt.set_parameter(parameter="nbinr", value=bin_block, namelist="&optimization")
        self.lrdmcopt.set_parameter(parameter="iboot", value=warmupblocks, namelist="&optimization")
        self.lrdmcopt.set_parameter(parameter="tpar", value=learning_rate, namelist="&optimization")
        self.lrdmcopt.set_parameter(parameter="parr", value=regularization, namelist="&optimization")

        self.lrdmcopt.set_parameter(parameter="iesdonebodyoff", value=iesdonebodyoff, namelist="&optimization")
        self.lrdmcopt.set_parameter(parameter="iesdtwobodyoff", value=iesdtwobodyoff, namelist="&optimization")
        self.lrdmcopt.set_parameter(parameter="twobodyoff", value=twobodyoff, namelist="&optimization")

        self.lrdmcopt.set_parameter(parameter="iesd", value=iesd, namelist="&parameters")
        self.lrdmcopt.set_parameter(parameter="iesfree", value=iesfree, namelist="&parameters")
        self.lrdmcopt.set_parameter(parameter="iessw", value=iessw, namelist="&parameters")
        self.lrdmcopt.set_parameter(parameter="iesup", value=iesup, namelist="&parameters")
        self.lrdmcopt.set_parameter(parameter="iesm", value=iesm, namelist="&parameters")

        self.lrdmcopt.set_parameter(parameter="alat", value=alat, namelist="&dmclrdmc")
        self.lrdmcopt.set_parameter(parameter="etry", value=etry, namelist="&dmclrdmc")

        typereg, npow = get_nonlocalmoves_setting(nonlocalmoves=nonlocalmoves)
        self.lrdmcopt.set_parameter(parameter="typereg", value=typereg, namelist="&dmclrdmc")
        self.lrdmcopt.set_parameter(parameter="npow", value=npow, namelist="&dmclrdmc")

        if self.twist_average:  # not 0 (= not False)!!
            logger.error("twist_average=true is not implemented.")
            raise NotImplemtedError

    def run_all(self, cont=False, input_name="datasfn_opt.input", output_name="out_fn_opt", average_parameters=True):
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        if average_parameters: self.average(input_name=input_name, output_name=output_name)

    def generate_input(self, cont=False, input_name="datasfn_opt.input"):
        io_fort10=IO_fort10(fort10=self.fort10)
        io_fort10.io_flag=0
        if cont: self.lrdmcopt.set_parameter("iopt", 0, "$systems")
        self.lrdmcopt.generate_input(input_name=input_name)

    def run(self, input_name="datasfn_opt.input", output_name="out_fn_opt"):
        self.lrdmcopt.run(input_name=input_name, output_name=output_name)
        flags=self.lrdmcopt.check_results(output_names=[output_name])
        assert all(flags)

    def store_result(self, output_names=["out_fn_opt"]):
        self.energy, self.energy_error = self.get_energy(output_names=output_names)
        self.estimated_time_for_1_generation = self.get_estimated_time_for_1_generation(output_names=output_names)

    def average(self, optwarmupsteps=10, graph_plot=False, input_name="datasfn_opt.input", output_name="out_fn_opt"):
        flags = self.lrdmcopt.check_results(output_names=[output_name])
        assert all(flags)
        self.lrdmcopt.average_optimized_parameters(equil_steps=optwarmupsteps, input_file_used=input_name, graph_plot=graph_plot)

    def get_energy(self, output_names=["out_fn_opt"]):
        return self.lrdmcopt.get_energy(output_names=output_names)

    def get_estimated_time_for_1_generation(self, output_names=["out_fn_opt"]):
        return self.lrdmcopt.get_estimated_time_for_1_generation(output_names=output_names)

    def check_results(self, output_names=["out_fn_opt"]):
        return self.lrdmcopt.check_results(output_names=output_names)

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("INFO")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "lrdmcopt"))

    lrdmcopt_genius=LRDMCopt_genius(
                fort10="fort.10",
                lrdmcoptsteps=100,
                optwarmupsteps=10,
                steps=500,
                bin_block=2,
                warmupblocks=1,
                optimizier="sr",
                learning_rate=0.02,
                regularization=0.001,
                alat=-0.20,
                etry=-1.0,
                nonlocalmoves="tmove",
                opt_onebody=True,
                opt_twobody=True,
                opt_det_mat=False,
                opt_jas_mat=False,
                opt_det_basis_exp=False,
                opt_jas_basis_exp=False,
                opt_det_basis_coeff=False,
                opt_jas_basis_coeff=False,
                twist_average=False,
                )
    lrdmcopt_genius.run()