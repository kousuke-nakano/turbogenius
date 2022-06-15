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
from pyturbo.vmcopt import VMCopt
from pyturbo.io_fort10 import IO_fort10

#turbo-genius modules
from utils_workflows.env import turbo_genius_root
from utils_workflows.utility import get_optimizer_flags
from tools_genius import copy_jastrow
from geniusIO import GeniusIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)

from turbo_genius_cli import cli, decorate_grpost, header
@cli.command(short_help = "vmcopt_genius")
@decorate_grpost
@click.option("-vmcoptsteps", "vmcoptsteps",
              help= 'Specify vmcoptsteps',
              default = 1000,
              type = int)
@click.option("-optwarmup", "optwarmupsteps",
              help= 'Specify optwarmupsteps',
              default = 1,
              type = int)
@click.option("-steps", "steps",
              help= 'Specify steps per one iteration',
              default = 20,
              type = int)
@click.option("-bin", "bin_block",
              help= 'Specify bin_block',
              default = 1,
              type = int)
@click.option("-warmup", "warmupblocks",
              help= 'Specify warmupblocks',
              default = 0,
              type = int)
@click.option("-nw", "num_walkers",
              help= 'Specify num_walkers',
              default = -1,
              type = int)
@click.option("-maxtime", "maxtime",
              help= 'Specify maxtime',
              default = 3600,
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
@click.option("-opt_onebody", "opt_onebody",
              help= 'flag for opt_onebody',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-opt_twobody", "opt_twobody",
              help= 'flag for opt_twobody',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-opt_det_mat", "opt_det_mat",
              help= 'flag for opt_det_mat',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-opt_jas_mat", "opt_jas_mat",
              help= 'flag for opt_jas_mat',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-opt_det_basis_exp", "opt_det_basis_exp",
              help= 'flag for opt_det_basis_exp',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-opt_jas_basis_exp", "opt_jas_basis_exp",
              help= 'flag for opt_jas_basis_exp',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-opt_det_basis_coeff", "opt_det_basis_coeff",
              help= 'flag for opt_det_basis_coeff',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-opt_jas_basis_coeff", "opt_jas_basis_coeff",
              help= 'flag for opt_jas_basis_coeff',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-twist", "twist_average",
              help= 'flag for twist_average',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-kpts", "kpoints",
              help= 'kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]',
              nargs=6,
              default=[0,0,0,0,0,0],
              type=int
              )
@click.option("-plot", "plot_graph",
              help= 'flag for plotting graph',
              is_flag=True,
              default=False,
              type=bool
              )
@click.option("-interactive", "plot_interactive",
              help= 'flag for interactive plotting graph',
              is_flag=True,
              default=False,
              type=bool
              )
@header
def vmcopt(
            g,r,post,
            operation,
            log_level,
            vmcoptsteps=100,
            optwarmupsteps=10,
            steps=10,
            bin_block=1,
            warmupblocks=0,
            num_walkers=-1,  # default -1 -> num of MPI process.
            maxtime=172800,
            optimizer="sr",
            learning_rate=0.35,
            regularization=0.001,
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
            plot_graph=False,
            plot_interactive=False
):
    pkl_name="vmcopt_genius_cli.pkl"
    root_dir=os.getcwd()
    pkl_file=os.path.join(root_dir, pkl_name)

    if g:
        vmcopt_genius=VMCopt_genius(
            vmcoptsteps=vmcoptsteps,
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            steps=steps,
            num_walkers=num_walkers,  # default -1 -> num of MPI process.
            maxtime=maxtime,
            optimizer=optimizer,
            learning_rate=learning_rate,
            regularization=regularization,
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
        vmcopt_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(vmcopt_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                vmcopt_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        vmcopt_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                vmcopt_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags=vmcopt_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        if plot_graph: vmcopt_genius.plot_energy_and_devmax(interactive=plot_interactive)
        if plot_graph: vmcopt_genius.plot_parameters_history(interactive=plot_interactive)
        vmcopt_genius.average(optwarmupsteps=optwarmupsteps, graph_plot=plot_graph)

class VMCopt_genius(GeniusIO):

    def __init__(self,
                 fort10="fort.10",
                 vmcoptsteps = 100,
                 steps = 10,
                 bin_block = 1,
                 warmupblocks = 0,
                 num_walkers=-1, # default -1 -> num of MPI process.
                 maxtime=172800,
                 optimizer="sr",
                 learning_rate=0.35,
                 regularization=0.001,
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

        self.fort10 = fort10
        self.twist_average = twist_average
        self.kpoints = kpoints

        self.optimizer = optimizer
        self.opt_onebody = opt_onebody
        self.opt_twobody = opt_twobody
        self.opt_det_mat = opt_det_mat
        self.opt_jas_mat = opt_jas_mat
        self.opt_det_basis_exp = opt_det_basis_exp
        self.opt_jas_basis_exp = opt_jas_basis_exp
        self.opt_det_basis_coeff = opt_det_basis_coeff
        self.opt_jas_basis_coeff = opt_jas_basis_coeff

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
                    qmc_type="vmc"
            )

        self.energy=None
        self.energy_error=None
        self.estimated_time_for_1_generation=None

        self.vmcopt=VMCopt.parse_from_default_namelist(in_fort10=fort10, twist_average=self.twist_average)

        self.vmcopt.set_parameter(parameter="itestr4", value=optimizer_number, namelist="&simulation")
        self.vmcopt.set_parameter(parameter="ngen", value=vmcoptsteps*steps, namelist="&simulation")
        self.vmcopt.set_parameter(parameter="maxtime", value=maxtime, namelist="&simulation")
        if num_walkers!=-1: self.vmcopt.set_parameter(parameter="nw", value=num_walkers, namelist="&simulation")
        self.vmcopt.set_parameter(parameter="nweight", value=steps, namelist="&optimization")
        self.vmcopt.set_parameter(parameter="nbinr", value=bin_block, namelist="&optimization")
        self.vmcopt.set_parameter(parameter="iboot", value=warmupblocks, namelist="&optimization")
        self.vmcopt.set_parameter(parameter="tpar", value=learning_rate, namelist="&optimization")
        self.vmcopt.set_parameter(parameter="parr", value=regularization, namelist="&optimization")

        self.vmcopt.set_parameter(parameter="iesdonebodyoff", value=iesdonebodyoff, namelist="&optimization")
        self.vmcopt.set_parameter(parameter="iesdtwobodyoff", value=iesdtwobodyoff, namelist="&optimization")
        self.vmcopt.set_parameter(parameter="twobodyoff", value=twobodyoff, namelist="&optimization")

        self.vmcopt.set_parameter(parameter="iesd", value=iesd, namelist="&parameters")
        self.vmcopt.set_parameter(parameter="iesfree", value=iesfree, namelist="&parameters")
        self.vmcopt.set_parameter(parameter="iessw", value=iessw, namelist="&parameters")
        self.vmcopt.set_parameter(parameter="iesup", value=iesup, namelist="&parameters")
        self.vmcopt.set_parameter(parameter="iesm", value=iesm, namelist="&parameters")

        # Does the VMC optimization changes the nodal surface? if not, it is better to switch off epscut option.
        if opt_det_mat or opt_det_basis_exp or opt_det_basis_coeff:
            self.vmcopt.comment_out(parameter="epscut")
        else:
            self.vmcopt.set_parameter(parameter="epscut", value=0.0, namelist="&vmc")

        # kpoints
        if self.twist_average: # not 0 (= not False)!!
            if self.twist_average == 1: # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.vmcopt.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                #self.vmcopt.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.vmcopt.set_parameter(parameter="kp_type", value=1, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="nk1", value=nkx, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="nk2", value=nky, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="nk3", value=nkz, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="k1", value=kx, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="k2", value=ky, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="k3", value=kz, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="skip_equivalence", value='.true.', namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
            elif self.twist_average == 2: # k-points are set from the user
                assert len(self.kpoints) == 2
                kpoints_up, kpoints_dn = self.kpoints
                assert len(kpoints_up) == len(kpoints_dn)
                for kup, kdn in zip(kpoints_up, kpoints_dn):
                    assert len(kup) == 4 # kx, ky, kz, wkp for up
                    assert len(kdn) == 4 # kx, ky, kz, wkp for dn
                self.vmcopt.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                #self.vmcopt.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.vmcopt.set_parameter(parameter="kp_type", value=2, namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="nk1", value=len(kpoints_up), namelist="&kpoints")
                self.vmcopt.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
                self.vmcopt.manual_kpoints=self.kpoints

            else:
                logger.error(f"twist_average = {self.twist_average} is not implemented.")
                raise NotImeplementedError

    def run_all(self, cont=False, input_name="datasmin.input", output_name="out_min", average_parameters=True):
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        if average_parameters: self.average(input_name=input_name, output_name=output_name)

    def generate_input(self, cont=False, input_name="datasmin.input"):
        io_fort10=IO_fort10(fort10=self.fort10)
        io_fort10.io_flag=0
        if cont: self.vmcopt.set_parameter("iopt", 0, "$systems")
        self.vmcopt.generate_input(input_name=input_name)

    def run(self, input_name="datasmin.input", output_name="out_min"):
        self.vmcopt.run(input_name=input_name, output_name=output_name)
        flags=self.vmcopt.check_results(output_names=[output_name])
        assert all(flags)

    def check_results(self, output_names=["out_min"]):
        return self.vmcopt.check_results(output_names=output_names)

    def plot_energy_and_devmax(self, output_names=["out_min"], interactive=True):
        return self.vmcopt.plot_energy_and_devmax(output_names=output_names, interactive=interactive)

    def store_result(self, output_names=["out_min"]):
        self.energy, self.energy_error = self.get_energy(output_names=output_names)
        self.estimated_time_for_1_generation = self.get_estimated_time_for_1_generation(output_names=output_names)

    def average(self, optwarmupsteps, input_name="datasmin.input", output_name="out_min", graph_plot=False):
        current_dir=os.getcwd()
        if self.twist_average:
            if self.opt_det_mat or self.opt_det_basis_exp or self.opt_det_basis_coeff:
                logger.warning("The twist average with a JAGP ansatz, turbogenius assumes real_agp=.true. option.")
                twist_average_copyjas=False
            elif self.opt_jas_basis_exp or self.opt_jas_basis_coeff:
                logger.warning("Sorry opt_jas_basis_exp or opt_jas_basis_coeff are not supported.")
                raise NotImplementedError
            else:  # self.opt_jas_mat:
                logger.info("The twist average with a JDFT ansatz (i.e., optimize only jas. mat.) is now supported.")
                twist_average_copyjas = True
        else:
            logger.info("Twist-average flag is False. Open boundary condition or PBC at single k point.")
            twist_average_copyjas = False

        flags = self.vmcopt.check_results(output_names=[output_name])
        assert all(flags)
        self.vmcopt.average_optimized_parameters(equil_steps=optwarmupsteps, input_file_used=input_name, graph_plot=graph_plot)

        if twist_average_copyjas:
            logger.info("Additional commands are needed for averaging Jas. mat. with k average")
            logger.info("cp fort.10 turborvb.scratch/fort.10;")
            logger.info("cp fort.10 turborvb.scratch/fort.10_new;")
            logger.info("cp kp_info.dat turborvb.scratch/kp_info.dat;")
            logger.info("cp parminimized.d turborvb.scratch/parminimized.d;")
            logger.info("cd turborvb.scratch/")
            logger.info("copyjas.x kpoints")

            shutil.copy("fort.10", os.path.join("turborvb.scratch", "fort.10"))
            shutil.copy("fort.10", os.path.join("turborvb.scratch", "fort.10_new"))
            shutil.copy("kp_info.dat", os.path.join("turborvb.scratch", "kp_info.dat"))
            shutil.copy("parminimized.d", os.path.join("turborvb.scratch", "parminimized.d"))
            os.chdir("turborvb.scratch")
            copy_jastrow(args="kpoints")
            os.chdir(current_dir)

    def get_energy(self, output_names=["out_min"]):
        return self.vmcopt.get_energy(output_names=output_names)

    def get_estimated_time_for_1_generation(self, output_names=["out_min"]):
        return self.vmcopt.get_estimated_time_for_1_generation(output_names=output_names)

    def plot_parameters_history(self, interactive=True):
        self.vmcopt.plot_parameters_history(interactive=interactive)

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "vmcopt"))
    # removed to examples
