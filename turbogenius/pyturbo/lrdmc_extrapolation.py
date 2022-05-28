#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from io_fort10 import IO_fort10
from lrdmc import LRDMC
from utils.env import turbo_qmc_run_command
from utils.env import pyturbo_data_dir
from utils.utility import file_check, pygrep_lineno, get_line_from_file
from utils.execute import run
from namelist import Namelist
from fortranIO import FortranIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo').getChild(__name__)
#logger = getLogger(__name__)

class LRDMC_extrapolation(FortranIO):

    def __init__(self,
                 fort10="fort.10",
                 namelist=Namelist(),
                 alat_list=[0.10, 0.20, 0.30],
                 etry=0.0,
                 twist_average=False
                 ):


        self.fort10=fort10
        self.alat_list=alat_list
        self.etry=etry
        self.twist_average=twist_average

        file_check(fort10)
        if IO_fort10(fort10=fort10).pp_flag: file_check("pseudo.dat")

        self.namelist=namelist

        self.lrdmc = LRDMC(
                 in_fort10=self.fort10,
                 namelist=self.namelist,
                 twist_average=self.twist_average
                 )
        self.lrdmc.set_parameter("etry", etry, namelist="&dmclrdmc")

    def sanity_check(self):
        pass

    def generate_input(self, input_name="datasfn.input"):
        # check files exits
        file_check(self.fort10)
        io_fort10=IO_fort10(fort10=self.fort10)
        if io_fort10.pp_flag: file_check("pseudo.dat")

        # make directory and write fort.10 and input
        lrdmc_root_dir = os.getcwd()
        os.chdir(lrdmc_root_dir)
        for alat in self.alat_list:
            dir_alat = os.path.join(lrdmc_root_dir, f"alat_{alat}")
            os.makedirs(dir_alat, exist_ok=True)
            os.chdir(dir_alat)

            shutil.copyfile(os.path.join(lrdmc_root_dir, self.fort10), os.path.join(dir_alat, self.fort10))
            if io_fort10.pp_flag:
                shutil.copyfile(os.path.join(lrdmc_root_dir, "pseudo.dat"), os.path.join(dir_alat, "pseudo.dat"))

            self.lrdmc.set_parameter("alat", alat, "&dmclrdmc")

            self.lrdmc.generate_input(input_name=input_name)
            logger.info(f"datasfn.input file for alat = {alat} is generated.")

            os.chdir(lrdmc_root_dir)

    def run(self, alat_list=None, input_name="datasfn.input", output_name="out_fn"):
        lrdmc_root_dir = os.getcwd()

        if alat_list is None:
            alat_list=self.alat_list

        for alat in alat_list:
            dir_alat = os.path.join(lrdmc_root_dir, f"alat_{alat}")
            os.chdir(dir_alat)
            run(turbo_qmc_run_command, input_name=input_name, output_name=output_name)
            os.chdir(lrdmc_root_dir)

    def check_results(self, alat_list=None, output_name="out_fn"):
        lrdmc_root_dir = os.getcwd()
        flags=[]

        if alat_list is None:
            alat_list=self.alat_list

        for alat in alat_list:
            dir_alat = os.path.join(lrdmc_root_dir, f"alat_{alat}")
            os.chdir(dir_alat)
            flags+=self.lrdmc.check_results(output_names=[output_name])
            os.chdir(lrdmc_root_dir)

        return flags

    def get_extrapolated_energy(self, alat_list=None, degree_poly=2, input_name="datasfn.input", graph_plot=True, init=10, correct=10, bin=10, num_proc=1):
        lrdmc_root_dir = os.getcwd()

        if alat_list is None:
            alat_list=self.alat_list

        logger.info(f"the polynomial degree for fitting of energies with respect to alat^2 is {degree_poly}")

        evsa_line = f"{degree_poly}  {len(alat_list)}  4  1\n"
        evsa_gnu_line = f"# alat  energy  error\n"
        lrdmc_root_dir = os.getcwd()

        for alat in alat_list:
            dir_alat = os.path.join(lrdmc_root_dir, f"alat_{alat}")
            os.chdir(dir_alat)
            logger.debug(dir_alat)
            lrdmc=LRDMC.parse_from_file(file=input_name)
            logger.info(f"Postprocessing for alat = {alat}")
            energy, error = lrdmc.get_energy(init=init, correct=correct, bin=bin, num_proc=num_proc)
            evsa_line = evsa_line + f"{np.abs(alat)} {energy} {error}\n"
            evsa_gnu_line = evsa_gnu_line + f"{np.abs(alat)} {energy} {error}\n"
            os.chdir(lrdmc_root_dir)

        with open("evsa.in", 'w') as f:
            f.writelines(evsa_line)
        with open("evsa.gnu", 'w') as f:
            f.writelines(evsa_gnu_line)

        run(binary="funvsa.x", input_name="evsa.in", output_name="evsa.out")
        logger.info(f"Output from fitting process saved to evsa.out.")

        coeff_index = pygrep_lineno('evsa.out', 'Coefficient found')
        coeff_list = []
        coeff_error_list = []
        for poly in range(degree_poly + 1):
            coeff_list.append(float(get_line_from_file('evsa.out', coeff_index + 1 + poly).split()[1]))
            coeff_error_list.append(float(get_line_from_file('evsa.out', coeff_index + 1 + poly).split()[2]))
        coeff_list.reverse()  # because const -> x -> x2 -> ... in evsa.out
        coeff_error_list.reverse()
        energy_ext=coeff_list[-1]
        error_ext=coeff_error_list[-1]
        logger.info(f"Extrapolated lrdmc energy (a->0) = {energy_ext} Ha +- {error_ext} Ha")
        poly = np.poly1d(coeff_list)

        if graph_plot:
            with open("evsa.gnu", "r") as f:
                lines = f.readlines()
                alat = []
                energy = []
                energy_error = []
                for line in lines:
                    if '#' in line:
                        pass
                    else:
                        alat.append(float(line.split()[0]))
                        energy.append(float(line.split()[1]))
                        energy_error.append(float(line.split()[2]))

            alat_squared = np.array(alat) ** 2
            alat_squared_extrapolated = np.linspace(0, np.max(alat_squared) * 1.10, 500)
            energy = np.array(energy)
            energy_error = np.array(energy_error)

            plt.xlim([alat_squared_extrapolated.min(), alat_squared_extrapolated.max()])
            plt.annotate("E(alat->0) = {:.5f} Ha +- {:.5f} Ha".format(coeff_list[-1], coeff_error_list[-1]), xy=(0.05, 0.05),
                         xycoords="axes fraction")
            plt.errorbar(alat_squared, energy, yerr=energy_error, color="black", marker="o", fmt='')
            plt.plot(alat_squared_extrapolated, poly(alat_squared_extrapolated), color="blue", linestyle="dashed")
            plt.xlabel("alat$^2$ (Bohr$^2$)", fontname="Times New Roman", fontsize=14)
            plt.ylabel("Energy (Ha)", fontname="Times New Roman", fontsize=14)
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)  # No offset for y-axis
            plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)  # No offset for x-axis
            plt.savefig("Energy_vs_alat.png", bbox_inches='tight', pad_inches=0.2)
            logger.info("The graph of the extrapolation is saved as Energy_vs_alat.png")

        return energy_ext, error_ext

    @staticmethod
    def read_default_namelist():
        lrdmc_default_file = os.path.join(pyturbo_data_dir, "lrdmc", "datasfn.input")
        namelist = Namelist.parse_namelist_from_file(lrdmc_default_file)
        return namelist
    @staticmethod
    def read_namelist_from_file(file):
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist
    @classmethod
    def parse_from_default_namelist(cls, in_fort10="fort.10", alat_list=[], etry=0.0, twist_average=False):
        namelist = cls.read_default_namelist()
        return cls(
            fort10=in_fort10,
            namelist=namelist,
            alat_list=alat_list,
            etry=etry,
            twist_average=twist_average
        )
    @classmethod
    def parse_from_file(cls, file, in_fort10="fort.10", alat_list=[], etry=0.0, twist_average=False):
        namelist = Namelist.parse_namelist_from_file(file)
        return cls(
            fort10=in_fort10,
            namelist=namelist,
            alat_list=alat_list,
            etry=etry,
            twist_average=twist_average
        )

if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples