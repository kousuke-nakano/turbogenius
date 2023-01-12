#!python
# -*- coding: utf-8 -*-

"""

pyturbo: lrdmcopt related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# python modules
import os, sys
import shutil
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# turbo-genius modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from namelist import Namelist
from fortranIO import FortranIO
from io_fort10 import IO_fort10
from utils.env import pyturbo_data_dir, turborvb_bin_root
from utils.env import turbo_qmc_run_command
from utils.utility import file_check, file_check_flag, remove_file
from utils.utility import remove_new_parameter_lines_in_fort10
from utils.execute import run

# Logger
from logging import config, getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)
# logger = getLogger(__name__)


class LRDMCopt(FortranIO):
    def __init__(
        self,
        in_fort10="fort.10",
        namelist=Namelist(),
        pp_flag=False,
        twist_average=False,
    ):

        """
        input values
        """
        file_check(in_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag:
            file_check("pseudo.dat")

        self.in_fort10 = in_fort10
        self.namelist = namelist
        self.pp_flag = pp_flag
        self.twist_average = twist_average

    def __str__(self):

        output = [
            "TurboRVB lrdmcopt python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        pass

    def generate_input(self, input_name):
        self.namelist.write(input_name)
        logger.info(f"{input_name} has been generated. \n")

    def run(self, input_name="datasfn_opt.input", output_name="out_fn_opt"):
        run(
            turbo_qmc_run_command,
            input_name=input_name,
            output_name=output_name,
        )
        remove_file(file="stop.dat")

    def check_results(self, output_names=["out_fn_opt"]):
        flags = []
        for output_name in output_names:
            file_check(output_name)
            with open(output_name, "r") as f:
                lines = f.readlines()
            if any(
                [re.match(r".*TurboRVB.*profiling.*", line) for line in lines]
            ):
                flags.append(True)
            else:
                flags.append(False)

        return flags

    def plot_energy_and_devmax(
        self, output_names=["out_fn_opt"], interactive=False
    ):

        # plot the energies and devmax
        logger.info("Plotting the energies and devmax")
        out_min = []
        for output_name in output_names:
            with open(output_name, "r") as f:
                out_min += f.readlines()

        col = ["dum1", "dum2", "dum3", "Energy", "error"]
        energy_list = list(
            map(
                lambda x: x.split(),
                [i for i in out_min if re.match(r".*New.*Energy.*", i)],
            )
        )
        energy_pandas_str = pd.DataFrame(energy_list, columns=col)

        col = ["dum1", "dum2", "dum3", "dum4", "devmax", "num", "num2"]
        devmax_list = list(
            map(
                lambda x: x.split(),
                [
                    i
                    for i in out_min
                    if re.match(r".*devmax.*par.*Normal.*", i)
                ],
            )
        )
        devmax_pandas_str = pd.DataFrame(devmax_list, columns=col)

        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.direction"] = "in"
        plt.rcParams["xtick.major.width"] = 1.0
        plt.rcParams["ytick.major.width"] = 1.0
        plt.rcParams["font.size"] = 12
        plt.rcParams["axes.linewidth"] = 1.5
        fig, ax1 = plt.subplots()
        ax1.errorbar(
            energy_pandas_str.index,
            energy_pandas_str["Energy"].astype("float"),
            energy_pandas_str["error"].astype("float"),
            color="black",
            marker="o",
            linestyle="dashed",
            capsize=2,
            label="Energy",
        )

        ax2 = ax1.twinx()
        ax2.plot(
            devmax_pandas_str.index,
            devmax_pandas_str["devmax"].astype("float"),
            color="red",
            marker="x",
            linestyle="dashed",
            label="devmax",
        )
        ax2.hlines(
            4.5,
            devmax_pandas_str.index.min(),
            devmax_pandas_str.index.max(),
            "blue",
            linestyles="dashed",
            label="devmax < 4.5 (converged criteria)",
        )

        plt.xlabel("Optimiztion")
        ax1.set_ylabel("Energy (Ha)")
        ax2.set_ylabel("devmax")

        plt.gca().get_yaxis().get_major_formatter().set_useOffset(
            False
        )  # No offset for y-axis
        plt.gca().get_xaxis().set_major_locator(
            ticker.MaxNLocator(integer=True)
        )  # Interger for x-axis

        h1, l1 = ax1.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax1.legend(h1 + h2, l1 + l2, loc="upper right")
        if interactive:
            plt.waitforbuttonpress()
        plt.savefig(
            "plot_energy_and_devmax.png", bbox_inches="tight", pad_inches=0.2
        )
        plt.close()

    def get_energy(self, output_names=["out_fn_opt"]):

        out_min = []
        for output_name in output_names:
            with open(output_name, "r") as f:
                out_min += f.readlines()

        energy_list = list(
            map(
                float,
                [
                    i.split()[3]
                    for i in out_min
                    if re.match(r".*New.*Energy.*", i)
                ],
            )
        )
        error_list = list(
            map(
                float,
                [
                    i.split()[4]
                    for i in out_min
                    if re.match(r".*New.*Energy.*", i)
                ],
            )
        )

        return energy_list, error_list

    def get_devmax(self, output_names=["out_fn_opt"]):

        out_min = []
        for output_name in output_names:
            with open(output_name, "r") as f:
                out_min += f.readlines()

        devmax_list = list(
            map(
                float,
                [
                    i.split()[4]
                    for i in out_min
                    if re.match(r".*devmax.*par.*Normal.*", i)
                ],
            )
        )

        return devmax_list

    def get_estimated_time_for_1_generation(self, output_names=["out_fn_opt"]):

        out_min = []
        for output_name in output_names:
            with open(output_name, "r") as f:
                out_min += f.readlines()

        ave_time_1000_generations = list(
            map(
                float,
                [
                    i.split()[5]
                    for i in out_min
                    if re.match(
                        r".*Average.*time.*for.*1000.*generations.*", i
                    )
                ],
            )
        )
        ave_time_1_generation = (
            np.mean(np.array(ave_time_1000_generations)) / 1000
        )

        return ave_time_1_generation  # sec.

    def average_optimized_parameters(
        self,
        equil_steps=10,
        input_file_used="datasfnopt.input",
        graph_plot=False,
    ):
        if self.twist_average:
            raise NotImplementedError

        current_dir = os.getcwd()

        # averaging
        with open("forces.dat", "r") as f:
            lines = f.readlines()
        length_iter = len(lines)
        logger.info(f"The total number of iterations is {length_iter}")
        logger.info(f"Equil. steps {equil_steps}")
        assert length_iter > equil_steps  # this assertion is more general!!

        # remove previous parameters before average.
        logger.info("Removing previous averaged parameters")
        remove_new_parameter_lines_in_fort10()  # workaround!! to be refactored in the near future.

        # save parameters
        if graph_plot:
            cmd = f"(echo '1 1 0 0'; echo '0'; echo '100000') | {os.path.join(turborvb_bin_root, 'readalles.x')}"
            run(binary=cmd, output_name="out_readalles_for_plot_all")
            shutil.copyfile("./story.d", "./all_story.d")
            all_story_pandas_data = pd.read_csv(
                "./all_story.d",
                header=None,
                delim_whitespace=True,
                index_col=False,
            )

        # average and show the result
        logger.info("Averaging parameters using readalles.x")
        logger.info(
            f"The first {equil_steps} iterations are disregarded in the average."
        )
        cmd = f"(echo '1 {equil_steps + 1} 1 0'; echo '0'; echo '100000') | {os.path.join(turborvb_bin_root, 'readalles.x')}"
        run(binary=cmd, output_name="out_readalles_for_average")

        if graph_plot:
            shutil.copyfile("./story.d", "./average_story.d")
            ave_story_pandas_data = pd.read_csv(
                "./average_story.d",
                header=None,
                delim_whitespace=True,
                index_col=False,
            )

        if graph_plot:
            graph_save_dir = os.path.join(current_dir, "parameters_graphs")
            os.makedirs(graph_save_dir, exist_ok=True)

            for i in range(1, len(all_story_pandas_data.columns)):
                plt.figure(figsize=(10, 8))
                plt.plot(
                    all_story_pandas_data[0].astype("int"),
                    all_story_pandas_data[i].astype("float"),
                    color="black",
                    marker="o",
                    linestyle="dashed",
                    label="all",
                )
                plt.plot(
                    ave_story_pandas_data[0].astype("int"),
                    ave_story_pandas_data[i].astype("float"),
                    color="red",
                    marker="o",
                    linestyle="dashed",
                    label="averaged",
                )
                plt.xlabel(
                    "Iteration", fontname="Times New Roman", fontsize=14
                )
                plt.ylabel("Value", fontname="Times New Roman", fontsize=14)
                # plt.legend(frameon=True)
                plt.title("Parameter_No.{}".format(i))
                plt.legend(frameon=True)
                plt.gca().get_yaxis().get_major_formatter().set_useOffset(
                    False
                )  # No offset for y-axis
                plt.gca().get_xaxis().set_major_locator(
                    ticker.MaxNLocator(integer=True)
                )  # Interger for x-axis
                plt.savefig(
                    os.path.join(
                        graph_save_dir, "Parameter_No{}.png".format(i)
                    ),
                    bbox_inches="tight",
                    pad_inches=0.2,
                )
                plt.close()

        # update fort10 using the averaged parameters
        logger.info("Update fort10 using the averaged parameters")
        logger.info(f"The input file used for the opt. is {input_file_used}")
        lrdmcopt = LRDMCopt.parse_from_file(
            file=input_file_used,
            in_fort10=self.in_fort10,
            twist_average=self.twist_average,
        )
        lrdmcopt.set_parameter("iopt", 1)
        lrdmcopt.set_parameter("ngen", 0)

        # replace iopt from 0 -> 1
        io_fort10 = IO_fort10(self.in_fort10)
        io_fort10.f10header.io_flag = 1

        dir_ave_temp = os.path.join(current_dir, "ave_temp")
        os.makedirs(dir_ave_temp, exist_ok=True)

        shutil.copy(
            os.path.join(os.getcwd(), "fort.10"),
            os.path.join(dir_ave_temp, "fort.10"),
        )
        lrdmcopt.generate_input(os.path.join(dir_ave_temp, "ave.input"))

        # to be refactored with some flag
        if file_check_flag("pseudo.dat"):
            logger.info("calc. with pseudo potentials")
            shutil.copyfile(
                os.path.join(os.getcwd(), "pseudo.dat"),
                os.path.join(dir_ave_temp, "pseudo.dat"),
            )

        os.chdir(dir_ave_temp)
        cmd = f"{os.path.join(turborvb_bin_root, 'turborvb-serial.x')}"
        run(binary=cmd, input_name="ave.input", output_name="out_ave")

        os.chdir(current_dir)
        shutil.move(
            os.path.join(current_dir, "fort.10"),
            os.path.join(current_dir, "fort.10_bak"),
        )
        shutil.copyfile(
            os.path.join(dir_ave_temp, "fort.10"),
            os.path.join(current_dir, "fort.10"),
        )
        logger.info("The averaged fort10 is labeled as fort.10")
        logger.info("The original fort10 is saved as fort.10_bak")

    @staticmethod
    def read_default_namelist():
        lrdmcopt_default_file = os.path.join(
            pyturbo_data_dir, "lrdmcopt", "datasfn_opt.input"
        )
        namelist = Namelist.parse_namelist_from_file(lrdmcopt_default_file)
        return namelist

    @staticmethod
    def read_namelist_from_file(file):
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist

    @classmethod
    def parse_from_default_namelist(
        cls, in_fort10="fort.10", twist_average=False
    ):
        namelist = cls.read_default_namelist()
        return cls(
            in_fort10=in_fort10, namelist=namelist, twist_average=twist_average
        )

    @classmethod
    def parse_from_file(cls, file, in_fort10="fort.10", twist_average=False):
        namelist = Namelist.parse_namelist_from_file(file)
        return cls(
            in_fort10=in_fort10, namelist=namelist, twist_average=twist_average
        )


if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples
