#!python
# -*- coding: utf-8 -*-

"""

pyturbo: vmcopt related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# python modules
import os
import shutil
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from typing import Optional

# Logger
from logging import getLogger, StreamHandler, Formatter


# turbo-genius modules
from turbogenius.pyturbo.namelist import Namelist
from turbogenius.pyturbo.fortranIO import FortranIO
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.pyturbo.utils.env import pyturbo_data_dir
from turbogenius.pyturbo.utils.env import (
    turbo_qmc_run_command,
    turborvb_bin_root,
)
from turbogenius.pyturbo.utils.utility import (
    file_check,
    file_check_flag,
    remove_new_parameter_lines_in_fort10,
    remove_file,
)
from turbogenius.pyturbo.utils.execute import run


logger = getLogger("pyturbo").getChild(__name__)


class VMCopt(FortranIO):
    def __init__(
        self,
        in_fort10: str = "fort.10",
        namelist: Optional[Namelist] = None,
        twist_average: bool = False,  # False or 0: single-k, True or 1: Monkhorst-Pack, 2: manual k-grid
    ):
        if namelist is None:
            namelist = Namelist()

        """
        input values
        """
        self.in_fort10 = in_fort10
        self.namelist = namelist
        self.twist_average = twist_average

        file_check(in_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag:
            file_check("pseudo.dat")

        # manual k-grid! [[[kx, ky, kz, wkp for up], ....], [[# kx, ky, kz, wkp for dn], ...]]
        self.__manual_kpoints = []

    @property
    def manual_kpoints(self):
        return self.__manual_kpoints

    @manual_kpoints.setter
    def manual_kpoints(self, kpoints):
        assert len(kpoints) == 2
        kpoints_up, kpoints_dn = kpoints
        assert len(kpoints_up) == len(kpoints_dn)
        for kup, kdn in zip(kpoints_up, kpoints_dn):
            assert len(kup) == 4  # kx, ky, kz, wkp for up
            assert len(kdn) == 4  # kx, ky, kz, wkp for dn
        self.__manual_kpoints = kpoints
        self.namelist.set_parameter(
            parameter="yes_kpoints", value=".true.", namelist="&parameters"
        )
        # self.namelist.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
        self.namelist.set_parameter(
            parameter="kp_type", value=2, namelist="&kpoints"
        )
        self.namelist.set_parameter(
            parameter="nk1",
            value=len(kpoints_up) + len(kpoints_dn),
            namelist="&kpoints",
        )

    def __str__(self):

        output = [
            "TurboRVB vmcopt python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        pass

    def generate_input(self, input_name: str = "datasmin.input"):
        self.namelist.write(input_name)
        # check if twist_average is manual
        if self.twist_average == 2:  # k-points are set manually
            kpoints_up, kpoints_dn = self.manual_kpoints
            # read
            with open(input_name, "r") as f:
                lines = f.readlines()
            kpoint_index = [
                True if re.match(r".*&kpoints.*", line) else False
                for line in lines
            ].index(True)
            insert_lineno = -1
            for i, line in enumerate(lines[kpoint_index + 1 :]):
                if re.match(r".*/.*", line):
                    insert_lineno = kpoint_index + 1 + i + 1
                    break
                if re.match(r".*&.*", line):
                    insert_lineno = kpoint_index + 1 + i
                    break
            if insert_lineno == -1:
                logger.error("No suitable position for KPOINT is found")
                raise ValueError
            # write dn
            lines.insert(insert_lineno, "\n")
            for kx, ky, kz, wk in reversed(kpoints_dn):
                lines.insert(insert_lineno, f"{kx} {ky} {kz} {wk}\n")
            lines.insert(insert_lineno, "\n")
            # write up
            for kx, ky, kz, wk in reversed(kpoints_up):
                lines.insert(insert_lineno, f"{kx} {ky} {kz} {wk}\n")
            lines.insert(insert_lineno, "KPOINTS\n")
            lines.insert(insert_lineno, "\n")
            # saved
            with open(input_name, "w") as f:
                f.writelines(lines)
        logger.info(f"{os.path.basename(input_name)} has been generated.")

    def run(
        self, input_name: str = "datasmin.input", output_name: str = "out_min"
    ):
        run(
            turbo_qmc_run_command,
            input_name=input_name,
            output_name=output_name,
        )
        remove_file(file="stop.dat")

    def check_results(self, output_names: Optional[list] = None):
        if output_names is None:
            output_names = ["out_min"]
        flags = []
        for output_name in output_names:
            file_check(output_name)
            with open(output_name, "r") as f:
                lines = f.readlines()
            if any(
                [re.match(r".*Final.*tstep.*found.*", line) for line in lines]
            ):
                flags.append(True)
            else:
                flags.append(False)
        return flags

    def plot_energy_and_devmax(
        self, output_names: Optional[list] = None, interactive: bool = False
    ):
        if output_names is None:
            output_names = ["out_min"]
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

    def get_energy(self, output_names: Optional[list] = None):
        if output_names is None:
            output_names = ["out_min"]
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

    def get_devmax(self, output_names: Optional[list] = None):
        if output_names is None:
            output_names = ["out_min"]
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

    def get_estimated_time_for_1_generation(
        self, output_names: Optional[list] = None
    ):
        if output_names is None:
            output_names = ["out_min"]

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

    def plot_parameters_history(self, interactive: bool = True):
        # save parameters
        current_dir = os.getcwd()
        cmd = f"(echo '1 1 0 0'; echo '0'; echo '100000') | {os.path.join(turborvb_bin_root, 'readalles.x')}"
        run(binary=cmd, output_name="out_readalles_for_plot_all")
        shutil.copyfile("./story.d", "./all_story.d")
        all_story_pandas_data = pd.read_csv(
            "./all_story.d",
            header=None,
            delim_whitespace=True,
            index_col=False,
        )

        graph_save_dir = os.path.join(current_dir, "parameters_graphs")
        os.makedirs(graph_save_dir, exist_ok=True)

        if interactive:
            logger.info(
                "How to stop to show the graphs, please type ctrl+C, then close the currently shown graph."
            )

        for i in range(1, len(all_story_pandas_data.columns)):
            plt.rcParams["font.family"] = "sans-serif"
            plt.rcParams["xtick.direction"] = "in"
            plt.rcParams["ytick.direction"] = "in"
            plt.rcParams["xtick.major.width"] = 1.0
            plt.rcParams["ytick.major.width"] = 1.0
            plt.rcParams["font.size"] = 12
            plt.rcParams["axes.linewidth"] = 1.5
            plt.figure(figsize=(10, 8))
            plt.plot(
                all_story_pandas_data[0].astype("int"),
                all_story_pandas_data[i].astype("float"),
                color="black",
                marker="o",
                linestyle="dashed",
                label="all",
            )
            plt.xlabel("Iteration", fontname="Times New Roman", fontsize=14)
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
            if interactive:
                try:
                    plt.waitforbuttonpress()
                except KeyboardInterrupt:
                    logger.info("KeyboardInterrupt")
                    break
            plt.savefig(
                os.path.join(
                    graph_save_dir, "Parameter_No{}_all.png".format(i)
                ),
                bbox_inches="tight",
                pad_inches=0.2,
            )
            plt.close()

    def average_optimized_parameters(
        self,
        equil_steps: int = 10,
        input_file_used: str = "datasmin.input",
        graph_plot: bool = False,
    ):

        """
        if self.twist_average:
            logger.error("Sorry! average_optimized_parameters is not implemented for twist_average = true at present.")
            logger.error("Please do it manually.")
            raise NotImplementedError
        """

        current_dir = os.getcwd()

        # averaging
        with open("forces.dat", "r") as f:
            lines = f.readlines()
        length_iter = len(lines)
        logger.info(f"The total number of iterations is {length_iter}")
        assert length_iter > equil_steps

        # remove previous parameters before average.
        logger.info("Removing previous averaged parameters")
        remove_new_parameter_lines_in_fort10()  # workaround!! to be refactored in the near future.

        # average and show the result
        logger.info("Averaging parameters using readalles.x")
        logger.info(
            f"The first {equil_steps} iterations are disregarded in the average."
        )
        cmd = f"(echo '1 {equil_steps + 1} 1 0'; echo '0'; echo '100000') | {os.path.join(turborvb_bin_root, 'readalles.x')}"
        run(binary=cmd, output_name="out_readalles_for_average")
        shutil.copyfile("./story.d", "./average_story.d")
        ave_story_pandas_data = pd.read_csv(
            "./average_story.d",
            header=None,
            delim_whitespace=True,
            index_col=False,
        )

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
                        graph_save_dir, "Parameter_No{}_averaged.png".format(i)
                    ),
                    bbox_inches="tight",
                    pad_inches=0.2,
                )
                plt.close()

        # update fort10 using the averaged parameters
        logger.info("Update fort10 using the averaged parameters")
        logger.info(f"The input file used for the opt. is {input_file_used}")
        vmcopt = VMCopt.parse_from_file(
            file=input_file_used,
            in_fort10=self.in_fort10,
            twist_average=self.twist_average,
        )
        vmcopt.set_parameter("iopt", 1)
        vmcopt.set_parameter("ngen", 0)
        if self.twist_average:
            vmcopt.set_parameter("yes_kpoints", ".false.", "&parameters")

        # replace iopt from 0 -> 1
        io_fort10 = IO_fort10(self.in_fort10)
        io_fort10.f10header.io_flag = 1

        dir_ave_temp = os.path.join(current_dir, "ave_temp")
        os.makedirs(dir_ave_temp, exist_ok=True)

        shutil.copy(
            os.path.join(os.getcwd(), "fort.10"),
            os.path.join(dir_ave_temp, "fort.10"),
        )
        vmcopt.generate_input(os.path.join(dir_ave_temp, "ave.input"))

        # to be refactored with some flag
        if io_fort10.pp_flag:
            file_check_flag("pseudo.dat")
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
        vmcopt_default_file = os.path.join(
            pyturbo_data_dir, "vmcopt", "datasmin.input"
        )
        namelist = Namelist.parse_namelist_from_file(vmcopt_default_file)
        return namelist

    @staticmethod
    def read_namelist_from_file(file: str):
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist

    @classmethod
    def parse_from_default_namelist(
        cls, in_fort10: str = "fort.10", twist_average: bool = False
    ):
        namelist = cls.read_default_namelist()
        return cls(
            in_fort10=in_fort10, namelist=namelist, twist_average=twist_average
        )

    @classmethod
    def parse_from_file(
        cls, file, in_fort10: str = "fort.10", twist_average: bool = False
    ):
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
