#!python
# -*- coding: utf-8 -*-

"""

pyturbo: vmc related classes and methods

Todo:
    * docstrings are not completed.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.
    * implementing __str__ method.
    * implementing sanity_check method.

"""

# python modules
import os
import sys
import re
import numpy as np

# turbo-genius modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from namelist import Namelist
from fortranIO import FortranIO
from io_fort10 import IO_fort10
from utils.env import pyturbo_root, pyturbo_data_dir
from utils.env import (
    turbo_qmc_run_command,
    turbo_forcevmc_run_command,
    turbo_forcevmc_kpoints_run_command,
    turbo_forcevmc_kpoints_para_run_command,
)
from utils.utility import (
    file_check,
    file_check_flag,
    get_line_from_file,
    remove_file,
)
from utils.execute import run

# Logger
from logging import config, getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)
# logger = getLogger(__name__)


class VMC(FortranIO):
    def __init__(
        self,
        in_fort10="fort.10",
        namelist=Namelist(),
        twist_average=False,  # False or 0: single-k, True or 1: Monkhorst-Pack, 2: manual k-grid
    ):

        """
        input values
        """
        file_check(in_fort10)
        self.in_fort10 = in_fort10
        self.namelist = namelist
        self.twist_average = twist_average

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
            parameter="nk1", value=len(kpoints_up), namelist="&kpoints"
        )
        self.namelist.set_parameter(
            parameter="double_kpgrid", value=".true.", namelist="&kpoints"
        )

    def __str__(self):

        output = [
            "TurboRVB vmc python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        pass

    def generate_input(self, input_name="datasvmc.input"):
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
            lines.insert(insert_lineno, f"\n")
            for kx, ky, kz, wk in reversed(kpoints_dn):
                lines.insert(insert_lineno, f"{kx} {ky} {kz} {wk}\n")
            lines.insert(insert_lineno, f"\n")
            # write up
            for kx, ky, kz, wk in reversed(kpoints_up):
                lines.insert(insert_lineno, f"{kx} {ky} {kz} {wk}\n")
            lines.insert(insert_lineno, f"KPOINTS\n")
            lines.insert(insert_lineno, f"\n")
            # saved
            with open(input_name, "w") as f:
                f.writelines(lines)
        logger.info(f"{input_name} has been generated.")

    def run(self, input_name="datasvmc.input", output_name="out_vmc"):
        remove_file(file="pip0.d")
        remove_file(file="forces.dat")
        run(
            turbo_qmc_run_command,
            input_name=input_name,
            output_name=output_name,
        )

    def check_results(self, output_names=["out_vmc"]):
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

    def get_estimated_time_for_1_generation(self, output_names=["out_vmc"]):

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

    @staticmethod
    def read_energy(twist_average=False):
        if twist_average:
            line = get_line_from_file(file="pip0.d", line_no=0).split()
            energy = float(line[3])
            error = float(line[5])
        else:
            line = get_line_from_file(file="pip0.d", line_no=1).split()
            energy = float(line[2])
            error = float(line[3])
        return energy, error

    def get_energy(self, init=10, bin=10, num_proc=-1, rerun=False):
        force_compute_flag = False
        if rerun:
            force_compute_flag = True
        else:
            if file_check_flag("pip0.d"):
                with open("pip0.d") as f:
                    pip0 = f.read()
                if os.stat("pip0.d").st_size == 0:
                    force_compute_flag = True
                elif "ERROR" in pip0:
                    force_compute_flag = True
                else:
                    force_compute_flag = False
            else:  # if pip0.d does not exist. always computed.
                force_compute_flag = True

        if force_compute_flag:
            self.compute_energy_and_forces(
                init=init, bin=bin, num_proc=num_proc
            )

        energy, error = self.read_energy(twist_average=self.twist_average)

        logger.debug("energy={}, error={}".format(energy, error))
        return energy, error

    def get_forces(self, init=10, bin=10, num_proc=-1, rerun=False):

        force_compute_flag = False
        if rerun:
            force_compute_flag = True
        else:
            if file_check_flag("forces_vmc.dat"):
                with open("forces_vmc.dat") as f:
                    fff = f.read()
                if os.stat("forces_vmc.dat").st_size == 0:
                    force_compute_flag = True
                elif "ERROR" in fff:
                    force_compute_flag = True
                else:
                    force_compute_flag = False
            else:  # if forces_vmc.dat' does not exist. always computed.
                force_compute_flag = True

        if force_compute_flag:
            self.compute_energy_and_forces(
                init=init, bin=bin, num_proc=num_proc
            )

        file_check(self.in_fort10)
        fort10 = IO_fort10(self.in_fort10)
        force_file = "forces_vmc.dat"
        with open(force_file, "r") as f:
            force_dat = f.readlines()

        # force constant matrix
        force_matrix = np.zeros((fort10.f10header.natom, 3))
        force_matrix_error_bar = np.zeros((fort10.f10header.natom, 3))

        logger.info(
            f"The number of forces calculated = {fort10.f10forceconstraint.ieskinr}"
        )

        for num_force in range(fort10.f10forceconstraint.ieskinr):
            constraint_num = fort10.f10forceconstraint.constraint_num[
                num_force
            ]

            # reading force from force_vmc.dat
            if self.twist_average:  # with k point
                force_element = float(force_dat[num_force].split()[5])
                force_element_error_bar = float(
                    force_dat[num_force].split()[7]
                )
            else:  # gamma point
                # since the number of intervals depends on compiler...
                # start_index = force_dat.index([s for s in force_dat if re.match('.*Force component.*', s)][0])
                # end_index = force_dat.index([s for s in force_dat if re.match('.*<OH>.*-.*<O><H>.*', s)][0])
                # interval = end_index - start_index + 1
                # force_element = float(force_dat[1 + interval * num_force].split()[2])
                # force_element_error_bar = float(force_dat[1 + interval * num_force].split()[3])

                force_element_list = [
                    s for s in force_dat if re.match(".*Force.*=.*", s)
                ]
                force_element = float(force_element_list[num_force].split()[2])
                force_element_error_bar = float(
                    force_element_list[num_force].split()[3]
                )

            # multiplied according to the symmetry
            # note! They are "correlated", we shoulud divide the error bar by force_constraints_d["num_constraints"]
            # not by sqrt(force_constraints_d["num_constraints"])... Right?
            force_element = force_element / constraint_num
            force_element_error_bar = force_element_error_bar / constraint_num

            # add the result to force_matrix
            # check if the selected sym element list includes minus signs.
            nindex_list = [
                i
                for i, x in enumerate(
                    fort10.f10forceconstraint.constraint_index
                )
                if x == num_force
            ]

            for nindex in nindex_list:
                atom_label = fort10.f10forceconstraint.atom_label[nindex]
                direction = fort10.f10forceconstraint.direction[nindex]

                # caution!! force_element can be opposite? No, no. it is correct.
                if atom_label < 0:
                    force_element = -1 * force_element
                else:
                    pass

                force_matrix[
                    np.abs(atom_label) - 1, direction - 1
                ] = force_element
                force_matrix_error_bar[
                    np.abs(atom_label) - 1, direction - 1
                ] = force_element_error_bar

        # unit (Ha/au)
        return force_matrix, force_matrix_error_bar

    def compute_energy_and_forces(self, init=10, bin=10, pulay=1, num_proc=-1):
        if self.twist_average:
            if num_proc == -1:
                logger.warning(
                    f"num_proc is -1. The maximum possible cpus are used for computing energies and forces."
                )
                logger.warning(
                    f"num_proc is set to {os.cpu_count()}, which is obtained by os.cpu_count()"
                )
                num_proc = os.cpu_count()
            if num_proc > 1:
                command = turbo_forcevmc_kpoints_para_run_command
                #command = turbo_forcevmc_kpoints_run_command  # for the time being!!! because paperoga does not have sufficient memory.
            else:
                command = turbo_forcevmc_kpoints_run_command
            cmd = "{:s} {:d} {:d} {:d} {:d}".format(
                command, bin, init * -1, pulay, num_proc
            )
        else:  # twist_flag is False:
            cmd = "{:s} {:d} {:d} {:d}".format(
                turbo_forcevmc_run_command, bin, init * -1, pulay
            )
        logger.info(f"cmd={cmd}")
        run(binary=cmd, output_name="forcevmc.out")

    @staticmethod
    def read_default_namelist():
        vmc_default_file = os.path.join(
            pyturbo_data_dir, "vmc", "datasvmc.input"
        )
        namelist = Namelist.parse_namelist_from_file(vmc_default_file)
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

    # moved to examples.
