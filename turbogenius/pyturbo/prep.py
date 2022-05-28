#!python
# -*- coding: utf-8 -*-

"""prep module

Todo:
    Todo list

"""

#python modules
import os, sys
import re
import numpy as np

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from namelist import Namelist
from fortranIO import FortranIO
from utils.env import pyturbo_data_dir
from utils.env import turbo_prep_run_command
from utils.utility import file_check
from utils.execute import run
from io_fort10 import IO_fort10

from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo').getChild(__name__)

class Prep(FortranIO):
    """Prep

    Prep class

    Attributes:
        in_fort10 (str): filrname of fort.10, always fort.10
        namelist (namelist): namelist of prep.input

    Examples, how to use:
        See tests dir.

    """

    def __init__(self,
                 in_fort10='fort.10',
                 namelist=Namelist(),
                 nelocc_list=[],
                 neloccdn_list=[],
                 magnetic_moments_3d_array=[],  # dim = 3, shape = (nzs,nys,nxs)
                 twist_average=False,  # False or 0: single-k, True or 1: Monkhorst-Pack, 2: manual k-grid
                 ):

        """
        input values
        """
        file_check(in_fort10)
        if IO_fort10(fort10=in_fort10).pp_flag: file_check("pseudo.dat")

        self.in_fort10=in_fort10
        self.namelist=namelist
        self.nelocc_list=nelocc_list
        self.neloccdn_list=neloccdn_list
        self.magnetic_moments_3d_array=magnetic_moments_3d_array
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
        self.namelist.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
        self.namelist.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
        self.namelist.set_parameter(parameter="kp_type", value=2, namelist="&kpoints")
        self.namelist.set_parameter(parameter="nk1", value=len(kpoints_up) + len(kpoints_dn), namelist="&kpoints")
    @property
    def magnetic_moments_3d_array(self):
        return self.__magnetic_moments_3d_array
    @magnetic_moments_3d_array.setter
    def magnetic_moments_3d_array(self, matrix):
        if len(matrix) != 0:
            try:
                nxs = self.namelist.get_parameter(parameter="nxs")
                nys = self.namelist.get_parameter(parameter="nys")
                nzs = self.namelist.get_parameter(parameter="nzs")
            except KeyError:
                logger.error("nxs, nys, or nzs are not correctly set in namelist!")
                raise KeyError
            assert matrix.ndim == 3
            assert matrix.shape == (nzs, nys, nxs)
        self.__magnetic_moments_3d_array=matrix

    def __str__(self):

        output = [
            "TurboRVB prep python wrapper",
        ]
        return "\n".join(output)

    def sanity_check(self):
        pass

    def generate_input(self, input_name):
        """
        Args:
            input_name (str): Name of input file

        Returns:
            NA

        Yields:
            "input_name" is generated.
        """
        self.namelist.write(input_name)
        # check if twist_average is manual
        if self.twist_average == 2:  # k-points are set manually
            kpoints_up, kpoints_dn = self.manual_kpoints
            # read
            with open(input_name, 'r') as f:
                lines = f.readlines()
            kpoint_index = [True if re.match(r'.*&kpoints.*', line) else False for line in lines].index(True)
            insert_lineno = -1
            for i, line in enumerate(lines[kpoint_index + 1:]):
                if re.match(r'.*/.*', line): insert_lineno = kpoint_index + 1 + i + 1; break
                if re.match(r'.*&.*', line): insert_lineno = kpoint_index + 1 + i; break
            if insert_lineno == -1: logger.error("No suitable position for KPOINT is found"); raise ValueError
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
            with open(input_name, 'w') as f:
                f.writelines(lines)
        # write occupation list
        logger.debug(self.nelocc_list)
        logger.debug(self.neloccdn_list)
        if len(self.nelocc_list) != 0:
            with open(input_name, 'a') as f:
                line=" ".join(list(map(str,self.nelocc_list))) + "\n"
                f.write(line)
        if len(self.neloccdn_list) != 0:
            with open(input_name, 'a') as f:
                line="\n" +  " ".join(list(map(str, self.neloccdn_list))) + "\n"
                f.write(line)
        if len(self.magnetic_moments_3d_array) != 0:
            # write magnetic fields
            logger.info(f"Writing magnetic moments at the end of prep.input file ... \n")
            with open(input_name, 'a') as f:
                f.write("\n")
                nzs = self.namelist.get_parameter(parameter="nzs")
                for z_index in range(nzs):
                    np.savetxt(f, self.magnetic_moments_3d_array[z_index, :, :], "%d")
                    f.write("\n")
            logger.info(f"Magnetic moments written. \n")
        logger.info(f"{input_name} has been generated. \n")

    def run(self, input_name="prep.input", output_name="out_prep"):
        run(turbo_prep_run_command, input_name=input_name, output_name=output_name)

    def check_results(self, output_names=["out_prep"]):
        flags = []
        for output_name in output_names:
            file_check(output_name)
            with open(output_name, "r") as f:
                lines = f.readlines()
            if any([re.match(r".*OK.*Turbo-DFT.*converged.*with.*energy.*tollerance.*", line) for line in lines]):
                flags.append(True)
            else:
                flags.append(False)
        return flags

    @staticmethod
    def read_default_namelist(in_fort10="fort.10"):
        prep_default_file = os.path.join(pyturbo_data_dir, "prep", "prep.input")
        namelist = Namelist.parse_namelist_from_file(prep_default_file)

        #fort.10
        io_fort10 = IO_fort10(fort10=in_fort10)

        #contraction
        if io_fort10.det_contraction_flag:
            namelist.set_parameter(parameter="contracted_on", value=".true.", namelist="&dft")
        else:
            namelist.set_parameter(parameter="contracted_on", value=".false.", namelist="&dft")

        #occ
        typedft=namelist.get_parameter(parameter="typedft")
        optocc=namelist.get_parameter(parameter="optocc")
        io_fort10 = IO_fort10(fort10=in_fort10)
        if optocc == 0:
            if typedft==1: # lda
                nelocc = io_fort10.f10header.nelup
                namelist.set_parameter(parameter="nelocc", value=nelocc, namelist="&dft")
            elif typedft==4: # 4
                nelocc = self.io_fort10.f10header.nelup
                neloccdn = self.io_fort10.f10header.neldn
                namelist.set_parameter(parameter="nelocc", value=nelocc, namelist="&dft")
                namelist.set_parameter(parameter="neloccdo", value=neloccdn, namelist="&dft")

        return namelist
    @staticmethod
    def read_namelist_from_file(file):
        namelist = Namelist.parse_namelist_from_file(file)
        return namelist
    @classmethod
    def parse_from_default_namelist(cls, in_fort10="fort.10", magnetic_moments_3d_array = [], twist_average=False):
        namelist = cls.read_default_namelist(in_fort10=in_fort10)
        #fort10
        io_fort10 = IO_fort10(fort10=in_fort10)
        #occ
        typedft=namelist.get_parameter(parameter="typedft")
        optocc=namelist.get_parameter(parameter="optocc")
        if optocc == 0:
            if typedft==1: # lda
                nelocc_list = [2 for _ in range(io_fort10.f10header.neldn)] + [1 for _ in range(io_fort10.f10header.nelup - io_fort10.f10header.neldn)]
                neloccdn_list = []
            elif typedft==4: # 4
                nelocc_list = [1 for _ in range(io_fort10.f10header.nelup)]
                neloccdn_list = [1 for _ in range(io_fort10.f10header.neldn)]

        else:
            nelocc_list=[]
            neloccdn_list=[]

        return cls(
            in_fort10=in_fort10,
            namelist=namelist,
            nelocc_list = nelocc_list,
            neloccdn_list = neloccdn_list,
            magnetic_moments_3d_array = magnetic_moments_3d_array,
            twist_average=twist_average
        )
    @classmethod
    def parse_from_file(cls, file, in_fort10="fort.10"):
        namelist = Namelist.parse_namelist_from_file(file)
        logger.warning(f"nelocc_list and neloccdn_list are not read from {file}")
        logger.warning(f"magnetic_moments_3d_array is not read from {file}")
        logger.warning(f"KPOINTS is not read from {file}")
        return cls(
            in_fort10=in_fort10,
            namelist=namelist
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