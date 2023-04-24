#!python
# -*- coding: utf-8 -*-

"""

Correlated sampling related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# python modules
import os
from typing import Optional

# Logger
from logging import getLogger, StreamHandler, Formatter

# turbogenius modules
from turbogenius.pyturbo.vmc import VMC
from turbogenius.pyturbo.readforward import Readforward
from turbogenius.utils_workflows.env import turbo_genius_root
from turbogenius.geniusIO import GeniusIO

logger = getLogger("Turbo-Genius").getChild(__name__)


class Correlated_sampling_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo LRDMC class

    Attributes:
         in_fort10 (str): fort.10 WF file
         corr_fort10 (str): fort.10 WF file
         vmcsteps (int): total number of MCMC steps.
         bin_block (int): binning length
         warmupblocks (int): the number of disregarded blocks
         num_walkers (int): The number of walkers, -1 (default) = the number of MPI processes
         maxtime (int): Maxtime (sec.)
         twist_average (bool): Twist average flag, True or False
         kpoints (list): k Monkhorst-Pack grids, [kx,ky,kz,nx,ny,nz], kx,y,z-> grids, nx,y,z-> shift=0, noshift=1.
    """

    def __init__(
        self,
        in_fort10: str = "fort.10_in",
        corr_fort10: str = "fort.10_corr",
        vmcsteps: int = 100,
        bin_block: int = 10,
        warmupblocks: int = 2,
        num_walkers: int = -1,  # default -1 -> num of MPI process.
        maxtime: int = 172800,
        twist_average: bool = False,
        kpoints: Optional[list] = None,
    ):
        if kpoints is None:
            kpoints = [1, 1, 1, 0, 0, 0]

        self.in_fort10 = in_fort10
        self.corr_fort10 = corr_fort10

        self.vmcsteps = vmcsteps
        self.bin_block = bin_block
        self.warmupblocks = warmupblocks
        self.num_walkers = num_walkers
        self.maxtime = maxtime
        self.twist_average = twist_average
        self.kpoints = kpoints

        # VMC
        self.vmc = VMC.parse_from_default_namelist(
            in_fort10=in_fort10, twist_average=twist_average
        )
        if vmcsteps < 40 * bin_block + bin_block * warmupblocks:
            logger.warning(
                f"vmcsteps = {vmcsteps} is too small! < 40 * bin_block + bin_block * warmupblocks = {40 * bin_block + bin_block * warmupblocks}"
            )
            logger.warning(
                f"vmcsteps = {vmcsteps} is set to 40 * bin_block + bin_block * warmupblocks = {40 * bin_block + bin_block * warmupblocks}"
            )
            vmcsteps = 40 * bin_block + bin_block * warmupblocks
        self.vmc.set_parameter(
            parameter="ngen", value=vmcsteps, namelist="&simulation"
        )
        self.vmc.set_parameter(
            parameter="maxtime", value=maxtime, namelist="&simulation"
        )
        if num_walkers != -1:
            self.vmc.set_parameter(
                parameter="nw", value=num_walkers, namelist="&simulation"
            )

        self.vmc.set_parameter(parameter="iread", value=3, namelist="&readio")

        # kpoints
        if self.twist_average:  # not 0 (= not False)!!
            if self.twist_average == 1:  # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.vmc.set_parameter(
                    parameter="yes_kpoints",
                    value=".true.",
                    namelist="&parameters",
                )
                self.vmc.set_parameter(
                    parameter="yeswrite10",
                    value=".true.",
                    namelist="&optimization",
                )
                self.vmc.set_parameter(
                    parameter="kp_type", value=1, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="nkx", value=nkx, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="nky", value=nky, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="nkz", value=nkz, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="kx", value=kx, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="ky", value=ky, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="kz", value=kz, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="skip_equivalence",
                    value=".true.",
                    namelist="&kpoints",
                )
                self.vmc.set_parameter(
                    parameter="double_kpgrid",
                    value=".true.",
                    namelist="&kpoints",
                )
            elif self.twist_average == 2:  # k-points are set from the user
                assert len(self.kpoints) == 2
                kpoints_up, kpoints_dn = self.kpoints
                assert len(kpoints_up) == len(kpoints_dn)
                for kup, kdn in zip(kpoints_up, kpoints_dn):
                    assert len(kup) == 4  # kx, ky, kz, wkp for up
                    assert len(kdn) == 4  # kx, ky, kz, wkp for dn
                self.vmc.set_parameter(
                    parameter="yes_kpoints",
                    value=".true.",
                    namelist="&parameters",
                )
                self.vmc.set_parameter(
                    parameter="yeswrite10",
                    value=".true.",
                    namelist="&optimization",
                )
                self.vmc.set_parameter(
                    parameter="kp_type", value=2, namelist="&kpoints"
                )
                self.vmc.set_parameter(
                    parameter="nk1",
                    value=len(kpoints_up) + len(kpoints_dn),
                    namelist="&kpoints",
                )
                self.vmc.manual_kpoints = self.kpoints

            else:
                logger.error(
                    f"twist_average = {self.twist_average} is not implemented."
                )
                raise NotImplementedError

        # readforward
        self.readforward = Readforward.parse_from_default_namelist(
            in_fort10=in_fort10
        )
        self.readforward.set_parameter(
            parameter="bin_length", value=bin_block, namelist="&corrfun"
        )
        self.readforward.set_parameter(
            parameter="initial_bin", value=warmupblocks, namelist="&corrfun"
        )
        self.readforward.set_parameter(
            parameter="correlated_samp", value=".true.", namelist="&corrfun"
        )

    def run_all(
        self,
        input_name: str = "datasvmc.input",
        vmc_output_name: str = "out_vmc",
        readforward_output_name: str = "out_readforward",
    ) -> None:
        """
        Generate input files and run the command.

        Args:
            input_name (str): input file name
            vmc_output_name (str): vmc output file name
            readforward_output_name (str): readforward output file name

        """
        self.vmc.generate_input(input_name=input_name)
        self.vmc.run(input_name=input_name, output_name=vmc_output_name)
        self.readforward.generate_input(input_name=input_name)
        self.readforward.run(
            input_name=input_name, output_name=readforward_output_name
        )

    def generate_input(self, input_name: str = "datasvmc.input") -> None:
        """
        Generate input file.

        Args:
            input_name (str): input file name

        """
        self.vmc.generate_input(input_name=input_name)
        self.readforward.generate_input(input_name="readforward.input")

    def run(
        self,
        input_name: str = "datasvmc.input",
        vmc_output_name: str = "out_vmc",
        readforward_output_name: str = "out_readforward",
    ) -> None:
        """
        Run the command.

        Args:
            input_name (str): input file name
            vmc_output_name (str): vmc output file name
            readforward_output_name (str): readforward output file name
        """
        self.vmc.run(input_name=input_name, output_name=vmc_output_name)
        flags = self.vmc.check_results(output_names=[vmc_output_name])
        assert all(flags)
        self.readforward.run(
            input_name=input_name, output_name=readforward_output_name
        )
        flags = self.readforward.check_results(
            output_names=[readforward_output_name]
        )
        assert all(flags)

    def check_results(
        self,
        vmc_output_names: Optional[list] = None,
        readforward_output_names: Optional[list] = None,
    ) -> bool:
        """
        Check the result.

        Args:
            vmc_output_names (list): a list of output file names
            readforward_output_names (list): a list of output file names
        Returns:
            bool: True if all the runs were successful, False if an error is detected in the files.
        """
        if vmc_output_names is None:
            vmc_output_names = ["out_vmc"]
        if readforward_output_names is None:
            readforward_output_names = ["out_readforward"]
        return self.readforward.check_results(
            output_names=readforward_output_names
        ) + self.vmc.check_results(output_names=vmc_output_names)


if __name__ == "__main__":
    from turbogenius.convertfort10_genius import Convertfort10_genius

    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "convertfort10"))

    # moved to examples
    convertfort10_genius = Convertfort10_genius()
    convertfort10_genius.generate_input()
    convertfort10_genius.run()
