#!python
# -*- coding: utf-8 -*-

"""

lrdmc genius related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# python modules
import os
from typing import Optional

# Logger
from logging import getLogger, StreamHandler, Formatter

# turbogenius modules
from turbogenius.pyturbo.lrdmc import LRDMC
from turbogenius.utils_workflows.env import turbo_genius_root
from turbogenius.utils_workflows.utility import get_nonlocalmoves_setting
from turbogenius.geniusIO import GeniusIO

logger = getLogger("Turbo-Genius").getChild(__name__)


class LRDMC_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo LRDMC class

    Attributes:
         fort10 (str): fort.10 WF file
         lrdmcsteps (int): total number of MCMC steps.
         alat (float): Lattice space (Bohr)
         time_branching: interval between two branching steps. (a.u.)
         etry (float): Trial Energy (Ha)
         num_walkers (int): The number of walkers, -1 (default) = the number of MPI processes
         maxtime (int): Maxtime (sec.)
         twist_average (bool): Twist average flag, True or False
         kpoints (list): k Monkhorst-Pack grids, [kx,ky,kz,nx,ny,nz], kx,y,z-> grids, nx,y,z-> shift=0, noshift=1.
         force_calc_flag (bool): if True, compute energy and force, if False, compute only energy
         nonlocalmoves (str): Treatment of locality approximation, choose from "tmove", "dla", "dlatm"
    """

    def __init__(
        self,
        fort10: str = "fort.10",
        lrdmcsteps: int = 100,
        alat: float = -0.20,
        time_branching: float = 0.10,
        etry: float = 0.0,
        num_walkers: int = -1,  # default -1 -> num of MPI process.
        maxtime: int = 172800,
        twist_average: bool = False,
        kpoints: Optional[list] = None,
        force_calc_flag: bool = False,
        nonlocalmoves: str = "tmove",  # tmove, dla, dlatm
    ):
        if kpoints is None:
            kpoints = [1, 1, 1, 0, 0, 0]

        self.force_calc_flag = force_calc_flag
        self.twist_average = twist_average
        self.kpoints = kpoints

        self.estimated_time_for_1_generation = None

        self.energy = None
        self.energy_error = None
        self.forces = None  # np.array([[]]) # 3 * natom matrix
        self.forces_error = None  # np.array([[]])  # 3 * natom matrix

        self.lrdmc = LRDMC.parse_from_default_namelist(
            in_fort10=fort10, twist_average=twist_average
        )
        self.lrdmc.set_parameter(
            parameter="ngen", value=lrdmcsteps, namelist="&simulation"
        )
        self.lrdmc.set_parameter(
            parameter="maxtime", value=maxtime, namelist="&simulation"
        )
        if num_walkers != -1:
            self.lrdmc.set_parameter(
                parameter="nw", value=num_walkers, namelist="&simulation"
            )

        self.lrdmc.set_parameter(
            parameter="etry", value=etry, namelist="&dmclrdmc"
        )
        self.lrdmc.set_parameter(
            parameter="alat", value=alat, namelist="&dmclrdmc"
        )
        self.lrdmc.set_parameter(
            parameter="tbra", value=time_branching, namelist="&dmclrdmc"
        )
        typereg, npow = get_nonlocalmoves_setting(nonlocalmoves=nonlocalmoves)
        self.lrdmc.set_parameter(
            parameter="typereg", value=typereg, namelist="&dmclrdmc"
        )
        self.lrdmc.set_parameter(
            parameter="npow", value=npow, namelist="&dmclrdmc"
        )

        # Do you want to compute forces?
        if not self.force_calc_flag:
            pass
        else:
            self.lrdmc.set_parameter(
                parameter="ieskin", value=1, namelist="&parameters"
            )
            # to be arguments of the class
            self.lrdmc.set_parameter(
                parameter="parcutg", value=0, namelist="&dmclrdmc"
            )
            self.lrdmc.set_parameter(
                parameter="true_wagner", value=1, namelist="&dmclrdmc"
            )
            self.lrdmc.set_parameter(
                parameter="cutweight", value=-1.0e-4, namelist="&dmclrdmc"
            )

        # pseudo integration
        self.lrdmc.set_parameter(
            parameter="npsamax", value=4, namelist="&pseudo"
        )

        # kpoints
        if self.twist_average:  # not 0 (= not False)!!
            if self.twist_average == 1:  # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.lrdmc.set_parameter(
                    parameter="yes_kpoints",
                    value=".true.",
                    namelist="&parameters",
                )
                # self.lrdmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.lrdmc.set_parameter(
                    parameter="kp_type", value=1, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="nk1", value=nkx, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="nk2", value=nky, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="nk3", value=nkz, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="k1", value=kx, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="k2", value=ky, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="k3", value=kz, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="skip_equivalence",
                    value=".true.",
                    namelist="&kpoints",
                )
                self.lrdmc.set_parameter(
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
                self.lrdmc.set_parameter(
                    parameter="yes_kpoints",
                    value=".true.",
                    namelist="&parameters",
                )
                # self.lrdmc.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.lrdmc.set_parameter(
                    parameter="kp_type", value=2, namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="nk1", value=len(kpoints_up), namelist="&kpoints"
                )
                self.lrdmc.set_parameter(
                    parameter="double_kpgrid",
                    value=".true.",
                    namelist="&kpoints",
                )
                self.lrdmc.manual_kpoints = self.kpoints

            else:
                logger.error(
                    f"twist_average = {self.twist_average} is not implemented."
                )
                raise NotImplementedError

    def run_all(
        self,
        bin_block: int = 10,
        warmupblocks: int = 2,
        correcting_factor: int = 2,
        cont: bool = False,
        input_name: str = "datasfn.input",
        output_name: str = "out_fn",
    ) -> None:
        """
        Generate input files and run the command.

        Args:
            bin_block (int): binning length
            warmupblocks (int): the number of disregarded blocks,
            correcting_factor (int): correcting factors
            cont (bool): if True, continuation run (i.e., iopt=0), if False, starting from scratch (i.e., iopt=1).
            input_name (str): input file name
            output_name (str): output file name

        """
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        self.compute_energy_and_forces(
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            correcting_factor=correcting_factor,
        )

    def generate_input(
        self, cont: bool = False, input_name: str = "datasfn.input"
    ) -> None:
        """
        Generate input file.

        Args:
            cont (bool): if True, continuation run (i.e., iopt=0), if False, starting from scratch (i.e., iopt=1).
            input_name (str): input file name

        """
        if cont:
            self.lrdmc.set_parameter("iopt", 0, "&simulation")
        self.lrdmc.generate_input(input_name=input_name)

    def run(
        self, input_name: str = "datasfn.input", output_name: str = "out_fn"
    ) -> None:
        """
        Run the command.

        Args:
            input_name (str): input file name
            output_name (str): output file name
        """
        self.lrdmc.run(input_name=input_name, output_name=output_name)
        flags = self.lrdmc.check_results(output_names=[output_name])
        assert all(flags)

    def store_result(
        self,
        bin_block: int = 10,
        warmupblocks: int = 2,
        correcting_factor: int = 2,
        output_names: Optional[list] = None,
        rerun: bool = False,
    ) -> None:
        """
        Store results. This procedure stores estimated_time_for_1_generation, energy, and energy_error.
        This method is needed for storing data and access to them later.

        Args:
            bin_block (int): binning length
            warmupblocks (int): the number of disregarded blocks
            correcting_factor (int): correcting factors
            output_names (list): a list of output file names
            rerun (bool): if true, compute energy and force again even if there are energy and force files.
        """
        if output_names is None:
            output_names = ["out_fn"]
        self.estimated_time_for_1_generation = (
            self.get_estimated_time_for_1_generation(output_names=output_names)
        )
        self.energy, self.energy_error = self.lrdmc.get_energy(
            init=warmupblocks,
            correct=correcting_factor,
            bin=bin_block,
            rerun=rerun,
        )

    def compute_energy_and_forces(
        self,
        bin_block: int = 10,
        warmupblocks: int = 2,
        correcting_factor: int = 2,
        rerun: bool = False,
    ) -> None:
        """
        Compute energy and forces

        Args:
            bin_block (int): binning length
            warmupblocks (int): the number of disregarded blocks
            correcting_factor (int): correcting factors
            rerun (bool): if true, compute energy and force again even if there are energy and force files.
        """
        self.energy, self.energy_error = self.lrdmc.get_energy(
            init=warmupblocks,
            correct=correcting_factor,
            bin=bin_block,
            rerun=rerun,
        )

    def get_estimated_time_for_1_generation(
        self, output_names: Optional[list] = None
    ) -> float:
        """
        This procedure stores estimated_time_for_1_generation.

        Args:
            output_names (list): a list of output file names

        Return:
            float: estimated_time_for_1_generation.
        """
        if output_names is None:
            output_names = ["out_fn"]
        return self.lrdmc.get_estimated_time_for_1_generation(
            output_names=output_names
        )

    def check_results(self, output_names: Optional[list] = None) -> bool:
        """
        Check the result.

        Args:
            output_names (list): a list of output file names
        Return:
            bool: True if all the runs were successful, False if an error is detected in the files.
        """
        if output_names is None:
            output_names = ["out_fn"]
        return self.lrdmc.check_results(output_names=output_names)


if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("INFO")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "lrdmc"))

    # moved to examples
