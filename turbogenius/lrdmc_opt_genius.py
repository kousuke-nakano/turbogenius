#!python
# -*- coding: utf-8 -*-
"""

lrdmc_opt genius related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# python modules
import os

# Logger
from logging import getLogger, StreamHandler, Formatter

# turbogenius modules
from turbogenius.pyturbo.lrdmcopt import LRDMCopt
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.utils_workflows.env import turbo_genius_root
from turbogenius.utils_workflows.utility import (
    get_optimizer_flags,
    get_nonlocalmoves_setting,
)
from turbogenius.geniusIO import GeniusIO

logger = getLogger("Turbo-Genius").getChild(__name__)


class LRDMCopt_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo LRDMCopt class

    Attributes:
         fort10 (str): fort.10 WF file
         lrdmcoptsteps (int): total number of optimization steps
         steps (int): number of MCMC steps per optimization step
         bin_block (int): binning length
         warmupblocks (int): the number of disregarded blocks,
         num_walkers (int): The number of walkers, -1 (default) = the number of MPI processes
         maxtime (int): Maxtime (sec.)
         optimizer (str): Choose optimizer, selected from sr:stochastic reconfiguration or lr:linear method.
         learning_rate (float): optimization step size, default values=sr:0.05, lr:0.35
         regularization (float): regularization parameter
         alat (float): Lattice space (Bohr)
         time_branching (float): interval between two branching steps (a.u.)
         etry (float): Trial Energy (Ha)
         nonlocalmoves (str): Treatment of locality approximation, choose from "tmove", "dla", "dlatm"
         opt_onebody (bool): flag to optimize onebody Jastrow
         opt_twobody (bool): flag to optimize twobody Jastrow
         opt_det_mat (bool): flag to optimize matrix elements in the determinant part
         opt_jas_mat (bool): flag to optimize matrix elements in the Jastrow part
         opt_det_basis_exp (bool): flag to optimize exponents of the determinant basis sets
         opt_jas_basis_exp (bool): flag to optimize exponents of the Jastrow basis sets
         opt_det_basis_coeff (bool): flag to optimize coefficients of the determinant basis sets
         opt_jas_basis_coeff (bool): flag to optimize coefficients of the Jastrow basis sets
         twist_average (bool): Twist average flag, True or False
         kpoints (list): k Monkhorst-Pack grids, [kx,ky,kz,nx,ny,nz], kx,y,z-> grids, nx,y,z-> shift=0, noshift=1.
    """

    def __init__(
        self,
        fort10: str = "fort.10",
        lrdmcoptsteps: int = 100,
        steps: int = 10,
        bin_block: int = 1,
        warmupblocks: int = 1,
        num_walkers: int = -1,  # default -1 -> num of MPI process.
        maxtime: int = 172800,
        optimizer: str = "sr",
        learning_rate: float = 0.02,
        regularization: float = 0.001,
        alat: float = -0.20,
        time_branching: float = 0.10,
        etry: float = 0.0,
        nonlocalmoves: str = "tmove",  # tmove, dla, dlatm
        opt_onebody: bool = True,
        opt_twobody: bool = True,
        opt_det_mat: bool = False,
        opt_jas_mat: bool = True,
        opt_det_basis_exp: bool = False,
        opt_jas_basis_exp: bool = False,
        opt_det_basis_coeff: bool = False,
        opt_jas_basis_coeff: bool = False,
        twist_average: bool = False,
        kpoints: list = [1, 1, 1, 0, 0, 0],
    ):

        self.fort10 = fort10
        self.twist_average = twist_average
        self.kpoints = kpoints

        self.energy = None
        self.energy_error = None
        self.estimated_time_for_1_generation = None

        (
            optimizer_number,
            iesdonebodyoff,
            iesdtwobodyoff,
            twobodyoff,
            iesd,
            iesfree,
            iessw,
            iesup,
            iesm,
        ) = get_optimizer_flags(
            optimizer=optimizer,
            opt_onebody=opt_onebody,
            opt_twobody=opt_twobody,
            opt_det_mat=opt_det_mat,
            opt_jas_mat=opt_jas_mat,
            opt_det_basis_exp=opt_det_basis_exp,
            opt_jas_basis_exp=opt_jas_basis_exp,
            opt_det_basis_coeff=opt_det_basis_coeff,
            opt_jas_basis_coeff=opt_jas_basis_coeff,
            qmc_type="lrdmc",
        )

        self.lrdmcopt = LRDMCopt.parse_from_default_namelist(
            in_fort10=fort10, twist_average=twist_average
        )

        self.lrdmcopt.set_parameter(
            parameter="itestr4", value=optimizer_number, namelist="&simulation"
        )
        self.lrdmcopt.set_parameter(
            parameter="ngen",
            value=lrdmcoptsteps * steps,
            namelist="&simulation",
        )
        self.lrdmcopt.set_parameter(
            parameter="maxtime", value=maxtime, namelist="&simulation"
        )
        if num_walkers != -1:
            self.lrdmcopt.set_parameter(
                parameter="nw", value=num_walkers, namelist="&simulation"
            )
        self.lrdmcopt.set_parameter(
            parameter="nweight", value=steps, namelist="&optimization"
        )
        self.lrdmcopt.set_parameter(
            parameter="nbinr", value=bin_block, namelist="&optimization"
        )
        self.lrdmcopt.set_parameter(
            parameter="iboot", value=warmupblocks, namelist="&optimization"
        )
        self.lrdmcopt.set_parameter(
            parameter="tpar", value=learning_rate, namelist="&optimization"
        )
        self.lrdmcopt.set_parameter(
            parameter="parr", value=regularization, namelist="&optimization"
        )

        self.lrdmcopt.set_parameter(
            parameter="iesdonebodyoff",
            value=iesdonebodyoff,
            namelist="&optimization",
        )
        self.lrdmcopt.set_parameter(
            parameter="iesdtwobodyoff",
            value=iesdtwobodyoff,
            namelist="&optimization",
        )
        self.lrdmcopt.set_parameter(
            parameter="twobodyoff", value=twobodyoff, namelist="&optimization"
        )

        self.lrdmcopt.set_parameter(
            parameter="iesd", value=iesd, namelist="&parameters"
        )
        self.lrdmcopt.set_parameter(
            parameter="iesfree", value=iesfree, namelist="&parameters"
        )
        self.lrdmcopt.set_parameter(
            parameter="iessw", value=iessw, namelist="&parameters"
        )
        self.lrdmcopt.set_parameter(
            parameter="iesup", value=iesup, namelist="&parameters"
        )
        self.lrdmcopt.set_parameter(
            parameter="iesm", value=iesm, namelist="&parameters"
        )

        self.lrdmcopt.set_parameter(
            parameter="alat", value=alat, namelist="&dmclrdmc"
        )
        self.lrdmcopt.set_parameter(
            parameter="etry", value=etry, namelist="&dmclrdmc"
        )
        self.lrdmcopt.set_parameter(
            parameter="tbra", value=time_branching, namelist="&dmclrdmc"
        )

        typereg, npow = get_nonlocalmoves_setting(nonlocalmoves=nonlocalmoves)
        self.lrdmcopt.set_parameter(
            parameter="typereg", value=typereg, namelist="&dmclrdmc"
        )
        self.lrdmcopt.set_parameter(
            parameter="npow", value=npow, namelist="&dmclrdmc"
        )

        # regularization
        # to be arguments of the class
        self.lrdmcopt.set_parameter(
            parameter="parcutg", value=0, namelist="&dmclrdmc"
        )
        self.lrdmcopt.set_parameter(
            parameter="true_wagner", value=1, namelist="&dmclrdmc"
        )
        self.lrdmcopt.set_parameter(
            parameter="cutweight", value=-1.0e-4, namelist="&dmclrdmc"
        )

        # pseudo integration
        self.lrdmcopt.set_parameter(
            parameter="npsamax", value=4, namelist="&pseudo"
        )

        # kpoints
        if self.twist_average:  # not 0 (= not False)!!
            if self.twist_average == 1:  # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.lrdmcopt.set_parameter(
                    parameter="yes_kpoints",
                    value=".true.",
                    namelist="&parameters",
                )
                # self.lrdmcopt.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.lrdmcopt.set_parameter(
                    parameter="kp_type", value=1, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="nk1", value=nkx, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="nk2", value=nky, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="nk3", value=nkz, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="k1", value=kx, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="k2", value=ky, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="k3", value=kz, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="skip_equivalence",
                    value=".true.",
                    namelist="&kpoints",
                )
                self.lrdmcopt.set_parameter(
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
                self.lrdmcopt.set_parameter(
                    parameter="yes_kpoints",
                    value=".true.",
                    namelist="&parameters",
                )
                # self.lrdmcopt.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.lrdmcopt.set_parameter(
                    parameter="kp_type", value=2, namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="nk1", value=len(kpoints_up), namelist="&kpoints"
                )
                self.lrdmcopt.set_parameter(
                    parameter="double_kpgrid",
                    value=".true.",
                    namelist="&kpoints",
                )
                self.lrdmcopt.manual_kpoints = self.kpoints

            else:
                logger.error(
                    f"twist_average = {self.twist_average} is not implemented."
                )
                raise NotImplementedError

    def run_all(
        self,
        cont: bool = False,
        input_name: str = "datasfn_opt.input",
        output_name: str = "out_fn_opt",
        average_parameters: bool = True,
    ) -> None:
        """
        Generate input files and run the command.

        Args:
            cont (bool): if True, continuation run (i.e., iopt=0), if False, starting from scratch (i.e., iopt=1).
            input_name (str): input file name
            output_name (str): output file name
            average_parameters (bool): if True, average the optimized parameters

        """
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)
        if average_parameters:
            self.average(input_name=input_name, output_name=output_name)

    def generate_input(
        self, cont: bool = False, input_name: str = "datasfn_opt.input"
    ) -> None:
        """
        Generate input file.

        Args:
            cont (bool): if True, continuation run (i.e., iopt=0), if False, starting from scratch (i.e., iopt=1).
            input_name (str): input file name

        """
        io_fort10 = IO_fort10(fort10=self.fort10)
        io_fort10.io_flag = 0
        if cont:
            self.lrdmcopt.set_parameter("iopt", 0, "&simulation")
        self.lrdmcopt.generate_input(input_name=input_name)

    def run(
        self,
        input_name: str = "datasfn_opt.input",
        output_name: str = "out_fn_opt",
    ) -> None:
        """
        Run the command.

        Args:
            input_name (str): input file name
            output_name (str): output file name
        """
        self.lrdmcopt.run(input_name=input_name, output_name=output_name)
        flags = self.lrdmcopt.check_results(output_names=[output_name])
        assert all(flags)

    def store_result(self, output_names: list = ["out_fn_opt"]) -> None:
        """
        Store results. energy, energy_error, and estimated_time_for_1_generation are stored in this class.

        Args:
            output_names (list): a list of output file names
        """
        self.energy, self.energy_error = self.get_energy(
            output_names=output_names
        )
        self.estimated_time_for_1_generation = (
            self.get_estimated_time_for_1_generation(output_names=output_names)
        )

    def plot_energy_and_devmax(
        self, output_names: list = ["out_fn_opt"], interactive: bool = True
    ):
        """
        plot energy and devmax

        Args:
            output_names (list): a list of output file names
            interactive (bool): flag for an interactive plot
        """
        self.lrdmcopt.plot_energy_and_devmax(
            output_names=output_names, interactive=interactive
        )

    def average(
        self,
        optwarmupsteps: int = 10,
        graph_plot: bool = False,
        input_name: str = "datasfn_opt.input",
        output_names: list = ["out_fn_opt"],
    ) -> None:
        """
        Average parameters of fort.10

        Args:
            optwarmupblocks (int): the number of disregarded optimization steps
            input_name (str): the input file used in the latest calculation
            output_names (list): a list of output file names
            graph_plot (bool): Flag for plotting a graph
        """
        flags = self.lrdmcopt.check_results(output_names=output_names)
        assert all(flags)
        self.lrdmcopt.average_optimized_parameters(
            equil_steps=optwarmupsteps,
            input_file_used=input_name,
            graph_plot=graph_plot,
        )

    def get_energy(self, output_names: list = ["out_fn_opt"]) -> list:
        """
        return energy list

        Args:
            output_names (list): a list of output file names

        Return:
            list: a list of history of energies.
        """
        return self.lrdmcopt.get_energy(output_names=output_names)

    def get_estimated_time_for_1_generation(
        self, output_names: list = ["out_fn_opt"]
    ) -> float:
        """
        This procedure stores estimated_time_for_1_generation.

        Args:
            output_names (list): a list of output file names

        Return:
            float: estimated_time_for_1_generation.
        """
        return self.lrdmcopt.get_estimated_time_for_1_generation(
            output_names=output_names
        )

    def check_results(self, output_names: list = ["out_fn_opt"]) -> bool:
        """
        Check the result.

        Args:
            output_names (list): a list of output file names
        Returns:
            bool: True if all the runs were successful, False if an error is detected in the files.
        """
        return self.lrdmcopt.check_results(output_names=output_names)

    def plot_parameters_history(self, interactive: bool = True) -> None:
        """
        plot history of optimized variational parameters

        Args:
            interactive (bool): flag for an interactive plot
        """
        self.lrdmcopt.plot_parameters_history(interactive=interactive)


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

    os.chdir(os.path.join(turbo_genius_root, "tests", "lrdmcopt"))

    lrdmcopt_genius = LRDMCopt_genius(
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
