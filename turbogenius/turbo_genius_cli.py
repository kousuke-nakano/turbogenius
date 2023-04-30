#!usr/bin/env python
# coding: utf-8

"""

TurboGenius command line interface

"""

import os
import shutil
import pickle
import click
from typing import Union
from logging import getLogger, StreamHandler, Formatter

from turbogenius.makefort10_genius import Makefort10_genius
from turbogenius.convertfort10mol_genius import Convertfort10mol_genius
from turbogenius.convertfort10_genius import Convertfort10_genius
from turbogenius.convertpfaff_genius import Convertpfaff_genius
from turbogenius.prep_genius import DFT_genius
from turbogenius.vmc_genius import VMC_genius
from turbogenius.lrdmc_genius import LRDMC_genius
from turbogenius.vmc_opt_genius import VMCopt_genius
from turbogenius.lrdmc_opt_genius import LRDMCopt_genius
from turbogenius.readforward_genius import Readforward_genius
from turbogenius.correlated_sampling_genius import Correlated_sampling_genius
from turbogenius.wavefunction import Wavefunction
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.database_setup import (
    all_electron_basis_set_list,
    ccECP_basis_set_list,
    BFD_basis_set_list,
)
from turbogenius.database_setup import ecp_list

try:
    from turbogenius._version import version as turbogenius_version
except (ModuleNotFoundError, ImportError):
    turbogenius_version = "unknown"

logger = getLogger("Turbo-Genius").getChild(__name__)


# Thanks Oto :-)) by Kosuke Nakano
class OptionEatAll(click.Option):
    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop("save_other_options", True)
        nargs = kwargs.pop("nargs", -1)
        assert nargs == -1, "nargs, if set, must be -1 not {}".format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval


def header(f):
    @click.option(
        "-log",
        "log_level",
        default="INFO",
        help="logger level, DEBUG, INFO, ERROR",
    )
    def ret(*args, **kwargs):

        # set loggers
        logger_t = getLogger("Turbo-Genius")
        logger_t.setLevel(kwargs["log_level"])
        stream_handler = StreamHandler()
        # handler_format = Formatter('%(message)s')
        handler_format = Formatter(
            "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
        )
        stream_handler.setFormatter(handler_format)
        logger_t.addHandler(stream_handler)

        logger_p = getLogger("pyturbo")
        logger_p.setLevel(kwargs["log_level"])
        stream_handler_p = StreamHandler()
        # handler_format_p = Formatter('%(message)s')
        handler_format_p = Formatter(
            "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
        )
        stream_handler_p.setFormatter(handler_format_p)
        logger_p.addHandler(stream_handler_p)

        logger_t.info(f"Turbo-Genius {turbogenius_version}")

        operation = ",".join(
            [
                val
                for val, com in zip(
                    ["generate", "run", "postprocess"], ["g", "r", "post"]
                )
                if com in kwargs and kwargs[com]
            ]
        )
        kwargs["operation"] = operation
        if any([x in kwargs for x in ["r", "g", "post"]]):
            if operation == "":
                logger_t.info(
                    "Jobtype is needed please type --help for more information"
                )
                logger_t.info("")

            else:
                logger_t.info(f"Selected jobtypes = {kwargs['operation']}")
        f(*args, **kwargs)

    ret.__name__ = f.__name__
    ret.__doc__ = f.__doc__

    return ret


def decorate_grpost(f):
    f = click.option("-g", "g", is_flag=True, help="Generate an input file")(f)
    f = click.option("-r", "r", is_flag=True, help="Run a program")(f)
    f = click.option("-post", "post", is_flag=True, help="Postprocess")(f)
    return f


@click.group()
def cli():
    """Turbo-Genius command-line tool"""


# --------------------------------------------------------------
# makefort10 command
# --------------------------------------------------------------
@cli.command(short_help="makefort10_genius")
@decorate_grpost
@click.option(
    "-str",
    "structure_file",
    help="Specify a structure file.",
    default=None,
    type=str,
)
@click.option(
    "-s",
    "supercell",
    help="Specify a supercell.",
    nargs=3,
    default=[1, 1, 1],
    type=int,
)
@click.option(
    "-detbasis",
    "det_basis_sets",
    help=f"Specify a basis set for the determinant part. "
    f"For all-electrons:{all_electron_basis_set_list}. "
    f"For PPs: ccECP:{ccECP_basis_set_list}, BFD:{BFD_basis_set_list}",
    default="cc-pVTZ",
    type=click.Choice(
        all_electron_basis_set_list + ccECP_basis_set_list + BFD_basis_set_list
    ),
)
@click.option(
    "-jasbasis",
    "jas_basis_sets",
    help=f"Specify a basis set for the Jastrow part {all_electron_basis_set_list}",
    default="cc-pVDZ",
    type=click.Choice(
        all_electron_basis_set_list + ccECP_basis_set_list + BFD_basis_set_list
    ),
)
@click.option(
    "-detcont",
    "det_contracted_flag",
    help="Contraction flag for the determinant part",
    is_flag=True,
    type=bool,
)
@click.option(
    "-jascont",
    "jas_contracted_flag",
    help="Contraction flag for the jastrow part",
    is_flag=True,
    type=bool,
)
@click.option(
    "-jasallele",
    "all_electron_jas_basis_set",
    help="all-electron flag for the jastrow basis set",
    is_flag=True,
    type=bool,
)
@click.option(
    "-pp",
    "pseudo_potential",
    help=f"Pseudopotential, {ecp_list} implemented.",
    default=None,
    type=click.Choice(ecp_list),
)
@click.option(
    "-detcutbasis",
    "det_cut_basis_option",
    help="Cutting the determinant basis set according to a default criteria",
    is_flag=True,
    type=bool,
)
@click.option(
    "-jascutbasis",
    "jas_cut_basis_option",
    help="Cutting the jastrow basis set according to a default criteria",
    is_flag=True,
    type=bool,
)
@click.option(
    "-complex",
    "complex",
    help="Flag for complex WF",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-phaseup",
    "phaseup",
    help="phase, e.g. 0.0 0.0 0.0",
    nargs=3,
    default=[0.0, 0.0, 0.0],
    type=float,
)
@click.option(
    "-phasedn",
    "phasedn",
    help="phase, e.g. 0.0 0.0 0.0",
    nargs=3,
    default=[0.0, 0.0, 0.0],
    type=float,
)
@click.option(
    "-neldiff",
    "neldiff",
    help="Diff. between up and dn electrons",
    default=0,
    type=int,
)
@header
def makefort10(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    structure_file: str,
    supercell: list,
    det_basis_sets: str,
    jas_basis_sets: str,
    det_contracted_flag: bool,
    jas_contracted_flag: bool,
    all_electron_jas_basis_set: bool,
    pseudo_potential: Union[str, None],
    det_cut_basis_option: bool,
    jas_cut_basis_option: bool,
    complex: bool,
    phaseup: list,
    phasedn: list,
    neldiff: int,
) -> None:
    """makefort10

    makefort10 class launched by turbogenius_cli

        Args:
            See Makefort10_genius arguments.

    """
    pkl_name = "makefort10_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if pseudo_potential is not None:
        logger.warning(
            "If you want to use all-electron basis for the Jastrow part, plz. use --jasallele option"
        )
    if g:
        os.chdir(root_dir)
        makefort10_genius = Makefort10_genius(
            structure_file=structure_file,
            supercell=supercell,
            det_basis_set=det_basis_sets,
            jas_basis_set=jas_basis_sets,
            det_contracted_flag=det_contracted_flag,
            jas_contracted_flag=jas_contracted_flag,
            all_electron_jas_basis_set=all_electron_jas_basis_set,
            pseudo_potential=pseudo_potential,
            det_cut_basis_option=det_cut_basis_option,
            jas_cut_basis_option=jas_cut_basis_option,
            complex=complex,
            phase_up=phaseup,
            phase_dn=phasedn,
            neldiff=neldiff,
        )
        makefort10_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(makefort10_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                makefort10_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        makefort10_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                makefort10_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags = makefort10_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return

        logger.info("Rename the generated fort.10_new as fort.10")
        shutil.move("fort.10_new", "fort.10")


# --------------------------------------------------------------
# convertfort10mol command
# --------------------------------------------------------------
@cli.command(short_help="convertfort10mol_genius")
@decorate_grpost
@click.option(
    "--random_mo",
    "add_random_mo",
    help="flag for adding random MOs",
    is_flag=True,
    default=True,
    type=bool,
)
@click.option("--add_mo", "additional_mo", help="additional MOs", default=0, type=int)
@click.option(
    "--grid_size",
    "grid_size",
    help="specify grid_size",
    default=0.10,
    type=float,
)
@header
def convertfort10mol(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    add_random_mo: bool,
    additional_mo: int,
    grid_size: float,
):
    pkl_name = "convertfort10mol_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        convertfort10mol_genius = Convertfort10mol_genius(
            add_random_mo=add_random_mo,
            additional_mo=additional_mo,
            grid_size=grid_size,
        )
        convertfort10mol_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(convertfort10mol_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                convertfort10mol_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        convertfort10mol_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                convertfort10mol_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags = convertfort10mol_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return

        logger.info("Rename the generated fort.10_new as fort.10")
        shutil.move("fort.10_new", "fort.10")


# --------------------------------------------------------------
# convertfort10 command
# --------------------------------------------------------------
@cli.command(short_help="convertfort10_genius")
@decorate_grpost
@click.option("-grid", "grid_size", help="Specify grid size", default=0.1, type=float)
@header
def convertfort10(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    grid_size: list,
):
    pkl_name = "convertfort10_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        convertfort10_genius = Convertfort10_genius(
            grid_size=grid_size,
        )
        convertfort10_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(convertfort10_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                convertfort10_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        convertfort10_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                convertfort10_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags = convertfort10_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        convertfort10_genius.check_results()


# --------------------------------------------------------------
# prep command
# --------------------------------------------------------------
@cli.command(short_help="prep_genius")
@decorate_grpost
@click.option(
    "-grid",
    "grid_size",
    help="Specify grid size",
    nargs=3,
    default=[0.1, 0.1, 0.1],
    type=float,
)
@click.option(
    "-lbox",
    "lbox",
    help="Specify lbox (only for molecules)",
    nargs=3,
    default=[15.0, 15.0, 15.0],
    type=float,
)
@click.option(
    "-smear", "smearing", help="Specify the smearing", default=0.0, type=float
)
@click.option("-f", "h_field", help="Specify the h_field", default=0.0, type=float)
@click.option(
    "-m",
    "magnetic_moment_list",
    help="Specify the magnetic_moment_list for each atom (u(up), d(dn), or n(no))",
    # nargs=1000,
    default=(),
    type=list,
    cls=OptionEatAll,
)
@click.option("-max", "maxtime", help="Specify maxtime (sec.)", default=3600, type=int)
@click.option("-xc", "xc", help="Specify xc (lda or lsda)", default="lda", type=str)
@click.option(
    "-kpts",
    "kpoints",
    help="kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]",
    nargs=6,
    default=[0, 0, 0, 0, 0, 0],
    type=int,
)
@header
def prep(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: bool,
    lbox: list,
    grid_size: list,
    smearing: float,
    h_field: float,
    magnetic_moment_list: list,
    xc: str,
    maxtime: int,
    kpoints: list,
) -> None:
    """

    prep class launched by turbogenius_cli

        Args:
            See Prep_genius arguments.

    """
    pkl_name = "dft_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        if sum(kpoints) == 0:
            twist_average = False
        else:
            twist_average = True

        for i, m in enumerate(magnetic_moment_list):
            if m in {"u", "up"}:
                magnetic_moment_list[i] = 1
            elif m in {"d", "dn"}:
                magnetic_moment_list[i] = -1
            elif m in {"n", "no"}:
                magnetic_moment_list[i] = 0
            else:
                raise ValueError

        dft_genius = DFT_genius(
            maxtime=maxtime,
            grid_size=grid_size,
            lbox=lbox,
            smearing=smearing,
            h_field=h_field,
            magnetic_moment_list=magnetic_moment_list,
            xc=xc,
            twist_average=twist_average,
            kpoints=kpoints,
        )
        dft_genius.generate_input()
        with open(pkl_file, "wb") as f:
            pickle.dump(dft_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                dft_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        dft_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                dft_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags = dft_genius.check_results()
        if all(flags):
            logger.info("Job was successful. DFT is converged.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        dft_genius.check_results()


# --------------------------------------------------------------
# VMCopt command
# --------------------------------------------------------------
@cli.command(short_help="vmcopt_genius")
@decorate_grpost
@click.option(
    "-vmcoptsteps",
    "vmcoptsteps",
    help="Specify vmcoptsteps",
    default=1000,
    type=int,
)
@click.option(
    "-optwarmup",
    "optwarmupsteps",
    help="Specify optwarmupsteps",
    default=1,
    type=int,
)
@click.option(
    "-steps",
    "steps",
    help="Specify steps per one iteration",
    default=20,
    type=int,
)
@click.option("-bin", "bin_block", help="Specify bin_block", default=1, type=int)
@click.option(
    "-warmup", "warmupblocks", help="Specify warmupblocks", default=0, type=int
)
@click.option("-nw", "num_walkers", help="Specify num_walkers", default=-1, type=int)
@click.option("-maxtime", "maxtime", help="Specify maxtime", default=3600, type=int)
@click.option(
    "-optimizer",
    "optimizer",
    help="Specify optimizer, sr or lr",
    default="lr",
    type=str,
)
@click.option(
    "-learn",
    "learning_rate",
    help="Specify learning_rate",
    default=0.35,
    type=float,
)
@click.option(
    "-reg",
    "regularization",
    help="Specify regularization",
    default=0.001,
    type=float,
)
@click.option(
    "-opt_onebody",
    "opt_onebody",
    help="flag for opt_onebody",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_twobody",
    "opt_twobody",
    help="flag for opt_twobody",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_det_mat",
    "opt_det_mat",
    help="flag for opt_det_mat",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_jas_mat",
    "opt_jas_mat",
    help="flag for opt_jas_mat",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_det_basis_exp",
    "opt_det_basis_exp",
    help="flag for opt_det_basis_exp",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_jas_basis_exp",
    "opt_jas_basis_exp",
    help="flag for opt_jas_basis_exp",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_det_basis_coeff",
    "opt_det_basis_coeff",
    help="flag for opt_det_basis_coeff",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_jas_basis_coeff",
    "opt_jas_basis_coeff",
    help="flag for opt_jas_basis_coeff",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-opt_structure",
    "opt_structure",
    help="flag for opt_structure",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-strlearn",
    "str_learning_rate",
    help="Specify str_learning_rate",
    default=1.0e-6,
    type=float,
)
@click.option(
    "-twist",
    "twist_average",
    help="flag for twist_average",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-kpts",
    "kpoints",
    help="kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]",
    nargs=6,
    default=[0, 0, 0, 0, 0, 0],
    type=int,
)
@click.option(
    "-plot",
    "plot_graph",
    help="flag for plotting graph",
    is_flag=True,
    default=False,
    type=bool,
)
@click.option(
    "-interactive",
    "plot_interactive",
    help="flag for interactive plotting graph",
    is_flag=True,
    default=False,
    type=bool,
)
@header
def vmcopt(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    vmcoptsteps: int = 100,
    optwarmupsteps: int = 10,
    steps: int = 10,
    bin_block: int = 1,
    warmupblocks: int = 0,
    num_walkers: int = -1,  # default -1 -> num of MPI process.
    maxtime: int = 172800,
    optimizer: str = "sr",
    learning_rate: float = 0.35,
    regularization: float = 0.001,
    opt_onebody: bool = True,
    opt_twobody: bool = True,
    opt_det_mat: bool = False,
    opt_jas_mat: bool = True,
    opt_det_basis_exp: bool = False,
    opt_jas_basis_exp: bool = False,
    opt_det_basis_coeff: bool = False,
    opt_jas_basis_coeff: bool = False,
    opt_structure: bool = False,
    str_learning_rate: float = 1.0e-6,
    twist_average: bool = False,
    kpoints: list = [1, 1, 1, 0, 0, 0],
    plot_graph: bool = False,
    plot_interactive: bool = False,
) -> None:
    pkl_name = "vmcopt_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        vmcopt_genius = VMCopt_genius(
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
            opt_structure=opt_structure,
            str_learning_rate=str_learning_rate,
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
                vmcopt_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        vmcopt_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                vmcopt_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags = vmcopt_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        vmcopt_genius.plot_energy_and_devmax(interactive=plot_interactive)
        if plot_graph:
            vmcopt_genius.plot_parameters_history(interactive=plot_interactive)
        vmcopt_genius.average(optwarmupsteps=optwarmupsteps, graph_plot=plot_graph)


# --------------------------------------------------------------
# VMC command
# --------------------------------------------------------------
@cli.command(short_help="vmc_genius")
@decorate_grpost
@click.option("-steps", "vmcsteps", help="Specify vmcsteps", default=1000, type=int)
@click.option("-bin", "bin_block", help="Specify bin_block", default=1, type=int)
@click.option(
    "-warmup", "warmupblocks", help="Specify warmupblocks", default=1, type=int
)
@click.option("-nw", "num_walkers", help="Specify num_walkers", default=-1, type=int)
@click.option("-maxtime", "maxtime", help="Specify maxtime", default=3600, type=int)
@click.option(
    "-twist",
    "twist_average",
    help="flag for twist_average",
    is_flag=True,
    type=bool,
)
@click.option(
    "-kpts",
    "kpoints",
    help="kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]",
    nargs=6,
    default=[0, 0, 0, 0, 0, 0],
    type=int,
)
@click.option(
    "-force",
    "force_calc_flag",
    help="flag for force_calc_flag",
    is_flag=True,
    type=bool,
)
@header
def vmc(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    vmcsteps: int,
    bin_block: int,
    warmupblocks: int,
    num_walkers: int,
    maxtime: int,
    twist_average: bool,
    kpoints: list,
    force_calc_flag: bool,
) -> None:
    pkl_name = "vmc_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        vmc_genius = VMC_genius(
            vmcsteps=vmcsteps,
            num_walkers=num_walkers,
            maxtime=maxtime,
            twist_average=twist_average,
            kpoints=kpoints,
            force_calc_flag=force_calc_flag,
        )
        vmc_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(vmc_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                vmc_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        vmc_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                vmc_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags = vmc_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        vmc_genius.check_results()
        vmc_genius.compute_energy_and_forces(
            bin_block=bin_block, warmupblocks=warmupblocks
        )


# --------------------------------------------------------------
# readforward command
# --------------------------------------------------------------
@cli.command(short_help="readforward_genius")
@decorate_grpost
@click.option(
    "-corr",
    "corr_sampling",
    help="correlated sampling",
    is_flag=True,
    default=True,
)
@click.option("-bin", "bin_block", help="Specify bin_block", default=1, type=int)
@click.option(
    "-warmup", "warmupblocks", help="Specify warmupblocks", default=1, type=int
)
@header
def readforward(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    bin_block: int,
    warmupblocks: int,
    corr_sampling: bool,
):
    pkl_name = "readforward_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        readforward_genius = Readforward_genius(
            in_fort10="fort.10",
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            corr_sampling=corr_sampling,
        )
        readforward_genius.generate_input()
        with open(pkl_file, "wb") as f:
            pickle.dump(readforward_genius, f)
    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                readforward_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        readforward_genius.run()

    if post:
        if post:
            os.chdir(root_dir)
            try:
                with open(pkl_file, "rb") as f:
                    readforward_genius = pickle.load(f)
            except FileNotFoundError:
                logger.error("Did you generate your input file using turbogenius?")
                raise FileNotFoundError
            flags = readforward_genius.check_results()
            if all(flags):
                logger.info("Job was successful.")
            else:
                logger.info("Job was failure. See the output file.")
                return
            readforward_genius.check_results()


# --------------------------------------------------------------
# correlated sampling command
# --------------------------------------------------------------
@cli.command(short_help="correlated_sampling_genius")
@decorate_grpost
@click.option("-steps", "vmcsteps", help="Specify vmcsteps", default=1000, type=int)
@click.option("-bin", "bin_block", help="Specify bin_block", default=2, type=int)
@click.option(
    "-warmup", "warmupblocks", help="Specify warmupblocks", default=1, type=int
)
@click.option("-nw", "num_walkers", help="Specify num_walkers", default=-1, type=int)
@click.option("-maxtime", "maxtime", help="Specify maxtime", default=3600, type=int)
@click.option(
    "-twist",
    "twist_average",
    help="flag for twist_average",
    is_flag=True,
    type=bool,
)
@header
def correlated_sampling(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    vmcsteps: int,
    bin_block: int,
    warmupblocks: int,
    num_walkers: int,
    maxtime: int,
    twist_average: bool,
) -> bool:
    pkl_name = "correlated_sampling_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        readforward_genius = Correlated_sampling_genius(
            in_fort10="fort.10",
            corr_fort10="fort.10_corr",
            vmcsteps=vmcsteps,
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            num_walkers=num_walkers,
            maxtime=maxtime,
            twist_average=twist_average,
        )
        readforward_genius.generate_input()
        with open(pkl_file, "wb") as f:
            pickle.dump(readforward_genius, f)
    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                readforward_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        readforward_genius.run()

    if post:
        if post:
            os.chdir(root_dir)
            try:
                with open(pkl_file, "rb") as f:
                    readforward_genius = pickle.load(f)
            except FileNotFoundError:
                logger.error("Did you generate your input file using turbogenius?")
                raise FileNotFoundError
            flags = readforward_genius.check_results()
            if all(flags):
                logger.info("Job was successful.")
            else:
                logger.info("Job was failure. See the output file.")
                return
            readforward_genius.check_results()


# --------------------------------------------------------------
# lrdmcopt command
# --------------------------------------------------------------
@cli.command(short_help="lrdmc_genius")
@decorate_grpost
@click.option(
    "-lrdmcoptsteps",
    "lrdmcoptsteps",
    help="Specify lrdmcoptsteps",
    default=1000,
    type=int,
)
@click.option(
    "-optwarmup",
    "optwarmupsteps",
    help="Specify optwarmupsteps",
    default=1,
    type=int,
)
@click.option(
    "-steps",
    "steps",
    help="Specify steps per one iteration",
    default=1,
    type=int,
)
@click.option("-bin", "bin_block", help="Specify bin_block", default=1, type=int)
@click.option(
    "-warmup", "warmupblocks", help="Specify warmupblocks", default=1, type=int
)
@click.option("-nw", "num_walkers", help="Specify num_walkers", default=-1, type=int)
@click.option("-maxtime", "maxtime", help="Specify maxtime", default=1, type=int)
@click.option(
    "-optimizer",
    "optimizer",
    help="Specify optimizer, sr or lr",
    default="lr",
    type=str,
)
@click.option(
    "-learn",
    "learning_rate",
    help="Specify learning_rate",
    default=0.35,
    type=float,
)
@click.option(
    "-reg",
    "regularization",
    help="Specify regularization",
    default=0.001,
    type=float,
)
@click.option("-alat", "alat", help="Specify alat", default=-0.20, type=float)
@click.option("-etry", "etry", help="Specify etry", default=0.0, type=float)
@click.option(
    "-nonlocal",
    "nonlocalmoves",
    help="Specify nonlocalmoves, tmove, dla, dlatm",
    default="tmove",
    type=str,
)
@click.option(
    "-opt_onebody",
    "opt_onebody",
    help="flag for opt_onebody",
    is_flag=True,
    type=bool,
)
@click.option(
    "-opt_twobody",
    "opt_twobody",
    help="flag for opt_twobody",
    is_flag=True,
    type=bool,
)
@click.option(
    "-opt_det_mat",
    "opt_det_mat",
    help="flag for opt_det_mat",
    is_flag=False,
    type=bool,
)
@click.option(
    "-opt_jas_mat",
    "opt_jas_mat",
    help="flag for opt_jas_mat",
    is_flag=True,
    type=bool,
)
@click.option(
    "-opt_det_basis_exp",
    "opt_det_basis_exp",
    help="flag for opt_det_basis_exp",
    is_flag=False,
    type=bool,
)
@click.option(
    "-opt_jas_basis_exp",
    "opt_jas_basis_exp",
    help="flag for opt_jas_basis_exp",
    is_flag=False,
    type=bool,
)
@click.option(
    "-opt_det_basis_coeff",
    "opt_det_basis_coeff",
    help="flag for opt_det_basis_coeff",
    is_flag=False,
    type=bool,
)
@click.option(
    "-opt_jas_basis_coeff",
    "opt_jas_basis_coeff",
    help="flag for opt_jas_basis_coeff",
    is_flag=False,
    type=bool,
)
@click.option(
    "-twist",
    "twist_average",
    help="flag for twist_average",
    is_flag=False,
    type=bool,
)
@click.option(
    "-kpts",
    "kpoints",
    help="kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]",
    nargs=6,
    default=[0, 0, 0, 0, 0, 0],
    type=int,
)
@header
def lrdmcopt(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    lrdmcoptsteps: int,
    steps: int,
    bin_block: int,
    warmupblocks: int,
    num_walkers: int,  # default -1 -> num of MPI process.
    maxtime: int,
    optimizer: str,
    learning_rate: float,
    regularization: float,
    alat: float,
    etry: float,
    nonlocalmoves: str,  # tmove, dla, dlatm
    opt_onebody: bool,
    opt_twobody: bool,
    opt_det_mat: bool,
    opt_jas_mat: bool,
    opt_det_basis_exp: bool,
    opt_jas_basis_exp: bool,
    opt_det_basis_coeff: bool,
    opt_jas_basis_coeff: bool,
    twist_average: bool,
    kpoints: list,
):
    if g:
        lrdmcopt_genius = LRDMCopt_genius(
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


# --------------------------------------------------------------
# lrdmc command
# --------------------------------------------------------------
@cli.command(short_help="lrdmc_genius")
@decorate_grpost
@click.option("-steps", "lrdmcsteps", help="Specify lrdmcsteps", default=1000, type=int)
@click.option("-bin", "bin_block", help="Specify bin_block", default=1, type=int)
@click.option(
    "-corr",
    "correcting_factor",
    help="Specify correcting_factor",
    default=1,
    type=int,
)
@click.option("-alat", "alat", help="Specify alat", default=-0.20, type=float)
@click.option("-etry", "etry", help="Specify etry", default=0.0, type=float)
@click.option(
    "-warmup", "warmupblocks", help="Specify warmupblocks", default=1, type=int
)
@click.option("-nw", "num_walkers", help="Specify num_walkers", default=-1, type=int)
@click.option("-maxtime", "maxtime", help="Specify maxtime", default=3600, type=int)
@click.option(
    "-twist",
    "twist_average",
    help="flag for twist_average",
    is_flag=True,
    type=bool,
)
@click.option(
    "-force",
    "force_calc_flag",
    help="flag for force_calc_flag",
    is_flag=True,
    type=bool,
)
@click.option(
    "-nonlocal",
    "nonlocalmoves",
    help="Specify nonlocalmoves, tmove, dla, dlatm",
    default="tmove",
    type=str,
)
@header
def lrdmc(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    lrdmcsteps: int,
    bin_block: int,
    warmupblocks: int,
    correcting_factor: int,
    alat: float,
    etry: float,
    num_walkers: int,  # default -1 -> num of MPI process.
    maxtime: int,
    twist_average: bool,
    force_calc_flag: bool,
    nonlocalmoves: str,  # tmove, dla, dlatm
):
    pkl_name = "lrdmc_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        lrdmc_genius = LRDMC_genius(
            lrdmcsteps=lrdmcsteps,
            alat=alat,
            etry=etry,
            num_walkers=num_walkers,  # default -1 -> num of MPI process.
            maxtime=maxtime,
            twist_average=twist_average,
            force_calc_flag=force_calc_flag,
            nonlocalmoves=nonlocalmoves,  # tmove, dla, dlatm
        )
        lrdmc_genius.generate_input()

        with open(pkl_file, "wb") as f:
            pickle.dump(lrdmc_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                lrdmc_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        lrdmc_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                lrdmc_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags = lrdmc_genius.check_results()
        if all(flags):
            logger.info("Job was successful.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        lrdmc_genius.check_results()
        lrdmc_genius.compute_energy_and_forces(
            bin_block=bin_block,
            warmupblocks=warmupblocks,
            correcting_factor=correcting_factor,
        )


# --------------------------------------------------------------
# convertpfaff
# --------------------------------------------------------------
@cli.command(short_help="readforward_genius")
@decorate_grpost
@click.option("-rotate", "rotate_flag", help="rotate_flag", is_flag=True, default=False)
@click.option(
    "-angle",
    "rotate_angle",
    help="Specify rotate_angle",
    default=0.0,
    type=float,
)
@click.option(
    "-scale",
    "scale_mean_field",
    help="Specify scale_mean_field",
    default=1000,
    type=int,
)
@header
def convertpfaff(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    rotate_flag: bool,
    rotate_angle: float,
    scale_mean_field: int,
):
    pkl_name = "convertpfaff_genius_cli.pkl"
    root_dir = os.getcwd()
    pkl_file = os.path.join(root_dir, pkl_name)

    if g:
        convertpfaff_genius = Convertpfaff_genius(
            in_fort10="fort.10_in", out_fort10="fort.10_out"
        )
        convertpfaff_genius.generate_input()
        with open(pkl_file, "wb") as f:
            pickle.dump(convertpfaff_genius, f)
    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                convertpfaff_genius = pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        convertpfaff_genius.run(
            rotate_flag=rotate_flag,
            rotate_angle=rotate_angle,
            scale_mean_field=scale_mean_field,
        )

    if post:
        if post:
            os.chdir(root_dir)
            try:
                with open(pkl_file, "rb") as f:
                    convertpfaff_genius = pickle.load(f)
            except FileNotFoundError:
                logger.error("Did you generate your input file using turbogenius?")
                raise FileNotFoundError
            flags = convertpfaff_genius.check_results()
            if all(flags):
                logger.info("Job was successful.")
            else:
                logger.info("Job was failure. See the output file.")
                return
            convertpfaff_genius.check_results()


# --------------------------------------------------------------
# convert WFs
# --------------------------------------------------------------
@cli.command(short_help="convert wavefunction")
@decorate_grpost
@click.option(
    "-to",
    "to_ansatz",
    help="Specify to ansatz",
    default="agps",
    type=click.Choice(["sd", "agps", "agpu", "pf"], case_sensitive=False),
)
@click.option(
    "-hyb",
    "hybrid_orbitals",
    help="Specify the number of added hybrid orbitals for each atom (0, 5, ...)",
    default=(),
    type=list,
    cls=OptionEatAll,
)
@click.option(
    "-nosym",
    "nosymmetry",
    help="flag for nosymmetry)",
    is_flag=True,
    default=False,
)
@click.option(
    "-rot",
    "rotate_angle",
    help="Specify rotate angle (-0.5 - + 0.5)",
    default=0.00,
    type=float,
)
@click.option("-grid", "grid_size", help="Specify grid size", default=0.10, type=float)
@header
def convertwf(
    g: bool,
    r: bool,
    post: bool,
    operation: bool,
    log_level: str,
    to_ansatz: str,
    hybrid_orbitals: list,
    grid_size: float,
    nosymmetry: bool,
    rotate_angle: float,
):
    # pkl_name = "wavefunction_cli.pkl"
    root_dir = os.getcwd()
    os.chdir(root_dir)
    number_of_additional_hybrid_orbitals = list(map(int, hybrid_orbitals))

    wavefunction = Wavefunction()
    wavefunction.from_fort10(fort10="fort.10")

    if to_ansatz == "agps":
        logger.info("convert the wf to agps")
        wavefunction.to_agps(
            grid_size=grid_size,
            additional_hyb=number_of_additional_hybrid_orbitals,
            nosym=nosymmetry,
            clean_flag=True,
        )

    elif to_ansatz == "agpu":
        logger.info("convert the wf to agpu")
        wavefunction.to_agpu(
            grid_size=grid_size,
            additional_hyb=number_of_additional_hybrid_orbitals,
            nosym=nosymmetry,
            clean_flag=True,
        )

    elif to_ansatz == "sd":
        logger.info("convert the wf to jsd (i.e. adding MOs)")
        wavefunction.to_sd(grid_size=grid_size, clean_flag=False)

    elif to_ansatz == "pf":
        logger.info("convert the wf to pfaffian")
        logger.info(f"rotate_angle = {rotate_angle}")
        wavefunction.to_pf(
            grid_size=grid_size,
            rotate_angle=rotate_angle,
            nosym=nosymmetry,
            clean_flag=True,
        )


# --------------------------------------------------------------
# command-line tools
# --------------------------------------------------------------
@cli.command(short_help="visualize fort.10 structure by ASE")
@header
def view(
    operation: bool,
    log_level: str,
) -> None:
    """
    Visualize a molecule or crystal structure written in fort.10

    """
    io_fort10 = IO_fort10(fort10="fort.10")
    structure = io_fort10.f10structure.structure
    structure.view()


@cli.command(short_help="convert fort.10 to a structure file by ASE-read")
@click.option(
    "-s",
    "structure",
    help="Specify structure file name (e.g., xxxx.xsf)",
    default="ase.xsf",
    type=str,
)
@header
def writestr(operation: bool, log_level: str, structure: str):
    """
    Write a structure file (i.e., convert fort.10 to a structure file, e.g., xyz)

    Args:
        structure (str): structure file name. all formats supported by ASE are acceptable.

    """
    io_fort10 = IO_fort10(fort10="fort.10")
    structure_ = io_fort10.f10structure.structure
    structure_.write(structure)


if __name__ == "__main__":
    cli()
