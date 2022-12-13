#!python
# -*- coding: utf-8 -*-

"""

Prep related classes and methods

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

#python modules
import os, sys
import numpy as np
import click
import pickle

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pyturbo.prep import Prep
from pyturbo.io_fort10 import IO_fort10

#turbo-genius modules
from utils_workflows.env import turbo_genius_root
from geniusIO import GeniusIO

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)

from turbo_genius_cli import cli, decorate_grpost, header, OptionEatAll
@cli.command(short_help = "prep_genius")
@decorate_grpost
@click.option("-grid", "grid_size",
              help= 'Specify grid size',
              nargs = 3,
              default = [0.1, 0.1, 0.1],
              type = float)
@click.option("-lbox", "lbox",
              help= 'Specify lbox (only for molecules)',
              nargs = 3,
              default = [15.0, 15.0, 15.0],
              type = float)
@click.option("-smear", "smearing",
              help= 'Specify the smearing',
              default = 0.0,
              type = float)
@click.option("-f", "h_field",
              help= 'Specify the h_field',
              default = 0.0,
              type = float)
@click.option("-m", "magnetic_moment_list",
              help= 'Specify the magnetic_moment_list for each atom (u(up), d(dn), or n(no))',
              #nargs=1000,
              default=(),
              type = list,
              cls=OptionEatAll
              )
@click.option("-max", "maxtime",
              help= 'Specify maxtime (sec.)',
              default = 3600,
              type = int)
@click.option("-xc", "xc",
              help= 'Specify xc (lda or lsda)',
              default = 'lda',
              type = str)
@click.option("-kpts", "kpoints",
              help= 'kpts, Specify Monkhorst-Pack grids and shifts, [nkx,nky,nkz,kx,ky,kz]',
              nargs=6,
              default=[0,0,0,0,0,0],
              type=int
              )
@header
def prep(
            g:bool,r:bool,post:bool,
            operation:bool,
            log_level:bool,
            lbox:list,
            grid_size:list,
            smearing:float,
            h_field:float,
            magnetic_moment_list:list,
            xc:str,
            maxtime:int,
            kpoints:list,
)->None:
    """

    prep class launched by turbogenius_cli

        Args:
            See Prep_genius arguments.

    """
    pkl_name="dft_genius_cli.pkl"
    root_dir=os.getcwd()
    pkl_file=os.path.join(root_dir, pkl_name)

    if g:
        if sum(kpoints)==0:
            twist_average = False
        else:
            twist_average = True

        for i, m in enumerate(magnetic_moment_list):
            if m in {"u", "up"}:
                magnetic_moment_list[i]=1
            elif m in {"d", "dn"}:
                magnetic_moment_list[i]=-1
            elif m in {"n", "no"}:
                magnetic_moment_list[i]=0
            else:
                raise ValueError

        dft_genius=DFT_genius(
            maxtime=maxtime,
            grid_size=grid_size,
            lbox=lbox,
            smearing=smearing,
            h_field=h_field,
            magnetic_moment_list=magnetic_moment_list,
            xc=xc,
            twist_average=twist_average,
            kpoints=kpoints
        )
        dft_genius.generate_input()
        with open(pkl_file, "wb") as f:
            pickle.dump(dft_genius, f)

    if r:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                dft_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        dft_genius.run()

    if post:
        os.chdir(root_dir)
        try:
            with open(pkl_file, "rb") as f:
                dft_genius=pickle.load(f)
        except FileNotFoundError:
            logger.error("Did you generate your input file using turbogenius?")
            raise FileNotFoundError
        flags=dft_genius.check_results()
        if all(flags):
            logger.info("Job was successful. DFT is converged.")
        else:
            logger.info("Job was failure. See the output file.")
            return
        dft_genius.check_results()

class DFT_genius(GeniusIO):
    """

    This class is a wrapper of pyturbo prep class

    Attributes:
         fort.10 (str): fort.10 WF file
         grid_size (list):  3 floats, grid sizes [x,y,z]
         lbox (list):  3 floats, Box sizes [x,y,z] (angstrom)
         smearing (float): smearing parameter (Ha)
         maxtime (int): maximum time (sec.)
         memlarge (bool): use more memory to speed up
         maxit (int): maximum iterations
         epsdft (float): Tolerance in the convergence of total energy
         h_field (float): magnetic field putting on each grid.
         magnetic_moment_list (list): magnetic moment list, for all atoms.
         xc (str): Exchange correlation functionals, lda or lsda
         twist_average (bool): Twist average flag, True or False.
         independent_kpoints (bool): Independent kpoint calculation, True or False
         kpoints (list): k Monkhorst-Pack grids, [kx,ky,kz,nx,ny,nz], kx,y,z-> grids, nx,y,z-> shift=0, noshift=1.
    """
    def __init__(self,
                 fort10:str="fort.10",
                 grid_size:list=[0.1, 0.1, 0.1],
                 lbox:list=[15.0, 15.0, 15.0],
                 smearing:float=0.0,
                 maxtime:int=172800,
                 memlarge:bool=False,
                 maxit:int=50,
                 epsdft:float=1.0e-5,
                 h_field:float=0.0,
                 magnetic_moment_list:list=[],
                 xc:str='lda', # lda or lsda
                 twist_average:bool=False,
                 independent_kpoints:bool=False,
                 kpoints:list=[1,1,1,0,0,0]
                 ):

        self.fort10 = fort10
        self.grid_a, self.grid_b, self.grid_c = grid_size
        self.lbox_a, self.lbox_b, self.lbox_c = lbox
        self.smearing = smearing
        self.maxtime = maxtime
        self.memlarge = memlarge
        self.maxit = maxit
        self.epsdft = epsdft
        self.h_field = h_field
        self.magnetic_moment_list=magnetic_moment_list
        self.xc = xc
        self.twist_average = twist_average
        self.independent_kpoints = independent_kpoints
        self.kpoints = kpoints

        io_fort10 = IO_fort10(self.fort10)
        #self.io_fort10 = IO_fort10(self.fort10)
        # this should not be an attribute!! because fort.10 is sometimes very large.

        if io_fort10.f10structure.pbc_flag:
            # for crystals, Lx, Ly, and Lz are cells
            self.Lx=io_fort10.f10structure.norm_vec_a
            self.Ly=io_fort10.f10structure.norm_vec_b
            self.Lz=io_fort10.f10structure.norm_vec_c
            logger.info("Lbox is the norms of the lattice vectors")
            logger.info(f"Lx={self.Lx}, Ly={self.Ly}, Lz={self.Lz}")
            self.ax=self.grid_a
            self.ay=self.grid_b
            self.az=self.grid_c
            self.nx = int(self.Lx/self.ax)
            self.ny = int(self.Ly/self.ay)
            self.nz = int(self.Lz/self.az)
        else:
            # +- 7.5 bohr from the edges.
            pos = io_fort10.f10structure.positions
            self.Lx = np.max(pos[:, 0]) - np.min(pos[:, 0]) + self.lbox_a
            self.Ly = np.max(pos[:, 1]) - np.min(pos[:, 1]) + self.lbox_b
            self.Lz = np.max(pos[:, 2]) - np.min(pos[:, 2]) + self.lbox_c
            logger.info("Actual lboxes are set to +- (the input lboxes)/2.0 (bohr) from the edges of the molecules.")
            logger.info(f"Lx={self.Lx}, Ly={self.Ly}, Lz={self.Lz}")
            self.ax = self.grid_a
            self.ay = self.grid_b
            self.az = self.grid_c
            self.nx = int(self.Lx / self.ax)
            self.ny = int(self.Ly / self.ay)
            self.nz = int(self.Lz / self.az)
            logger.info(f"nx={self.nx}, ny={self.ny}, nz={self.nz}")

        # self prep class!!
        self.prep=Prep.parse_from_default_namelist(in_fort10=fort10, twist_average=self.twist_average)
        
        self.prep.set_parameter(parameter="maxtime", value=maxtime, namelist="&simulation")

        # set L_box
        if io_fort10.f10structure.pbc_flag:
            self.prep.set_parameter(parameter="nx", value=self.nx, namelist="&molecul")
            self.prep.set_parameter(parameter="ny", value=self.ny, namelist="&molecul")
            self.prep.set_parameter(parameter="nz", value=self.nz, namelist="&molecul")
            self.prep.comment_out(parameter="ax")
            self.prep.comment_out(parameter="ay")
            self.prep.comment_out(parameter="az")
        else:
            self.prep.set_parameter(parameter="ax", value=self.ax, namelist="&molecul")
            self.prep.set_parameter(parameter="ay", value=self.ay, namelist="&molecul")
            self.prep.set_parameter(parameter="az", value=self.az, namelist="&molecul")
            self.prep.set_parameter(parameter="nx", value=self.nx, namelist="&molecul")
            self.prep.set_parameter(parameter="ny", value=self.ny, namelist="&molecul")
            self.prep.set_parameter(parameter="nz", value=self.nz, namelist="&molecul")

        ## &dft part

        #contraction
        if io_fort10.det_contraction_flag:
            self.prep.set_parameter(parameter="contracted_on", value=".true.", namelist="&dft")
        else:
            self.prep.set_parameter(parameter="contracted_on", value=".false.", namelist="&dft")

        #xc
        assert self.xc in {'lda', 'lsda'}
        if self.xc == 'lda':
            self.prep.set_parameter(parameter="typedft", value=1, namelist="&dft")
        elif self.xc == 'lsda':
            self.prep.set_parameter(parameter="typedft", value=4, namelist="&dft")
        else:
            logger.error(f"self.xc ={self.xc} is not implemented in TurboRVB.")
            raise NotImplemented

        # memory
        if self.memlarge:
            self.prep.set_parameter(parameter="memlarge", value='.true.', namelist="&dft")
        else:
            self.prep.set_parameter(parameter="memlarge", value='.false.', namelist="&dft")

        # maxit
        self.prep.set_parameter(parameter="maxit", value=self.maxit, namelist="&dft")

        # dft convergence thr
        self.prep.set_parameter(parameter="epsdft", value=self.epsdft, namelist="&dft")

        # smearing!
        if self.smearing == 0.0:
            self.prep.set_parameter(parameter="optocc", value=0, namelist="&dft")
            self.prep.set_parameter(parameter="epsshell", value=self.smearing, namelist="&dft")

            if self.xc == 'lda':
                nelocc=io_fort10.f10header.nelup
                nelocc_list=[2 for _ in range(io_fort10.f10header.neldn)] + [1 for _ in range(io_fort10.f10header.nelup - io_fort10.f10header.neldn)]

                self.prep.set_parameter(parameter="nelocc", value=nelocc, namelist="&dft")
                self.prep.nelocc_list=nelocc_list

            else:
                nelocc = io_fort10.f10header.nelup
                neloccdn = io_fort10.f10header.neldn
                nelocc_list=[1 for _ in range(io_fort10.f10header.nelup)]
                neloccdn_list=[1 for _ in range(io_fort10.f10header.neldn)]
                self.prep.set_parameter(parameter="nelocc", value=nelocc, namelist="&dft")
                self.prep.nelocc_list=nelocc_list

                self.prep.set_parameter(parameter="neloccdo", value=neloccdn, namelist="&dft")
                self.prep.neloccdn_list=neloccdn_list

        else:
            self.prep.set_parameter(parameter="optocc", value=1, namelist="&dft")
            self.prep.set_parameter(parameter="epsshell", value=self.smearing, namelist="&dft")

        if self.h_field != 0.0:
            self.nxs = int(self.nx / 10)
            self.nys = int(self.ny / 10)
            self.nzs = int(self.nz / 10)
            self.nxs = 1 # test
            self.nys = 1 # test
            self.nzs = 2 # test
            self.prep.set_parameter(parameter="h_field", value=self.h_field, namelist="&dft")
            self.prep.set_parameter(parameter="nxs", value=self.nxs, namelist="&dft")
            self.prep.set_parameter(parameter="nys", value=self.nys, namelist="&dft")
            self.prep.set_parameter(parameter="nzs", value=self.nzs, namelist="&dft")

            assert len(self.magnetic_moment_list) == io_fort10.f10header.natom

            self.prep.magnetic_moments_3d_array=self.get_mangetic_moments_3d_array()

        # kpoints
        if self.twist_average: # not 0 (= not False)!!

            if self.independent_kpoints: # True case, independent k calculation, i.e., decoupled_run = True, in Turbo.
                self.prep.set_parameter(parameter="decoupled_run", value=".true.", namelist="&parameters")
            else:
                self.prep.set_parameter(parameter="decoupled_run", value=".false.", namelist="&parameters") #default value

            if self.twist_average == 1: # True case, Monkhorst-Pack algorithm
                assert len(self.kpoints) == 6
                nkx, nky, nkz, kx, ky, kz = self.kpoints
                self.prep.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                self.prep.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.prep.set_parameter(parameter="kp_type", value=1, namelist="&kpoints")
                self.prep.set_parameter(parameter="nk1", value=nkx, namelist="&kpoints")
                self.prep.set_parameter(parameter="nk2", value=nky, namelist="&kpoints")
                self.prep.set_parameter(parameter="nk3", value=nkz, namelist="&kpoints")
                self.prep.set_parameter(parameter="k1", value=kx, namelist="&kpoints")
                self.prep.set_parameter(parameter="k2", value=ky, namelist="&kpoints")
                self.prep.set_parameter(parameter="k3", value=kz, namelist="&kpoints")
                self.prep.set_parameter(parameter="skip_equivalence", value='.true.', namelist="&kpoints")
                """ comment on 9th Dec. double_kpgrid does not play any role with kp_type=1. See kpoints.f90
                    if one wants to use the same phases (up, dn), one should set finite phase to fort.10.
                    TurboRVB automatically detects and sets 'opposite_phase' .true.
                if self.twist_up_dn_opposite_signs:
                    self.prep.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
                else:
                    self.prep.set_parameter(parameter="double_kpgrid", value='.false.', namelist="&kpoints")
                """
            elif self.twist_average == 2: # k-points are set from the user
                assert len(self.kpoints) == 2
                kpoints_up, kpoints_dn = self.kpoints
                assert len(kpoints_up) == len(kpoints_dn)
                for kup, kdn in zip(kpoints_up, kpoints_dn):
                    assert len(kup) == 4 # kx, ky, kz, wkp for up
                    assert len(kdn) == 4 # kx, ky, kz, wkp for dn
                self.prep.set_parameter(parameter="yes_kpoints", value=".true.", namelist="&parameters")
                self.prep.set_parameter(parameter="yeswrite10", value=".true.", namelist="&optimization")
                self.prep.set_parameter(parameter="kp_type", value=2, namelist="&kpoints")
                self.prep.set_parameter(parameter="nk1", value=len(kpoints_up), namelist="&kpoints")
                self.prep.set_parameter(parameter="double_kpgrid", value='.true.', namelist="&kpoints")
                self.prep.manual_kpoints=self.kpoints

            else:
                logger.error(f"twist_average = {self.twist_average} is not implemented.")
                raise NotImplemented

        self.energy=None

    def run_all(self, cont:bool=False, input_name:str="prep.input", output_name:str="out_prep")->None:
        """
            Generate input files and run the command.

            Args:
                input_name (str): input file name
                output_name (str): output file name

        """
        self.generate_input(cont=cont, input_name=input_name)
        self.run(input_name=input_name, output_name=output_name)

    def generate_input(self, cont:bool=False, input_name:str="prep.input")->None:
        """
            Generate input file.

            Args:
                input_name (str): input file name

        """
        if cont: self.prep.set_parameter("iopt", 0, "$systems")
        self.prep.generate_input(input_name=input_name)

    def run(self, input_name:str="prep.input", output_name:str="out_prep")->None:
        """
            Run the command.

            Args:
                input_name (str): input file name
                output_name (str): output file name
        """
        self.prep.run(input_name=input_name, output_name=output_name)
        flags=self.prep.check_results(output_names=[output_name])
        assert all(flags)

    def check_results(self, output_names:list=["out_prep"])->None:
        """
            Check the result.

            Args:
                output_names (list): a list of output file names
            Returns:
                bool: True if all the runs were successful, False if an error is detected in the files.
        """
        return self.prep.check_results(output_names=output_names)

    def get_mangetic_moments_3d_array(self): # -> numpy.array(XXX)
        """
            Return magnetic moment array according to the TurboRVB format.

            Returns:
                numpy.array: numpy array containing the Magnetic moments according to the TurboRVB format.

        """

        io_fort10 = IO_fort10(self.fort10)

        if io_fort10.f10structure.ortho_flag:
            logger.error("Sorry, set magnetic moment is not implemented for a non-orthorhombic case.")
            raise NotImplementedError

        # effective radius of atoms for specifying magnetism. This does not effect the final results apparently.
        radius = 0.50 # bohr # ? 3.0 * (ax+ay+az)/3.0

        nzs, nys, nxs = self.nzs, self.nys, self.nxs
        Lz, Ly, Lx = self.Lz, self.Ly, self.Lx
        mangnetic_moments_3d_array = np.zeros((nzs, nys, nxs))

        #print(io_fort10.f10structure.positions)
        natom = io_fort10.f10header.natom
        x_coord = io_fort10.f10structure.positions[:, 0]
        y_coord = io_fort10.f10structure.positions[:, 1]
        z_coord = io_fort10.f10structure.positions[:, 2]
        logger.debug(x_coord)
        logger.debug(y_coord)
        logger.debug(z_coord)

        if io_fort10.f10structure.pbc_flag:
            # periodic case!! -> the origin is always 0.5,0.5,0.5
            x_shift= 1.0/2.0 * Lx
            y_shift= 1.0/2.0 * Ly
            z_shift= 1.0/2.0 * Lz
            logger.info(f"In periodic cases, the origin is always set to (0.5, 0.5, 0.5).")

        else:
            # isolated molecule case!!
            # compute the center of the molecule # to set the center [0,0,0]
            x_shift = np.mean(x_coord)
            y_shift = np.mean(y_coord)
            z_shift = np.mean(z_coord)
            logger.info(f"The center of the molecule is ({x_shift},{y_shift},{z_shift}) Bohr => shifted to 0")

        logger.info(f"Now, generating magnetic grids....")
        logger.info(self.magnetic_moment_list)

        for z_index in range(nzs):
            #nz_cell = Lz / nzs
            #z = (1 / 2 * nz_cell + z_index * nz_cell) - 1 / 2 * Lz
            z_left = - (1.0 / 2.0 * Lz) + Lz * z_index/nzs
            z_right = - (1.0 / 2.0 * Lz) + Lz * (z_index+1)/nzs
            z = (z_left + z_right)/2.0
            logger.debug(f"z_left={z_left}")
            logger.debug(f"z_right={z_right}")
            logger.debug(f"z={z}")

            for y_index in range(nys):
                #ny_cell = Ly / nys
                #y = (1 / 2 * ny_cell + y_index * ny_cell) - 1 / 2 * Ly
                y_left = - (1.0 / 2.0 * Ly) + Ly * y_index / nys
                y_right = - (1.0 / 2.0 * Ly) + Ly * (y_index + 1) / nys
                y = (y_left + y_right) / 2.0
                logger.debug(f"y_left={y_left}")
                logger.debug(f"y_right={y_right}")
                logger.debug(f"y={y}")
                #print(f"y={y}")

                for x_index in range(nxs):
                    #nx_cell = Lx / nxs
                    #x = (1 / 2 * nx_cell + x_index * nx_cell) - 1 / 2 * Lx
                    x_left = - (1.0 / 2.0 * Lx) + Lx * x_index / nxs
                    x_right = - (1.0 / 2.0 * Lx) + Lx * (x_index + 1) / nxs
                    x = (x_left + x_right) / 2.0
                    logger.debug(f"x_left={x_left}")
                    logger.debug(f"x_right={x_right}")
                    logger.debug(f"x={x}")

                    magnetic_moment_list = []

                    for num in range(natom):
                        x_coord_shifted = x_coord[num] - x_shift
                        y_coord_shifted = y_coord[num] - y_shift
                        z_coord_shifted = z_coord[num] - z_shift

                        logger.debug(x_coord_shifted)
                        logger.debug(y_coord_shifted)
                        logger.debug(z_coord_shifted)

                        if io_fort10.f10structure.pbc_flag:
                            while x_coord_shifted < -1.0 / 2.0 * Lx:
                                x_coord_shifted += Lx
                            while 1.0 / 2.0 * Lx < x_coord_shifted:
                                x_coord_shifted -= Lx

                            while y_coord_shifted < -1.0 / 2.0 * Ly:
                                y_coord_shifted += Ly
                            while 1.0 / 2.0 * Ly < y_coord_shifted:
                                y_coord_shifted -= Ly

                            while z_coord_shifted < -1.0 / 2.0 * Lz:
                                z_coord_shifted += Lz
                            while 1.0 / 2.0 * Lz < z_coord_shifted:
                                z_coord_shifted -= Lz

                        if x_coord_shifted < x_left:
                            x_diff = np.abs(x_left - x_coord_shifted)
                        elif x_right < x_coord_shifted:
                            x_diff = np.abs(x_coord_shifted - x_right)
                        else:
                            x_diff = 0
                        logger.debug(f"x_diff={x_diff}")

                        if y_coord_shifted < y_left:
                            y_diff = np.abs(y_left - y_coord_shifted)
                        elif y_right < y_coord_shifted:
                            y_diff = np.abs(y_coord_shifted - y_right)
                        else:
                            y_diff = 0
                        logger.debug(f"y_diff={y_diff}")

                        if z_coord_shifted < z_left:
                            z_diff = np.abs(z_left - z_coord_shifted)
                        elif z_right < z_coord_shifted:
                            z_diff = np.abs(z_coord_shifted - z_right)
                        else:
                            z_diff = 0
                        logger.debug(f"z_diff={z_diff}")

                        dist = np.sqrt(x_diff ** 2 + y_diff ** 2 + z_diff ** 2)

                        if dist <= radius:
                            # return 1
                            logger.debug("inside radius!! sys.exit()")
                            #sys.exit()
                            magnetic_moment_list.append(self.magnetic_moment_list[num])

                    # check if unique!!
                    # if not, it means that the same grid belongs to more than two atoms with diffeerent moments.
                    # try increase the mesh_factor (i.e., larger nxs, nyx, and nyz)
                    if len(set(magnetic_moment_list)) == 0:
                        mangnetic_moments_3d_array[z_index, y_index, x_index] = 0
                    elif len(set(magnetic_moment_list)) == 1:
                        mangnetic_moments_3d_array[z_index, y_index, x_index] = magnetic_moment_list[0]
                    else:
                        logger.error("the same grid belongs to more than two atoms.")
                        logger.error("try increase nx, ny, and nz)")
                        raise ValueError

        logger.debug(mangnetic_moments_3d_array)
        #sys.exit()

        return mangnetic_moments_3d_array

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "prep"))

    # moved to examples
    prep_genius=DFT_genius(grid_size=0.1)
    prep_genius.generate_input()
    
