#!python
# -*- coding: utf-8 -*-

"""

Wavefunction related classes and methods

Todo:
    * Refactoring this module using TurboWorkflows because several tasks in the conversion could be heavy calculations.
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# python modules
from __future__ import annotations
import os
import shutil
import numpy as np
from typing import Union

# Logger
from logging import getLogger, StreamHandler, Formatter

# turbogenius modules
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.pyturbo.pseudopotentials import Pseudopotentials
from turbogenius.pyturbo.makefort10 import Makefort10
from turbogenius.pyturbo.utils.utility import remove_file
from turbogenius.utils_workflows.env import turbo_genius_root
from turbogenius.trexio_to_turborvb import trexio_to_turborvb_wf
from turbogenius.makefort10_genius import Makefort10_genius
from turbogenius.convertfort10mol_genius import Convertfort10mol_genius
from turbogenius.convertfort10_genius import Convertfort10_genius
from turbogenius.convertpfaff_genius import Convertpfaff_genius
from turbogenius.tools_genius import copy_jastrow
from turbogenius.pyturbo.basis_set import Jas_Basis_sets


logger = getLogger("Turbo-Genius").getChild(__name__)


class Wavefunction:
    """

    This class is for manipulating fort.10 using pyturbo IO_fort.10 class.


    """

    def __init__(self):
        self.read_flag = False  # flag for read
        self.wf_format = None  # turborvb or trexio
        self.trexio_filename = None

    def read_from_fort10(self, fort10: str) -> None:
        """

        Read WF from fort.10 file

        Attributes:
            fort10 (str): name of a WF file (usually fort.10)

        """
        if fort10 != "fort.10":
            try:
                shutil.copy(fort10, "fort.10")
            except shutil.SameFileError:
                pass

        self.io_fort10 = IO_fort10(fort10="fort.10")
        self.read_flag = True
        self.wf_format = "turborvb"
        self.trexio_filename = None

    def read_from_trexio(
        self,
        trexio_filename: str,
        jas_basis_sets: Jas_Basis_sets = Jas_Basis_sets(),
        max_occ_conv: int = 0,
        mo_num_conv: int = -1,
        only_mol: bool = True,
        cleanup: bool = True,
    ) -> None:
        """
        Convert trexio file to TurboRVB WF file (fort.10)

        Args:
            trexio_file (str): TREXIO file name
            jas_basis_sets (Jas_basis_sets): Jastrow basis sets added to the TREXIO WF.
            max_occ_conv (int): maximum occ used for the conv, not used with mo_num
            mo_num_conv (int): num mo used for the conv, not used with max occ
            only_mol (bool): if True, only moleculer orbitals option = True in convertfort10mol
            cleanup (bool): clean up temporary files
        """

        if os.path.isfile("fort.10"):
            logger.warning("fort.10 in the current directory is overwritten!")

        trexio_to_turborvb_wf(
            trexio_file=trexio_filename,
            jas_basis_sets=jas_basis_sets,
            max_occ_conv=max_occ_conv,
            mo_num_conv=mo_num_conv,
            only_mol=only_mol,
            cleanup=cleanup,
        )

        self.io_fort10 = IO_fort10(fort10="fort.10")
        self.read_flag = True
        self.wf_format = "trexio"
        self.trexio_filename = trexio_filename

    def read_from_structure(
        self,
        structure_file: str,
        supercell: list = [1, 1, 1],
        det_basis_set: Union[str, list] = "cc-pVQZ",
        jas_basis_set: Union[str, list] = "cc-pVQZ",
        det_contracted_flag: bool = True,
        jas_contracted_flag: bool = True,
        all_electron_jas_basis_set: bool = True,
        pseudo_potential: Union[str, None] = None,
        det_cut_basis_option: bool = False,
        jas_cut_basis_option: bool = False,
        jastrow_type: int = -6,
        complex: bool = False,
        phase_up: list = [0.0, 0.0, 0.0],
        phase_dn: list = [0.0, 0.0, 0.0],
        same_phase_up_dn: bool = False,
        neldiff: int = 0,
    ):
        """

        This class is a wrapper of pyturbo makefort10 class

        Attributes:
            structure_file (str): File name of the input structure, formats suppored by ASE are supported.
            supercell (list):  3 integers, supercell sizes [x,y,z]
            det_basis_set (str or list): basis set for the determinant part: e.g., "cc-pVQZ" (str), a list of gamess format basis sets is accepatable.
            jas_basis_set (str or list): basis set for the Jastrow part: e.g., "cc-pVQZ" (str), a list of gamess format basis sets is accepatable.
            det_contracted_flag (bool): if True determinant basis set is contracted, if False determinant basis set is uncontracted.
            jas_contracted_flag (bool): if True Jastrow basis set is contracted, if False Jastrow basis set is uncontracted.
            all_electron_jas_basis_set (bool): if True Jastrow basis set is read from the specified all-electron basis, if False, pseudo potential ones.
            pseudo_potential (str, list or None): if None, all-electron calculations, if "str", the corresponding PP is read from the database.
            det_cut_basis_option (bool): if True, determinant basis set is cut according to the Andrea Zen's procedure.
            jas_cut_basis_option (bool): if True, Jastrow basis set is cut according to the Andrea Zen's procedure.
            jastrow_type (int): One- and Two- Jastrow type specified.
            complex (bool): if True, the WF is complex, if False, the WF is real.
            phase_up (list): 3-float numbers for the up-phase [x, y, z].
            phase_dn (list): 3-float numbers for the dn-phase [x, y, z].
            same_phase_up_dn (bool): forced phase up == phase dn (valid only for gamma point.) it is automatically detected for other points.
            neldiff (int): The number of difference between up and dn electrons.
        """

        if os.path.isfile("fort.10"):
            logger.warning("fort.10 in the current directory is overwritten!")

        makefort10 = Makefort10_genius(
            structure_file=structure_file,
            supercell=supercell,
            det_basis_set=det_basis_set,
            jas_basis_set=jas_basis_set,
            det_contracted_flag=det_contracted_flag,
            jas_contracted_flag=jas_contracted_flag,
            all_electron_jas_basis_set=all_electron_jas_basis_set,
            pseudo_potential=pseudo_potential,
            det_cut_basis_option=det_cut_basis_option,
            jas_cut_basis_option=jas_cut_basis_option,
            jastrow_type=jastrow_type,
            complex=complex,
            phase_up=phase_up,
            phase_dn=phase_dn,
            same_phase_up_dn=same_phase_up_dn,
            neldiff=neldiff,
        )

        makefort10.run_all()
        if os.path.isfile("fort.10"):
            os.remove("fort.10")
        shutil.move("fort.10_new", "fort.10")

        self.io_fort10 = IO_fort10(fort10="fort.10")
        self.read_flag = True
        self.wf_format = "turborvb"
        self.trexio_filename = None

        # agp -> sd
        self.to_sd(grid_size=0.50, clean_flag=True)

        # sd
        logger.warning(
            "The generated fort.10 has random MO coefficients!! plz. initialize them by DFT."
        )

    # to pfaffian
    def to_pf(
        self,
        grid_size: float = 0.10,
        rotate_angle: float = -0.125,
        additional_hyb: list = [],
        nosym: bool = False,
        clean_flag: bool = False,
    ) -> None:
        """
        Convert for.10 to the Pfaffian format

        Args:
            grid_size (float): grid size (bohr)
            rotate_angle (float): angle for the rotation (/180 degree)
            additional_hyb (list): a list of the numbers of added hybrid orbitals
            nosym (bool): flag for nosymmetry
            clean_flag (bool): cleaning temporary files

        """

        if not self.read_flag:
            logger.warning(
                "WF file is not read yet! Please read a WF file first."
            )
            return

        # check if the wf is polarized or not.
        if self.io_fort10.f10header.nelup == self.io_fort10.f10header.neldn:
            logger.info("non-spin polarized case")
            self.io_fort10.to_agp(
                triplet=True,
                pfaffian_flag=True,
                grid_size=grid_size,
                additional_hyb=additional_hyb,
                nosym=nosym,
                clean_flag=clean_flag,
                only_generate_template=True,
            )
            shutil.copy(self.io_fort10.fort10, "fort.10_bak")
            shutil.copy(self.io_fort10.fort10, "fort.10_in")
            convertpfaff_genius = Convertpfaff_genius(
                in_fort10="fort.10_in", out_fort10="fort.10_out"
            )
            logger.info("convert to pf w/o rotation.")
            convertpfaff_genius.run(rotate_flag=False)
            shutil.move("fort.10_new", "fort.10")
            shutil.copy("fort.10_in", "fort.10_new")
            copy_jastrow(fort10_to="fort.10", fort10_from="fort.10_new")
            shutil.copy("fort.10", "fort.10_in")
            shutil.copy("fort.10", "fort.10_out")
            logger.info("convert to pf w rotation.")
            convertpfaff_genius.run(
                rotate_flag=True, rotate_angle=rotate_angle
            )
            shutil.move("fort.10_new", "fort.10")
            shutil.copy("fort.10_in", "fort.10_new")
            copy_jastrow(fort10_to="fort.10", fort10_from="fort.10_new")

            if clean_flag:
                os.remove("fort.10_new")
                os.remove("fort.10_in")
                os.remove("fort.10_out")

        else:
            logger.info("spin polarized case")
            shutil.copy(self.io_fort10.fort10, "fort.10_bak")
            shutil.copy(self.io_fort10.fort10, "fort.10_in")
            self.io_fort10.to_agp(
                triplet=True,
                pfaffian_flag=True,
                grid_size=grid_size,
                additional_hyb=additional_hyb,
                nosym=nosym,
                clean_flag=clean_flag,
                only_generate_template=True,
            )
            convertpfaff_genius = Convertpfaff_genius(
                in_fort10="fort.10_in", out_fort10="fort.10_out"
            )
            convertpfaff_genius.run(rotate_flag=False)
            shutil.move("fort.10_new", "fort.10")
            shutil.copy("fort.10_in", "fort.10_new")
            copy_jastrow(fort10_to="fort.10", fort10_from="fort.10_new")
            if clean_flag:
                os.remove("fort.10_new")
                os.remove("fort.10_in")
                os.remove("fort.10_out")

    # to sd
    def to_sd(self, grid_size: float = 0.10, clean_flag: bool = False) -> None:
        """
        Convert for.10 to the Slater determinant format

        Args:
            grid_size (float): grid size (bohr)
            clean_flag (bool): cleaning temporary files

        """
        if not self.read_flag:
            logger.warning(
                "WF file is not read yet! Please read a WF file first."
            )
            return

        # mo=Ndn/2
        self.add_MOs(
            add_random_mo=False,
            grid_size=grid_size,
            additional_mo=0,
            clean_flag=clean_flag,
        )

    # to agps
    def to_agps(
        self,
        grid_size: float = 0.10,
        additional_hyb: list = [],
        nosym: float = False,
        clean_flag: float = False,
    ) -> None:
        """
        Convert for.10 to the symmetric AGP format

        Args:
            grid_size (float): grid size (bohr)
            additional_hyb (list): a list of the numbers of added hybrid orbitals
            nosym (bool): flag for nosymmetry
            clean_flag (bool): cleaning temporary files

        """
        if not self.read_flag:
            logger.warning(
                "WF file is not read yet! Please read a WF file first."
            )
            return

        self.to_agp(
            triplet=False,
            grid_size=grid_size,
            additional_hyb=additional_hyb,
            pfaffian_flag=False,
            nosym=nosym,
            clean_flag=clean_flag,
        )

    # to agpu
    def to_agpu(
        self,
        grid_size: float = 0.10,
        additional_hyb: list = [],
        nosym: bool = False,
        clean_flag: bool = False,
    ) -> None:
        """
        Convert for.10 to the non-symetric AGP format

        Args:
            grid_size (float): grid size (bohr)
            additional_hyb (list): a list of the numbers of added hybrid orbitals
            nosym (bool): flag for nosymmetry
            clean_flag (bool): cleaning temporary files

        """
        if not self.read_flag:
            logger.warning(
                "WF file is not read yet! Please read a WF file first."
            )
            return

        self.to_agp(
            triplet=True,
            grid_size=grid_size,
            additional_hyb=additional_hyb,
            pfaffian_flag=False,
            nosym=nosym,
            clean_flag=clean_flag,
        )

    # to agps(triplet=false) or agpu(triplet=true)
    def to_agp(
        self,
        triplet: bool = False,
        pfaffian_flag: bool = False,
        grid_size: float = 0.10,
        additional_hyb: list = [],
        nosym: bool = False,
        clean_flag: bool = False,
        only_generate_template: bool = False,
    ) -> None:
        """
        Convert for.10 to the symmetric or non-symmetric AGP format

        Args:
            triplet (bool): flag for including the triplet pairing function, i.e., agps(triplet=false) or agpu(triplet=true)
            pfaffian_flag (bool): Pfaffian
            grid_size (float): grid size (bohr)
            additional_hyb (list): a list of the numbers of added hybrid orbitals
            nosym (bool): flag for nosymmetry
            clean_flag (bool): cleaning temporary files

        """
        if not self.read_flag:
            logger.warning(
                "WF file is not read yet! Please read a WF file first."
            )
            return

        shutil.copy(self.io_fort10.fort10, "fort.10_in")

        if pfaffian_flag:
            logger.info("Generating pfaffian ansatz")
            logger.info("triplet is set to True")
            triplet = True
        # prepare makefort10
        # - to do - add jastrow, etc...
        # makefort10 class
        structure = self.io_fort10.f10structure.structure
        jastrow_type = self.io_fort10.f10header.jas_2body
        twobody_list = self.io_fort10.f10jastwobody.twobody_list
        onebody_list = self.io_fort10.f10jastwobody.onebody_list
        if triplet:
            logger.info("AGPu anstaz!!")
            logger.info("spin-depenedent jastrow is recommended.")
            if jastrow_type in {0}:
                pass
            elif jastrow_type in {-15}:
                jastrow_type = -22
                logger.info(f"jastrow_type is set to {jastrow_type}")
            elif jastrow_type in {-5, -6}:
                jastrow_type = -26
                logger.info(f"jastrow_type is set to {jastrow_type}")
            elif jastrow_type in {-22, -26, -27}:
                pass
            else:
                logger.error(
                    f"jastrow_type = {jastrow_type} is not supported. only 0,-5,-6,-15,-22,-26,-27 are expected."
                )
        else:
            logger.info("AGPs anstaz!!")
        namelist = Makefort10.read_default_namelist(
            structure=structure, jastrow_type=jastrow_type
        )
        for i, twobody in enumerate(twobody_list):
            namelist.set_parameter(
                f"twobodypar({i + 1})", twobody, "&electrons"
            )
        for i, onebody in enumerate(onebody_list):
            namelist.set_parameter(
                f"onebodypar({i + 1})", onebody, "&electrons"
            )
        namelist.comment_out("twobodypar")
        namelist.comment_out("onebodypar")
        if triplet:  # jagpu
            namelist.set_parameter("symmagp", ".false.", "&symmetries")
        else:  # jagps
            namelist.set_parameter("symmagp", ".true.", "&symmetries")
        if self.io_fort10.pp_flag:
            pseudo_potentials = Pseudopotentials.parse_pseudopotential_from_turborvb_pseudo_dat(
                file="pseudo.dat"
            )
            # compute z_core
            atomic_numbers = self.io_fort10.f10structure.atomic_numbers
            valence_electrons = self.io_fort10.f10structure.valence_electrons
            z_core = list(
                np.array(atomic_numbers) - np.array(valence_electrons)
            )
            pseudo_potentials.z_core = z_core
        else:
            pseudo_potentials = Pseudopotentials()
        det_basis_sets = self.io_fort10.f10detbasissets.det_basis_sets
        jas_basis_sets = self.io_fort10.f10jasbasissets.jas_basis_sets

        # add number of hybrid orbitals
        if len(additional_hyb) != 0:
            assert len(additional_hyb) == self.io_fort10.f10header.natom
            det_basis_sets.number_of_additional_hybrid_orbitals = (
                additional_hyb
            )

        if self.io_fort10.pbc_flag:
            (
                phase_up_1,
                phase_up_2,
                phase_up_3,
            ) = self.io_fort10.f10structure.phase_up
            (
                phase_dn_1,
                phase_dn_2,
                phase_dn_3,
            ) = self.io_fort10.f10structure.phase_dn
            namelist.set_parameter(
                parameter="phase(1)", value=phase_up_1, namelist="&system"
            )
            namelist.set_parameter(
                parameter="phase(2)", value=phase_up_2, namelist="&system"
            )
            namelist.set_parameter(
                parameter="phase(3)", value=phase_up_3, namelist="&system"
            )
            namelist.set_parameter(
                parameter="phasedo(1)", value=phase_dn_1, namelist="&system"
            )
            namelist.set_parameter(
                parameter="phasedo(2)", value=phase_dn_2, namelist="&system"
            )
            namelist.set_parameter(
                parameter="phasedo(3)", value=phase_dn_3, namelist="&system"
            )

        if self.io_fort10.f10header.complex_flag:
            namelist.set_parameter(
                parameter="complexfort10", value=".true.", namelist="&system"
            )

        if pfaffian_flag:
            namelist.set_parameter(
                parameter="yes_pfaff", value=".true.", namelist="&system"
            )

        if nosym:
            namelist.set_parameter(
                parameter="nosym", value=".true.", namelist="&symmetries"
            )

        # neldiff
        neldiff = (
            self.io_fort10.f10header.nelup - self.io_fort10.f10header.neldn
        )
        namelist.set_parameter("neldiff", neldiff, "&electrons")

        makefort10 = Makefort10(
            structure=structure,
            det_basis_sets=det_basis_sets,
            jas_basis_sets=jas_basis_sets,
            pseudo_potentials=pseudo_potentials,
            namelist=namelist,
        )

        makefort10.generate_input(input_name="makefort10.input")
        makefort10.run(input_name="makefort10.input", output_name="out_make")
        shutil.move("fort.10_new", "fort.10_out")

        if only_generate_template:
            logger.warning("A template AGP file, fort.10_out, is generated.")
        else:
            # convertfort10
            convertfort10_genius = Convertfort10_genius(
                in_fort10="fort.10_in",
                out_fort10="fort.10_out",
                grid_size=grid_size,
            )

            convertfort10_genius.generate_input(
                input_name="convertfort10.input"
            )
            convertfort10_genius.run(
                input_name="convertfort10.input", output_name="out_conv"
            )
            shutil.move(self.io_fort10.fort10, "fort.10_bak")
            shutil.move("fort.10_new", "fort.10")
            shutil.copy("fort.10_in", "fort.10_new")
            copy_jastrow(fort10_to="fort.10", fort10_from="fort.10_new")

        if clean_flag:
            os.remove("fort.10_new")
            os.remove("fort.10_in")
            os.remove("fort.10_out")

    def add_MOs(
        self,
        add_random_mo: bool = True,
        grid_size: float = 0.10,
        additional_mo: int = 0,
        clean_flag: bool = False,
    ) -> None:
        """
        Add molecular orbitals to fort.10

        Args:
            add_random_mo (bool): flag for randomized MOs)
            grid_size (float): grid size (bohr)
            additional_mo (int): the number of added MOs (beyond the HOMO)
            clean_flag (bool): cleaning temporary files

        """
        shutil.copy(self.io_fort10.fort10, "fort.10_in")
        convertfort10mol_genius = Convertfort10mol_genius(
            fort10="fort.10_in",
            add_random_mo=add_random_mo,
            grid_size=grid_size,
            additional_mo=additional_mo,
        )
        convertfort10mol_genius.generate_input(
            input_name="convertfort10mol.input"
        )
        convertfort10mol_genius.run(
            input_name="convertfort10mol.input", output_name="out_mol"
        )

        if clean_flag:
            remove_file("fort.10")
            shutil.move("fort.10_new", "fort.10")
            remove_file("fort.10_in")
            remove_file("convertfort10.input")
            remove_file("out_conv")

    def write_structure(self, structure_file: str):
        """
        Write a structure file (i.e., convert fort.10 to a structure file, e.g., xyz)

        Args:
            structure_file (str): structure file name. all formats supported by ASE are acceptable.

        """
        structure = self.io_fort10.f10structure.structure
        structure.write(file=structure_file)


if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    logger_p = getLogger("pyturbo")
    logger_p.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)
    logger_p.addHandler(stream_handler)

    os.chdir(
        os.path.join(
            turbo_genius_root, "tests", "turbogenius_test", "wavefunction_test"
        )
    )

    # moved to examples
    wavefunction = Wavefunction()
    # wavefunction.read_from_fort10(
    #    fort10=os.path.join("fort10_examples", "fort.10_AcOH")
    # )
    # wavefunction.write_structure(structure_file="test.xyz")
    wavefunction.read_from_structure(structure_file="test.xyz")
    wavefunction.write_structure(structure_file="test2.xyz")
    wavefunction.read_from_trexio(
        trexio_filename=os.path.join("fort10_examples", "methane_trexio.hdf5")
    )

    # wavefunction.add_MOs()
    # wavefunction.to_agp()
    # wavefunction.to_sd()
    # wavefunction.to_pf()
