#!python
# -*- coding: utf-8 -*-

#python modules
from __future__ import annotations
import os, sys
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import click
import glob
import pickle

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pyturbo.io_fort10 import IO_fort10
from pyturbo.pseudopotentials import Pseudopotentials
from pyturbo.makefort10 import Makefort10
from pyturbo.utils.utility import remove_file
#turbo-genius modules
from utils_workflows.env import turbo_genius_root, turbo_genius_data_dir
from convertfort10mol_genius import Convertfort10mol_genius
from convertfort10_genius import Convertfort10_genius
from convertpfaff_genius import  Convertpfaff_genius
from tools_genius import copy_jastrow

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)
#logger = getLogger(__name__)

from turbo_genius_cli import cli, decorate_grpost, header, OptionEatAll
@cli.command(short_help = "convert wavefunction")
@decorate_grpost
@click.option("-to", "to_ansatz",
              help= 'Specify to ansatz',
              default = "agps",
              type = str)
@click.option("-hyb", "hybrid_orbitals",
              help= 'Specify the number of added hybrid orbitals for each atom (0, 5, ...)',
              default=(),
              type = list,
              cls=OptionEatAll
              )
@click.option("-nosym", "nosymmetry",
              help= 'flag for nosymmetry)',
              is_flag=True,
              default=False,
              )
@click.option("-rot", "rotate_angle",
              help= 'Specify rotate angle (-0.5 - + 0.5)',
              default = 0.00,
              type = float)
@click.option("-grid", "grid_size",
              help= 'Specify grid size',
              default = 0.10,
              type = float)
@header
def convertwf(
            g,r,post,
            operation,
            log_level,
            to_ansatz,
            hybrid_orbitals,
            grid_size,
            nosymmetry,
            rotate_angle,
):
    pkl_name="wavefunction_cli.pkl"
    root_dir=os.getcwd()
    number_of_additional_hybrid_orbitals=list(map(int,hybrid_orbitals))

    wavefunction = Wavefunction(fort10='fort.10')

    if to_ansatz=='agps':
        logger.info("convert the wf to agps")
        wavefunction.to_agps(grid_size=grid_size, additional_hyb=number_of_additional_hybrid_orbitals, nosym=nosymmetry, clean_flag=False)

    elif to_ansatz=='agpu':
        logger.info("convert the wf to agpu")
        wavefunction.to_agpu(grid_size=grid_size, additional_hyb=number_of_additional_hybrid_orbitals, nosym=nosymmetry, clean_flag=False)

    elif to_ansatz == 'sd':
        logger.info("convert the wf to jsd (i.e. adding MOs)")
        wavefunction.to_sd(grid_size=grid_size, clean_flag=False)

    elif to_ansatz == 'pf':
        logger.info("convert the wf to pfaffian")
        logger.info(f"rotate_angle = {rotate_angle}")
        wavefunction.to_pf(grid_size=grid_size, rotate_angle=rotate_angle, nosym=nosymmetry, clean_flag=False)

class Wavefunction(IO_fort10):

    def __init__(self, fort10):
        super().__init__(fort10=fort10)

    # to pfaffian
    def to_pf(self, grid_size=0.10, rotate_angle=-0.125, additional_hyb=[], nosym=False, clean_flag=False):
        # check if the wf is polarized or not.
        if self.f10header.nelup == self.f10header.neldn:
            logger.info("non-spin polarized case")
            self.to_agp(triplet=True, pfaffian_flag=True, grid_size=grid_size, additional_hyb=additional_hyb, nosym=nosym, clean_flag=clean_flag, only_generate_template=True)
            shutil.copy(self.fort10, "fort.10_bak")
            shutil.copy(self.fort10, "fort.10_in")
            convertpfaff_genius = Convertpfaff_genius(in_fort10="fort.10_in", out_fort10="fort.10_out")
            logger.info("convert to pf w/o rotation.")
            convertpfaff_genius.run(rotate_flag=False)
            shutil.move("fort.10_new", "fort.10")
            shutil.copy("fort.10_in", "fort.10_new")
            copy_jastrow(fort10_to='fort.10', fort10_from='fort.10_new')
            shutil.copy("fort.10", "fort.10_in")
            shutil.copy("fort.10", "fort.10_out")
            logger.info("convert to pf w rotation.")
            convertpfaff_genius.run(rotate_flag=True, rotate_angle=rotate_angle)
            shutil.move("fort.10_new", "fort.10")
            shutil.copy("fort.10_in", "fort.10_new")
            copy_jastrow(fort10_to='fort.10', fort10_from='fort.10_new')

            if clean_flag:
                os.remove("fort.10_new")
                os.remove("fort.10_in")
                os.remove("fort.10_out")

        else:
            logger.info("spin polarized case")
            shutil.copy(self.fort10, "fort.10_bak")
            shutil.copy(self.fort10, "fort.10_in")
            self.to_agp(triplet=True, pfaffian_flag=True, grid_size=grid_size, additional_hyb=additional_hyb, nosym=nosym, clean_flag=clean_flag, only_generate_template=True)
            convertpfaff_genius = Convertpfaff_genius(in_fort10="fort.10_in", out_fort10="fort.10_out")
            convertpfaff_genius.run(rotate_flag=False)
            shutil.move("fort.10_new", "fort.10")
            shutil.copy("fort.10_in", "fort.10_new")
            copy_jastrow(fort10_to='fort.10', fort10_from='fort.10_new')
            if clean_flag:
                os.remove("fort.10_new")
                os.remove("fort.10_in")
                os.remove("fort.10_out")

    # to sd
    def to_sd(self, grid_size=0.10, clean_flag=False):
        # mo=Ndn/2
        self.add_MOs(add_random_mo=False, grid_size=grid_size, additional_mo=0, clean_flag=clean_flag)

    # to agps
    def to_agps(self, grid_size=0.10, additional_hyb=[], nosym=False, clean_flag=False):
        self.to_agp(triplet=False, grid_size=grid_size, additional_hyb=additional_hyb, pfaffian_flag=False, nosym=nosym, clean_flag=clean_flag)

    # to agpu
    def to_agpu(self, grid_size=0.10, additional_hyb=[], nosym=False, clean_flag=False):
        self.to_agp(triplet=True, grid_size=grid_size, additional_hyb=additional_hyb, pfaffian_flag=False, nosym=nosym, clean_flag=clean_flag)

    # to agps(triplet=false) or agpu(triplet=true)
    def to_agp(self, triplet=False, pfaffian_flag=False, grid_size=0.10, additional_hyb=[], nosym=False, clean_flag=False, only_generate_template=False):
        shutil.copy(self.fort10, "fort.10_in")

        if pfaffian_flag:
            logger.info("Generating pfaffian ansatz")
            logger.info(f"triplet is set to True")
            triplet = True
        #prepare makefort10
        # - to do - add jastrow, etc...
        # makefort10 class
        structure = self.f10structure.structure
        jastrow_type = self.f10header.jas_2body
        twobody_list = self.f10jastwobody.twobody_list
        onebody_list = self.f10jastwobody.onebody_list
        if triplet:
            logger.info("AGPu anstaz!!")
            logger.info("spin-depenedent jastrow is recommended.")
            if jastrow_type in {0}:
                pass
            elif jastrow_type in {-15}:
                jastrow_type=-22
                logger.info(f"jastrow_type is set to {jastrow_type}")
            elif jastrow_type in {-5,-6}:
                jastrow_type=-26
                logger.info(f"jastrow_type is set to {jastrow_type}")
            elif jastrow_type in {-22,-26,-27}:
                pass
            else:
                logger.error(f"jastrow_type = {jastrow_type} is not supported. only 0,-5,-6,-15,-22,-26,-27 are expected.")
        else:
            logger.info("AGPs anstaz!!")
        namelist = Makefort10.read_default_namelist(structure=structure, jastrow_type=jastrow_type)
        for i, twobody in enumerate(twobody_list):
            namelist.set_parameter(f"twobodypar({i + 1})", twobody, '&electrons')
        for i, onebody in enumerate(onebody_list):
            namelist.set_parameter(f"onebodypar({i + 1})", onebody, '&electrons')
        namelist.comment_out("twobodypar")
        namelist.comment_out("onebodypar")
        if triplet: # jagpu
            namelist.set_parameter(f"symmagp", '.false.', '&symmetries')
        else: # jagps
            namelist.set_parameter(f"symmagp", '.true.', '&symmetries')
        if self.pp_flag:
            pseudo_potentials = Pseudopotentials.parse_pseudopotential_from_turborvb_pseudo_dat(file="pseudo.dat")
            # compute z_core
            atomic_numbers = self.f10structure.atomic_numbers
            valence_electrons = self.f10structure.valence_electrons
            z_core = list(np.array(atomic_numbers) - np.array(valence_electrons))
            pseudo_potentials.z_core = z_core
        else:
            pseudo_potentials = Pseudopotentials()
        det_basis_sets = self.f10detbasissets.det_basis_sets
        jas_basis_sets = self.f10jasbasissets.jas_basis_sets

        # add number of hybrid orbitals
        if len(additional_hyb) != 0:
            assert len(additional_hyb) == self.f10header.natom
            det_basis_sets.number_of_additional_hybrid_orbitals=additional_hyb

        if self.pbc_flag:
            phase_up_1, phase_up_2, phase_up_3 = self.f10structure.phase_up
            phase_dn_1, phase_dn_2, phase_dn_3 = self.f10structure.phase_dn
            namelist.set_parameter(parameter="phase(1)", value=phase_up_1, namelist="&system")
            namelist.set_parameter(parameter="phase(2)", value=phase_up_2, namelist="&system")
            namelist.set_parameter(parameter="phase(3)", value=phase_up_3, namelist="&system")
            namelist.set_parameter(parameter="phasedo(1)", value=phase_dn_1, namelist="&system")
            namelist.set_parameter(parameter="phasedo(2)", value=phase_dn_2, namelist="&system")
            namelist.set_parameter(parameter="phasedo(3)", value=phase_dn_3, namelist="&system")

        if self.f10header.complex_flag:
            namelist.set_parameter(parameter="complexfort10", value='.true.', namelist="&system")

        if pfaffian_flag:
            namelist.set_parameter(parameter="yes_pfaff", value='.true.', namelist="&system")

        if nosym:
            namelist.set_parameter(parameter="nosym", value='.true.', namelist="&symmetries")

        #neldiff
        neldiff = self.f10header.nelup - self.f10header.neldn
        namelist.set_parameter(f"neldiff", neldiff, '&electrons')

        makefort10 = Makefort10(
            structure=structure,
            det_basis_sets=det_basis_sets,
            jas_basis_sets=jas_basis_sets,
            pseudo_potentials=pseudo_potentials,
            namelist=namelist
        )

        makefort10.generate_input(input_name='makefort10.input')
        makefort10.run(input_name='makefort10.input', output_name='out_make')
        shutil.move("fort.10_new", "fort.10_out")

        if only_generate_template:
            return

        #convertfort10
        convertfort10_genius=Convertfort10_genius(
                                    in_fort10="fort.10_in",
                                    out_fort10="fort.10_out",
                                    grid_size=grid_size
                                    )

        convertfort10_genius.generate_input(input_name="convertfort10.input")
        convertfort10_genius.run(input_name='convertfort10.input', output_name='out_conv')
        shutil.move(self.fort10, "fort.10_bak")
        shutil.move("fort.10_new", "fort.10")
        shutil.copy("fort.10_in", "fort.10_new")
        copy_jastrow(fort10_to='fort.10', fort10_from='fort.10_new')

        if clean_flag:
            os.remove("fort.10_new")
            os.remove("fort.10_in")
            os.remove("fort.10_out")

    def add_MOs(self, add_random_mo=True, grid_size=0.10, additional_mo=0, clean_flag=False):
        shutil.copy(self.fort10, "fort.10_in")
        convertfort10mol_genius=Convertfort10mol_genius(fort10="fort.10_in",
                                                        add_random_mo=add_random_mo,
                                                        grid_size=grid_size,
                                                        additional_mo=additional_mo)
        convertfort10mol_genius.generate_input(input_name='convertfort10mol.input')
        convertfort10mol_genius.run(input_name='convertfort10mol.input', output_name='out_mol')

        if clean_flag:
            remove_file('fort.10')
            shutil.move('fort.10_new', 'fort.10')
            remove_file('fort.10_in')
            remove_file('convertfort10.input')
            remove_file('out_conv')

    def write(self, file):
        structure=self.f10structure.structure
        structure.write(file=file)

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    logger_p = getLogger("pyturbo")
    logger_p.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)
    logger_p.addHandler(stream_handler)

    from utils_workflows.env import turbo_genius_root

    os.chdir(os.path.join(turbo_genius_root, "tests", "wavefunction"))

    # moved to examples
    wavefunction=Wavefunction(fort10="fort.10_agpu")
    wavefunction.write(file="test.xyz")
    #wavefunction.add_MOs()
    #wavefunction.to_agp()
    #wavefunction.to_sd()
    wavefunction.to_pf()