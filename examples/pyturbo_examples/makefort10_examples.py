#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys, shutil

#set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo')
logger.setLevel("DEBUG")
stream_handler = StreamHandler()
stream_handler.setLevel("DEBUG")
handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)

#pyturbo modules
from turbogenius.pyturbo.makefort10 import Makefort10
from turbogenius.pyturbo.structure import Structure
from turbogenius.pyturbo.pseudopotentials import Pseudopotentials
from turbogenius.pyturbo.basis_set import Jas_Basis_sets, Det_Basis_sets
from turbogenius.pyturbo.utils.downloader import BSE, BFD, ccECP
from turbogenius.utils_workflows.env import turbo_genius_root

prefix="makefort10"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))
makefort10_test_dir=os.path.join(example_root_dir, prefix)

#test ok
# Hydrogen dimer, ccECP
str_input = "H2.xyz"
ccecp=ccECP(basis_sets_output_dir=makefort10_test_dir, pseudo_potential_output_dir=makefort10_test_dir)
if not os.path.isfile("H_ccECP.pseudo"): ccecp.to_file(element_list=["H"], basis_list="cc-pVDZ")
pp_files = ["H_ccECP.pseudo", "H_ccECP.pseudo"]
basis_file="H-ccECP-cc-pVDZ" # gamess format, downloaded later, see below
basis_files=[basis_file, basis_file]

# Benzene, all-electrons and pseudo potentials, mixed.
# error!! new option of makefort10 up to second digit -> up to third/fourth digit.
str_input = "benzene.xyz"
bfd = BFD(basis_sets_output_dir=makefort10_test_dir, pseudo_potential_output_dir=makefort10_test_dir) # basis sets and pseudo potentials are downloaded
if not os.path.isfile("C_BFD.pseudo"): bfd.to_file(element_list=["C"], basis_list=["vqz"])
if not os.path.isfile("H_BFD.pseudo"): bfd.to_file(element_list=["H"], basis_list=["vqz"])
C_pp_file="C_BFD.pseudo"; C_basis_file="C_vqz.basis"
H_pp_file=None; H_basis_file="H_vqz.basis"
pp_files = [C_pp_file] * 6 + [H_pp_file] * 6
basis_files = [C_basis_file] * 6 + [H_basis_file] * 6
structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(structure=structure, jastrow_type=-6)
pseudo_potentials = Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist
)
makefort10.sanity_check()
makefort10.generate_input(input_name="makefort10.input", basis_sets_unique_element=False)
makefort10.set_parameter(parameter="nosym", value='.true', namelist="&symmetries")
logger.info(makefort10.get_parameter(parameter="nosym"))
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])



# Diamond, PP, orthorhombic case
str_input = "Diamond.cif"
bfd = BFD(basis_sets_output_dir=makefort10_test_dir, pseudo_potential_output_dir=makefort10_test_dir) # basis sets and pseudo potentials are downloaded
if not os.path.isfile("C_ecp"): bfd.to_file(element_list=["C"], basis_list=["vqz"])
pp_files = ["C_BFD.pseudo"] * 8
basis_files=["C_vqz.basis"] * 8
structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(structure=structure, jastrow_type=-6)
pseudo_potentials = Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist
)
makefort10.sanity_check()
makefort10.generate_input(input_name="makefort10.input", basis_sets_unique_element=False)
makefort10.set_parameter(parameter="nosym", value='.true', namelist="&symmetries")
logger.info(makefort10.get_parameter(parameter="nosym"))
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])



# Diamond, PP, orthorhombic case
str_input = "Diamond.cif"
bfd = BFD(basis_sets_output_dir=makefort10_test_dir, pseudo_potential_output_dir=makefort10_test_dir) # basis sets and pseudo potentials are downloaded
if not os.path.isfile("C_BFD.pseudo"): bfd.to_file(element_list=["C"], basis_list=["vqz"])
pp_files = [None] * 8
basis_files=["C_vqz.basis"] * 8
structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(structure=structure, jastrow_type=-6)
pseudo_potentials = Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist
)
makefort10.sanity_check()
makefort10.generate_input(input_name="makefort10.input", basis_sets_unique_element=False)
makefort10.set_parameter(parameter="nosym", value='.true', namelist="&symmetries")
logger.info(makefort10.get_parameter(parameter="nosym"))
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])


# SiO2, pseudo potentials, non-orthorhombic case
str_input = "silicon_oxide.cif"
structure = Structure.parse_structure_from_file(file=str_input)
bfd = BFD(basis_sets_output_dir=makefort10_test_dir, pseudo_potential_output_dir=makefort10_test_dir) # basis sets and pseudo potentials are downloaded
if not os.path.isfile("Si_BFD.pseudo"): bfd.to_file(element_list=["Si"], basis_list=["vtz"])
if not os.path.isfile("O_BFD.pseudo"): bfd.to_file(element_list=["O"], basis_list=["vtz"])
pp_files = ["Si_BFD.pseudo"] * 3 + ["O_BFD.pseudo"] * 6
basis_files=["Si_vtz.basis"] * 3 + ["O_vtz.basis"] * 6

structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(structure=structure, jastrow_type=-6)
pseudo_potentials = Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(files=basis_files)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist
)
makefort10.sanity_check()
makefort10.generate_input(input_name="makefort10.input", basis_sets_unique_element=False)
makefort10.set_parameter(parameter="nosym", value='.true', namelist="&symmetries")
logger.info(makefort10.get_parameter(parameter="nosym"))
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])
