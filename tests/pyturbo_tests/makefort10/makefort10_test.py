#!python
# -*- coding: utf-8 -*-
import os, sys
import shutil

# pyturbo modules
from turbogenius.pyturbo.makefort10 import Makefort10
from turbogenius.pyturbo.structure import Structure
from turbogenius.pyturbo.pseudopotentials import Pseudopotentials
from turbogenius.pyturbo.basis_set import Jas_Basis_sets, Det_Basis_sets
from turbogenius.pyturbo.utils.downloader import BSE, BFD, ccECP

root_dir = os.path.dirname(__file__)
os.chdir(root_dir)

# Hydrogen dimer, ccECP
str_input = "H2.xyz"
prefix, _ = os.path.splitext(str_input)
if os.path.isdir(os.path.join(root_dir, prefix)):
    shutil.rmtree(os.path.join(root_dir, prefix))
os.makedirs(os.path.join(root_dir, prefix))
shutil.copy(
    os.path.join(root_dir, str_input),
    os.path.join(root_dir, prefix, str_input),
)
os.chdir(os.path.join(root_dir, prefix))
ccecp = ccECP(
    basis_sets_output_dir=os.path.join(root_dir, prefix),
    pseudo_potential_output_dir=os.path.join(root_dir, prefix),
)
if not os.path.isfile("H_ccECP.pseudo"):
    ccecp.to_file(element_list=["H"], basis_list="cc-pVDZ")
pp_files = ["H_ccECP.pseudo", "H_ccECP.pseudo"]
basis_file = "H_cc-pVDZ.basis"
basis_files = [basis_file, basis_file]
structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(
    structure=structure, jastrow_type=-6
)
pseudo_potentials = (
    Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist,
)
assert makefort10.get_parameter(parameter="nosym") == ".false."
makefort10.set_parameter(
    parameter="nosym", value=".true.", namelist="&symmetries"
)
assert makefort10.get_parameter(parameter="nosym") == ".true."
makefort10.generate_input(
    input_name="makefort10.input", basis_sets_unique_element=False
)
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])
assert flags
os.chdir(root_dir)

# Benzene, all-electrons and pseudo potentials, mixed.
str_input = "benzene.xyz"
prefix, _ = os.path.splitext(str_input)
if os.path.isdir(os.path.join(root_dir, prefix)):
    shutil.rmtree(os.path.join(root_dir, prefix))
os.makedirs(os.path.join(root_dir, prefix))
shutil.copy(
    os.path.join(root_dir, str_input),
    os.path.join(root_dir, prefix, str_input),
)
os.chdir(os.path.join(root_dir, prefix))
bfd = BFD(
    basis_sets_output_dir=os.path.join(root_dir, prefix),
    pseudo_potential_output_dir=os.path.join(root_dir, prefix),
)  # basis sets and pseudo potentials are downloaded
if not os.path.isfile("C_BFD.pseudo"):
    bfd.to_file(element_list=["C"], basis_list=["vqz"])
if not os.path.isfile("H_BFD.pseudo"):
    bfd.to_file(element_list=["H"], basis_list=["vqz"])
C_pp_file = "C_BFD.pseudo"
C_basis_file = "C_vqz.basis"
H_pp_file = None
H_basis_file = "H_vqz.basis"
pp_files = [C_pp_file] * 6 + [H_pp_file] * 6
basis_files = [C_basis_file] * 6 + [H_basis_file] * 6
structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(
    structure=structure, jastrow_type=-6
)
pseudo_potentials = (
    Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist,
)
makefort10.generate_input(
    input_name="makefort10.input", basis_sets_unique_element=False
)
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])
assert flags
os.chdir(root_dir)

# Diamond, PP, orthorhombic case
str_input = "Diamond.cif"
prefix, _ = os.path.splitext(str_input)
if os.path.isdir(os.path.join(root_dir, prefix)):
    shutil.rmtree(os.path.join(root_dir, prefix))
os.makedirs(os.path.join(root_dir, prefix))
shutil.copy(
    os.path.join(root_dir, str_input),
    os.path.join(root_dir, prefix, str_input),
)
os.chdir(os.path.join(root_dir, prefix))
bfd = BFD(
    basis_sets_output_dir=os.path.join(root_dir, prefix),
    pseudo_potential_output_dir=os.path.join(root_dir, prefix),
)  # basis sets and pseudo potentials are downloaded
if not os.path.isfile("C_ecp"):
    bfd.to_file(element_list=["C"], basis_list=["vqz"])
pp_files = ["C_BFD.pseudo"] * 8
basis_files = ["C_vqz.basis"] * 8
structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(
    structure=structure, jastrow_type=-6
)
pseudo_potentials = (
    Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist,
)
makefort10.generate_input(
    input_name="makefort10.input", basis_sets_unique_element=False
)
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])
assert flags
os.chdir(root_dir)

# SiO2, pseudo potentials, non-orthorhombic case
str_input = "silicon_oxide.cif"
prefix, _ = os.path.splitext(str_input)
if os.path.isdir(os.path.join(root_dir, prefix)):
    shutil.rmtree(os.path.join(root_dir, prefix))
os.makedirs(os.path.join(root_dir, prefix))
shutil.copy(
    os.path.join(root_dir, str_input),
    os.path.join(root_dir, prefix, str_input),
)
os.chdir(os.path.join(root_dir, prefix))
structure = Structure.parse_structure_from_file(file=str_input)
bfd = BFD(
    basis_sets_output_dir=os.path.join(root_dir, prefix),
    pseudo_potential_output_dir=os.path.join(root_dir, prefix),
)  # basis sets and pseudo potentials are downloaded
if not os.path.isfile("Si_BFD.pseudo"):
    bfd.to_file(element_list=["Si"], basis_list=["vtz"])
if not os.path.isfile("O_BFD.pseudo"):
    bfd.to_file(element_list=["O"], basis_list=["vtz"])
pp_files = ["Si_BFD.pseudo"] * 3 + ["O_BFD.pseudo"] * 6
basis_files = ["Si_vtz.basis"] * 3 + ["O_vtz.basis"] * 6

structure = Structure.parse_structure_from_file(file=str_input)
namelist = Makefort10.read_default_namelist(
    structure=structure, jastrow_type=-6
)
pseudo_potentials = (
    Pseudopotentials.parse_pseudopotential_from_gamess_format_files(pp_files)
)
det_basis_sets = Det_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(
    files=basis_files
)
makefort10 = Makefort10(
    structure=structure,
    det_basis_sets=det_basis_sets,
    jas_basis_sets=jas_basis_sets,
    pseudo_potentials=pseudo_potentials,
    namelist=namelist,
)
makefort10.generate_input(
    input_name="makefort10.input", basis_sets_unique_element=False
)
makefort10.run(input_name="makefort10.input", output_name="out_make")
flags = makefort10.check_results(output_names=["out_make"])
assert flags
os.chdir(root_dir)
