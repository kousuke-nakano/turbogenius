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
from turbogenius.pyturbo.pseudopotentials import Pseudopotentials
from turbogenius.pyturbo.utils.downloader import BFD, ccECP
from turbogenius.utils_workflows.env import turbo_genius_root

#basis sets
prefix="pseudo_potentials_"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))
pseudopotentials_test_dir=os.path.join(example_root_dir, prefix)

# BFD pseudo potential
bfd=BFD(pseudo_potential_output_dir=pseudopotentials_test_dir)
bfd.to_file(element_list=["C"], basis_list=["vtz"])

# ccECP pseudo potential
#element_list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
element_list = ["H", "He", "Li", "C"]
ccECP = ccECP(pseudo_potential_output_dir=pseudopotentials_test_dir)
ccECP.to_file(basis_list=["cc-pVDZ"], element_list=element_list)

## pseudo potential classes
os.chdir(pseudopotentials_test_dir)
file_name="C_ccECP.pseudo"
pseudopotentials=Pseudopotentials.parse_pseudopotential_from_gamess_format_files(files=[file_name, file_name])
pseudopotentials.set_cutoffs()
pseudopotentials.write_pseudopotential_turborvb_file()