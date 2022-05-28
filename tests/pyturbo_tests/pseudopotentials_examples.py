#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys

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
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from turbogenius.pyturbo.pseudopotentials import Pseudopotentials
from turbogenius.pyturbo.utils.downloader import BFD, ccECP

#basis sets

## basis set downloads
pseudopotentials_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "pseudopotentials")

os.chdir(pseudopotentials_test_dir)
# BFD pseudo potential
bfd=BFD(pseudo_potential_output_dir=pseudopotentials_test_dir)
bfd.to_file(element_list=["C"], basis_list=["vtz"])

# ccECP pseudo potential
element_list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]
ccECP = ccECP(pseudo_potential_output_dir=pseudopotentials_test_dir)
ccECP.to_file(basis_list=["cc-pVDZ"], element_list=element_list)

## pseudo potential classes
os.chdir(pseudopotentials_test_dir)
file_name="C_ccECP.pseudo"
pseudopotentials=Pseudopotentials.parse_pseudopotential_from_gamess_format_files(files=[file_name, file_name])
pseudopotentials.set_cutoffs()
pseudopotentials.write_pseudopotential_turborvb_file()