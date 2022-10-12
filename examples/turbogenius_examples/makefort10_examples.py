#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys
import shutil

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
from turbogenius.makefort10_genius import Makefort10_genius
from turbogenius.utils_workflows.env import turbo_genius_root

# makefort10
prefix="makefort10"
example_root_dir=os.path.join(turbo_genius_root, "examples", "turbogenius_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))
makefort10_test_dir=os.path.join(example_root_dir, prefix)

#test ok
# Hydrogen dimer, ccECP
structure_file = "H2.xyz"
makefort10_genius = Makefort10_genius(
                                      structure_file=structure_file,
                                      det_basis_set="cc-pVTZ",
                                      jas_basis_set="cc-pVTZ",
                                      det_contracted_flag=True,
                                      jas_contracted_flag=True,
                                      all_electron_jas_basis_set=True,
                                      pseudo_potential="ccECP",
                                      det_cut_basis_option=True,
                                      jas_cut_basis_option=True,
                                      complex=False
                                      )
makefort10_genius.generate_input()
makefort10_genius.run()
makefort10_genius.check_results()
shutil.move(os.path.join(makefort10_test_dir,"fort.10_new"), os.path.join(makefort10_test_dir,"fort.10_H2"))
shutil.move(os.path.join(makefort10_test_dir,"makefort10.input"), os.path.join(makefort10_test_dir,"makefort10_H2.input"))
shutil.move(os.path.join(makefort10_test_dir,"structure.xsf"), os.path.join(makefort10_test_dir,"H2_gen_structure.xsf"))

# Benzene, all-electrons
structure_file = "benzene.xyz"
makefort10_genius = Makefort10_genius(
                                      structure_file=structure_file,
                                      det_basis_set="cc-pVQZ",
                                      jas_basis_set="cc-pVDZ",
                                      det_contracted_flag=True,
                                      jas_contracted_flag=True,
                                      all_electron_jas_basis_set=True,
                                      pseudo_potential=None,
                                      det_cut_basis_option=True,
                                      jas_cut_basis_option=True,
                                      complex=False
                                      )
makefort10_genius.generate_input()
makefort10_genius.run()
makefort10_genius.check_results()
shutil.move(os.path.join(makefort10_test_dir,"makefort10.input"), os.path.join(makefort10_test_dir,"makefort10_benzene.input"))
shutil.move(os.path.join(makefort10_test_dir,"fort.10_new"), os.path.join(makefort10_test_dir,"fort.10_benzene"))
shutil.move(os.path.join(makefort10_test_dir,"structure.xsf"), os.path.join(makefort10_test_dir,"benzene_gen_structure.xsf"))

# Diamond, PP, orthorhombic case, BFD
structure_file = "Diamond.cif"
makefort10_genius = Makefort10_genius(
                                      structure_file=structure_file,
                                      det_basis_set="cc-pVQZ",
                                      jas_basis_set="cc-pVDZ",
                                      det_contracted_flag=True,
                                      jas_contracted_flag=True,
                                      all_electron_jas_basis_set=True,
                                      pseudo_potential=None,
                                      det_cut_basis_option=True,
                                      jas_cut_basis_option=True,
                                      complex=False
                                      )
makefort10_genius.generate_input()
makefort10_genius.run()
makefort10_genius.check_results()
shutil.move(os.path.join(makefort10_test_dir,"makefort10.input"), os.path.join(makefort10_test_dir,"makefort10_diamond.input"))
shutil.move(os.path.join(makefort10_test_dir,"fort.10_new"), os.path.join(makefort10_test_dir,"fort.10_diamond"))
shutil.move(os.path.join(makefort10_test_dir,"structure.xsf"), os.path.join(makefort10_test_dir,"diamond_gen_structure.xsf"))