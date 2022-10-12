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
from turbogenius.pyturbo.structure import Structure
from turbogenius.utils_workflows.env import turbo_genius_root

#Structures
prefix="structures"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

# H2.xyz
structure = Structure.parse_structure_from_file(file="H2.xyz")
# Diamond.vasp
structure = Structure.parse_structure_from_file(file="Diamond.vasp")
# Diamond.cif
structure = Structure.parse_structure_from_file(file="Diamond.cif")
# benzene.xyz
structure = Structure.parse_structure_from_file(file="benzene.xyz")
# silicon_oxide.cif
structure = Structure.parse_structure_from_file(file="silicon_oxide.cif")

# bulk Al (from ase instance)
from ase.spacegroup import crystal
a = 4.05
atom = crystal('Al', [(0, 0, 0)], spacegroup=225, cellpar=[a, a, a, 90, 90, 90])
structure = Structure.parse_structure_from_ase_atom(atom)
#view(structure.get_ase_atom())