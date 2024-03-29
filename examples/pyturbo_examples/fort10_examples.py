#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys, shutil

#set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo')
logger.setLevel("INFO")
stream_handler = StreamHandler()
stream_handler.setLevel("INFO")
handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)

#pyturbo modules
from turbogenius.pyturbo.io_fort10 import IO_fort10
from turbogenius.utils_workflows.env import turbo_genius_root

# fort10
prefix="fort10"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

fort10 = IO_fort10("fort.10_hydrogen")
#fort10 = IO_fort10("fort.10_benzene")
#fort10 = IO_fort10("fort.10_SiO2")

fort10.f10header.io_flag = -1 # => fort.10 is automatically overwritten!!

logger.debug(fort10.f10header.io_flag)
fort10.f10header.io_flag = -1 * fort10.f10header.io_flag # => fort.10 is automatically overwritten!!
logger.debug(fort10.f10header.io_flag)

#logger.debug(fort10.f10structure.vec_a)
#logger.debug(fort10.f10structure.vec_b)
#logger.debug(fort10.f10structure.vec_c)
#logger.debug(fort10.f10structure.atomic_numbers)
#logger.debug(fort10.f10structure.positions)
#logger.debug(fort10.f10forceconstraint.atom_label)
#logger.debug(fort10.f10forceconstraint.direction)
#logger.debug(fort10.f10detbasissets.mo_coefficient)
#logger.debug(fort10.f10detbasissets.exponent)
#logger.debug(fort10.f10detbasissets.coefficient)
#logger.debug(fort10.f10jasbasissets.coefficient)
#logger.debug(fort10.f10detocc.occ)
#logger.debug(fort10.f10jasocc.occ)
logger.info(fort10.f10detmatrix.row)
logger.info(fort10.f10jasmatrix.row)
logger.info(fort10.f10jasmatrix.col)
logger.info(fort10.f10jasmatrix.coeff)
coeff=fort10.f10jasmatrix.coeff
coeff[3]=1.0
fort10.f10jasmatrix.coeff=coeff
#logger.debug(fort10.f10jasmatrix.row)
#logger.debug(fort10.f10detmat_sym.row)
#logger.debug(fort10.f10jasmat_sym.row)
#logger.debug(fort10.f10detbasis_sym.basis_index)
#logger.debug(fort10.f10jasbasis_sym.basis_index)
#logger.debug(fort10.f10jastwobody.twobody_list)
#logger.debug(fort10.f10jastwobody.onebody_list)
