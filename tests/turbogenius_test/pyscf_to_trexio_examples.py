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
from turbogenius.pyscf_trexio.pyscf_wrapper import Pyscf_wrapper
from turbogenius.pyscf_trexio.pyscf_to_trexio import pyscf_to_trexio

pyscf_to_trexio_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "pyscf_to_trexio")

os.chdir(pyscf_to_trexio_test_dir)

# N2
structure_file="N2.xyz"
pyscf=Pyscf_wrapper(
        structure_file=structure_file,
)
pyscf.run_pyscf(pyscf_output="N2_pyscf.output")
pyscf_to_trexio(trexio_filename="N2_trexio.hdf5")

# CO2
structure_file="CO2.xyz"
pyscf=Pyscf_wrapper(
        structure_file=structure_file
)
pyscf.run_pyscf(pyscf_output="CO2_pyscf.output")
pyscf_to_trexio(trexio_filename="CO2_trexio.hdf5")
