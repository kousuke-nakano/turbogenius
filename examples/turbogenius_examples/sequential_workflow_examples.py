#!/usr/bin/env python
# coding: utf-8

# TREXIO -> JDFT WF (fort.10) -> VMCopt (only jastrow) -> VMC (JDFT)

#Logger
from logging import config, getLogger, StreamHandler, Formatter
#log level (CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET)
turbo_genius_log_level="INFO"
pyturbo_log_level="INFO"

logger = getLogger('Turbo-Genius')
logger.setLevel(turbo_genius_log_level)
stream_handler = StreamHandler()
stream_handler.setLevel(turbo_genius_log_level)
handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)

logger_p = getLogger('pyturbo')
logger_p.setLevel(pyturbo_log_level)
stream_handler_p = StreamHandler()
stream_handler_p.setLevel(pyturbo_log_level)
handler_format_p = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler_p.setFormatter(handler_format_p)
logger_p.addHandler(stream_handler_p)

# python packages
import numpy as np
import os, sys
import shutil
import pickle

# turbo-genius packages
from turbogenius.trexio_to_turborvb import trexio_to_turborvb_wf
from turbogenius.vmc_opt_genius import VMCopt_genius
from turbogenius.vmc_genius import VMC_genius

# pyturbo package
from turbogenius.pyturbo.basis_set import Jas_Basis_sets

# test calc. dir
from turbogenius.utils_workflows.env import turbo_genius_root

# convertfort10mol
prefix="sequential_workflow"
example_root_dir=os.path.join(turbo_genius_root, "examples", "turbogenius_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

# trexio filename
trexio_filename = "H2_trexio.hf"

# jastrow optimization (vmcopt)
vmcopt_optsteps = 100
vmcopt_steps = 50
vmcopt_bin_block = 1
vmcopt_warmupblocks = 0
vmcopt_num_walkers = 40
vmcopt_optimizer = "lr"
vmcopt_learning_rate = 0.35
vmcopt_regularization = 0.001
vmcopt_warmupoptsteps = 80

# vmc
vmc_steps= 300
vmc_bin_block = 10
vmc_warmupblocks = 5
vmc_num_walkers = 40

###############################################
# Start a workflow
###############################################
root_dir=os.getcwd()

#******************
#! TREXIO -> TurboRVB WF
#******************
trexio_dir=os.path.join(root_dir, "01trexio")
os.makedirs(trexio_dir, exist_ok=True)
shutil.copy(os.path.join(root_dir,trexio_filename), os.path.join(trexio_dir,trexio_filename))
os.chdir(trexio_dir)

H_jastrow_basis="""
        S  1
        1       1.873529  1.00000000
        S  1
        1       0.343709  1.00000000
        S  1
        1       0.139013  1.00000000
        P  1
        1       0.740212  1.00000000
"""

H2_jas_basis_sets=Jas_Basis_sets.parse_basis_sets_from_texts([H_jastrow_basis, H_jastrow_basis], format="gamess")
trexio_to_turborvb_wf(trexio_file=os.path.join(trexio_dir, trexio_filename), jas_basis_sets=H2_jas_basis_sets)

os.chdir(root_dir)

#******************
#! Jastrow optimization
#******************

vmcopt_dir=os.path.join(root_dir, "02vmcopt")
os.makedirs(vmcopt_dir, exist_ok=True)
copy_files=["fort.10", "pseudo.dat"]
for file in copy_files:
    shutil.copy(os.path.join(trexio_dir,file), os.path.join(vmcopt_dir,file))
os.chdir(vmcopt_dir)

vmcopt_genius = VMCopt_genius(
    vmcoptsteps=vmcopt_optsteps,
    steps=vmcopt_steps,
    bin_block=vmcopt_bin_block,
    warmupblocks=vmcopt_warmupblocks,
    num_walkers=vmc_num_walkers,
    optimizer=vmcopt_optimizer,
    learning_rate=vmcopt_learning_rate,
    regularization=vmcopt_regularization,
    opt_onebody=True,
    opt_twobody=True,
    opt_det_mat=False,
    opt_jas_mat=True,
    opt_det_basis_exp=False,
    opt_jas_basis_exp=False,
    opt_det_basis_coeff=False,
    opt_jas_basis_coeff=False,
    twist_average=False,
)

vmcopt_genius.generate_input()
vmcopt_genius.run()
vmcopt_genius.average(optwarmupsteps=vmcopt_warmupoptsteps, graph_plot=True)

os.chdir(root_dir)

#******************
#! VMC
#******************
vmc_dir=os.path.join(root_dir, "03vmc")
os.makedirs(vmc_dir, exist_ok=True)

copy_files=["fort.10", "pseudo.dat"]
for file in copy_files:
    shutil.copy(os.path.join(vmcopt_dir,file), os.path.join(vmc_dir, file))
os.chdir(vmc_dir)

vmc_genius=VMC_genius(
                vmcsteps=vmc_steps,
                num_walkers=vmc_num_walkers,
                )

vmc_genius.generate_input()
vmc_genius.run()
vmc_genius.compute_energy_and_forces(bin_block=vmc_bin_block, warmupblocks=vmc_warmupblocks)

energy, error= vmc_genius.energy, vmc_genius.energy_error
logger.info(f"VMC-JDFT energy = {energy:.5f} +- {error:3f} Ha")

os.chdir(root_dir)
