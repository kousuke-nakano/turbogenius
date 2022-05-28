#!python -u
# -*- coding: utf-8 -*-

from __future__ import print_function

# python modules
import os, sys
import subprocess

# set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo').getChild(__name__)

# pyturbo module
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# turbo-genius related path lists
pyturbo_root=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../'))
pyturbo_source_dir=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../'))
pyturbo_data_dir=os.path.join(pyturbo_source_dir, 'data')
pyturbo_tmp_dir=os.path.join(os.path.abspath(os.environ['HOME']), '.pyturbo_tmp')

# generate pyturbo temp. dir.
os.makedirs(pyturbo_tmp_dir, exist_ok=True)

# turborvb related binary paths
if "TURBORVB_ROOT" in os.environ:
    turborvb_root = os.environ["TURBORVB_ROOT"]
    turborvb_bin_root = os.path.join(turborvb_root, "bin")
else:
    try:
        sys_env = os.environ.copy()
        cmd = "which readalles.x"
        turborvb_root = os.path.dirname(os.path.dirname(subprocess.check_output(cmd, shell=True, env=sys_env))).decode()
        turborvb_bin_root = os.path.join(turborvb_root, "bin")
    except:
        raise ValueError("Set TURBORVB_ROOT (e.g., export TURBORVB_ROOT=XXX in ~.bashrc)")
if "TURBOMAKEFORT10_RUN_COMMAND" in os.environ:
    turbo_makefort10_run_command = os.environ["TURBOMAKEFORT10_RUN_COMMAND"]
else:
    turbo_makefort10_run_command = "makefort10.x"
if "TURBOCONVERTFORT10_RUN_COMMAND" in os.environ:
    turbo_convertfort10_run_command = os.environ["TURBOCONVERTFORT10_RUN_COMMAND"]
else:
    turbo_convertfort10_run_command = "convertfort10.x"
if "TURBOCONVERTFORT10MOL_RUN_COMMAND" in os.environ:
    turbo_convertfort10mol_run_command = os.environ["TURBOCONVERTFORT10MOL_RUN_COMMAND"]
else:
    turbo_convertfort10mol_run_command = "convertfort10mol.x"
if "TURBOPREP_RUN_COMMAND" in os.environ:
    turbo_prep_run_command = os.environ["TURBOPREP_RUN_COMMAND"]
else:
    turbo_prep_run_command = "prep-serial.x"
if "TURBOREADFORWARD_RUN_COMMAND" in os.environ:
    turbo_readforward_run_command = os.environ["TURBOREADFORWARD_RUN_COMMAND"]
else:
    turbo_readforward_run_command = "readforward-serial.x"
if "TURBOVMC_RUN_COMMAND" in os.environ:
    turbo_qmc_run_command = os.environ["TURBOVMC_RUN_COMMAND"]
else:
    turbo_qmc_run_command = "turborvb-serial.x"
if "TURBOFORCEVMC_RUN_COMMAND" in os.environ:
    turbo_forcevmc_run_command = os.environ["TURBOFORCEVMC_RUN_COMMAND"]
else:
    turbo_forcevmc_run_command = "forcevmc.sh"
if "TURBOFORCEVMC_KPOINTS_RUN_COMMAND" in os.environ:
    turbo_forcevmc_kpoints_run_command = os.environ["TURBOFORCEVMC_KPOINTS_RUN_COMMAND"]
else:
    turbo_forcevmc_kpoints_run_command = "forcevmc_kpoints.sh"
if "TURBOFORCEFN_RUN_COMMAND" in os.environ:
    turbo_forcefn_run_command = os.environ["TURBOFORCEFN_RUN_COMMAND"]
else:
    turbo_forcefn_run_command = "forcefn.sh"
if "TURBOFORCEFN_KPOINTS_RUN_COMMAND" in os.environ:
    turbo_forcefn_kpoints_run_command = os.environ["TURBOFORCEFN_KPOINTS_RUN_COMMAND"]
else:
    turbo_forcefn_kpoints_run_command = "forcefn_kpoints.sh"
if "TURBOCOPYJAS_RUN_COMMAND" in os.environ:
    turbo_copyjas_command = os.environ["TURBOCOPYJAS_RUN_COMMAND"]
else:
    turbo_copyjas_command = "copyjas.x"
if "TURBOCONVERTPFAFF_RUN_COMMAND" in os.environ:
    turbo_convertfortpfaff_run_command = os.environ["TURBOCONVERTPFAFF_RUN_COMMAND"]
else:
    turbo_convertfortpfaff_run_command = "convertpfaff.x"
