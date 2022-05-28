#!/usr/bin/env python
# coding: utf-8

# pyscf -> TREXIO -> JDFT WF(fort.10) -> VMCopt(only jastrow) -> VMC(JDFT) -> LRDMC(JDFT, a->0) a simple workflow (H2 dimer)

#Logger
from logging import config, getLogger, StreamHandler, Formatter
#log level (CRITICAL, ERROR, WARNING, INFO, DEBUG, NOTSET)
turbo_genius_log_level="INFO"
pyturbo_log_level="INFO"

logger = getLogger('Turbo-Genius')
logger.setLevel(turbo_genius_log_level)
stream_handler = StreamHandler()
stream_handler.setLevel(turbo_genius_log_level)
#handler_format = Formatter('%(message)s')
handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)

logger_p = getLogger('pyturbo')
logger_p.setLevel(pyturbo_log_level)
stream_handler_p = StreamHandler()
stream_handler_p.setLevel(pyturbo_log_level)
#handler_format = Formatter('%(message)s')
handler_format_p = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler_p.setFormatter(handler_format_p)
logger_p.addHandler(stream_handler_p)

# python packages
import numpy as np
import os, sys
import shutil
import pickle
import numpy as np

# turbo-genius packages
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from turbogenius.pyscf_trexio.pyscf_wrapper import Pyscf_wrapper
from turbogenius.trexio_to_turborvb import trexio_to_turborvb_wf
from turbogenius.trexio_wrapper import Trexio_wrapper_r
from turbogenius.pyscf_trexio.pyscf_to_trexio import pyscf_to_trexio
from turbogenius.vmc_opt_genius import VMCopt_genius
from turbogenius.vmc_genius import VMC_genius
from turbogenius.lrdmc_genius import LRDMC_genius, LRDMC_ext_genius

# pyturbo package
from turbogenius.pyturbo.basis_set import Det_Basis_sets, Jas_Basis_sets

# test calc. dir
sequential_workflow_test_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequential_workflow")
os.chdir(sequential_workflow_test_dir)

#**************************************
# Define all input parameters here!!
#**************************************
## pyscf
pyscf_rerun=False
pyscf_pkl="pyscf_genius_dummy.pkl"
structure_file="H2.xyz"
charge = 0
spin = 0
basis = "ccecp-ccpvtz"
ecp = 'ccecp'
scf_method = "DFT"  # HF or DFT
dft_xc = "LDA_X,LDA_C_PZ"
pyscf_output = "out_H2_pyscf"
trexio_filename = "H2_trexio.hdf5"

## jastrow optimization (vmc)
jasopt_rerun=False
jasopt_pkl="vmcopt_genius.pkl"
jasopt_target_error_bar=5.0e-3 # Ha
jasopt_trial_optsteps = 30
jasopt_trial_steps = 10
jasopt_production_optsteps = 100
jasopt_bin_block = 1
jasopt_warmupblocks = 0
jasopt_num_walkers = 40 # default -1 -> num of MPI process.
jasopt_optimizer = "lr"
jasopt_learning_rate = 0.35
jasopt_regularization = 0.001

## vmc
vmc_jdft_rerun=False
vmc_jdft_pkl="vmc_genius_jdft.pkl"
vmc_target_error_bar=1.0e-3 # Ha
vmc_trial_steps= 300
vmc_bin_block = 10
vmc_warmupblocks = 5
vmc_num_walkers = 40 # default -1 -> num of MPI process.
vmc_twist_average=False
vmc_force_calc_flag=False

## lrdmc
lrdmc_jdft_rerun=True
lrdmc_jdft_pkl="lrdmcext_genius_jdft.pkl"
lrdmc_target_error_bar=1.0e-3 # Ha
lrdmc_trial_steps= 300
lrdmc_bin_block = 10
lrdmc_warmupblocks = 10
lrdmc_correcting_factor = 10
lrdmc_num_walkers = 40 # default -1 -> num of MPI process.
lrdmc_alat_list = [-0.20, -0.30, -0.40]
lrdmc_polynominal_order = 2
lrdmc_twist_average = False
lrdmc_force_calc_flag = False
lrdmc_nonlocalmoves = "dlatm"  # tmove, dla, dlatm

###############################################
# Start a workflow
###############################################
root_dir=os.getcwd()
logger.info(f"Project root dir = {root_dir}")
pkl_dir=os.path.join(root_dir, "pkl")
os.makedirs(pkl_dir, exist_ok=True)

#******************
#! pyscf
#******************
pyscf_dir=os.path.join(root_dir, "pyscf")
trexio_dir=os.path.join(root_dir, "trexio_to_turbo")

if pyscf_rerun or not os.path.isfile(os.path.join(pkl_dir, pyscf_pkl)):
    logger.info(f"Start: pyscf calc.")
    os.makedirs(pyscf_dir, exist_ok=True)
    os.chdir(pyscf_dir)

    pyscf_wrapper=Pyscf_wrapper(structure_file=os.path.join(root_dir, structure_file))
    pyscf_wrapper.run_pyscf(
        charge=charge,
        spin=spin,
        basis=basis,
        ecp=ecp,
        scf_method=scf_method,
        dft_xc=dft_xc,
        pyscf_output=pyscf_output
    )
    pyscf_to_trexio(trexio_filename=os.path.join(pyscf_dir, trexio_filename))

    logger.info(f"End: pyscf calc.")
    os.chdir(root_dir)

    #! TREXIO -> TurboRVB conversion
    logger.info(f"Start: trexio -> turborvb conversion")
    os.makedirs(trexio_dir, exist_ok=True)
    os.chdir(trexio_dir)

    ## set jastrow basis (gammes format)
    jastrow_basis_dict={
        'H':
            """
            S  1
              1       1.873529  1.00000000
            S  1
              1       0.343709  1.00000000
            S  1
              1       0.139013  1.00000000
            P  1
              1       0.740212  1.00000000
            """
    }

    trexio_r=Trexio_wrapper_r(trexio_file=os.path.join(pyscf_dir, trexio_filename))
    jastrow_basis_list=[jastrow_basis_dict[element] for element in trexio_r.labels_r]
    jas_basis_sets=Jas_Basis_sets.parse_basis_sets_from_texts(jastrow_basis_list, format="gamess")

    ## trexio -> turborvb_wf
    trexio_to_turborvb_wf(trexio_file=os.path.join(pyscf_dir, trexio_filename), jas_basis_sets=jas_basis_sets, max_occ_conv=1.0e-3)

    logger.info(f"End: trexio -> turborvb conversion")

    # unfortunately, pyscf cannot be pickled due to the following error:
    # cannot serialize '_io.TextIOWrapper' object
    # so, here a dummy pickled file is dumped.
    with open(os.path.join(pkl_dir, pyscf_pkl), "wb") as f:
        pickle.dump("Dummy", f)

    os.chdir(root_dir)

#******************
#! Jastrow optimization
#******************

jasopt_dir=os.path.join(root_dir, "jasopt")
if jasopt_rerun or not os.path.isfile(os.path.join(pkl_dir, jasopt_pkl)):
    logger.info(f"Start: Jastrow optimization")
    os.makedirs(jasopt_dir, exist_ok=True)
    os.chdir(jasopt_dir)

    jasopt_test_dir=os.path.join(jasopt_dir, "jasopt_test")
    os.makedirs(jasopt_test_dir, exist_ok=True)
    os.chdir(jasopt_test_dir)
    copy_files=["fort.10", "pseudo.dat"]
    for file in copy_files:
        shutil.copy(os.path.join(trexio_dir,file), os.path.join(jasopt_test_dir,file))

    ## estimated necesary steps per optimization to achieve the target error bar.
    jasopt_optwarmupsteps=0
    vmcopt_genius = VMCopt_genius(
        vmcoptsteps=jasopt_trial_optsteps,
        optwarmupsteps=jasopt_optwarmupsteps,
        steps=jasopt_trial_steps,
        bin_block=jasopt_bin_block,
        warmupblocks=jasopt_warmupblocks,
        num_walkers=jasopt_num_walkers,
        optimizer=jasopt_optimizer,
        learning_rate=jasopt_learning_rate,
        regularization=jasopt_regularization,
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
    energy, error=vmcopt_genius.get_energy()
    logger.info(f"The vmc energy at the final step is {energy[-1]:.5f} Ha")
    logger.info(f"The error bar of the vmc energy at the final step is {error[-1]:.5f} Ha per mcmc step={(jasopt_trial_steps- jasopt_optwarmupsteps*jasopt_bin_block)}")
    jasopt_steps_estimated_proper=int((jasopt_trial_steps- jasopt_optwarmupsteps*jasopt_bin_block) * (error[-1] / jasopt_target_error_bar) ** 2)
    logger.info(f"The target error bar per optstep is {jasopt_target_error_bar:.5f} Ha")
    logger.info(f"The estimated steps to achieve the target error bar is {jasopt_steps_estimated_proper:d} steps")
    estimated_time_for_1_generation=vmcopt_genius.get_estimated_time_for_1_generation()

    os.chdir(jasopt_dir)

    #! VMCopt (production)
    estimated_time= estimated_time_for_1_generation * jasopt_steps_estimated_proper * jasopt_production_optsteps
    jasopt_optwarmupsteps_default=int(0.80 * jasopt_production_optsteps)
    logger.info(f"Production run starts.")
    logger.info(f"Jastrow optimization, production run step={jasopt_production_optsteps}")
    logger.info(f"optwarmupsteps is set to {jasopt_optwarmupsteps_default} (the first 80% steps)")
    logger.info(f"The final 20% steps will be used for averaging parameters.")
    logger.info(f"Estimated time = {estimated_time:.0f} sec.")

    copy_files=["fort.10", "pseudo.dat"]
    for file in copy_files:
        shutil.copy(os.path.join(jasopt_test_dir,file), os.path.join(jasopt_dir,file))

    ## estimated necesary steps per optimization to achieve the target error bar.
    vmcopt_genius = VMCopt_genius(
        vmcoptsteps=jasopt_production_optsteps,
        optwarmupsteps=jasopt_optwarmupsteps_default,
        steps=jasopt_steps_estimated_proper,
        bin_block=jasopt_bin_block,
        warmupblocks=jasopt_warmupblocks,
        num_walkers=vmc_num_walkers,
        optimizer=jasopt_optimizer,
        learning_rate=jasopt_learning_rate,
        regularization=jasopt_regularization,
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
    vmcopt_genius.run_all()

    with open(os.path.join(pkl_dir, jasopt_pkl), "wb") as f:
        pickle.dump(vmcopt_genius, f)

    os.chdir(root_dir)

else:
    with open(os.path.join(pkl_dir, jasopt_pkl), "rb") as f:
        jasopt_pkl=pickle.load(f)

#******************
#! VMC
#******************
vmc_jdft_dir=os.path.join(root_dir, "vmc_jdft")
if vmc_jdft_rerun or not os.path.isfile(os.path.join(pkl_dir, vmc_jdft_pkl)):
    logger.info(f"Start: VMC(JDFT) calculation")
    os.makedirs(vmc_jdft_dir, exist_ok=True)
    os.chdir(vmc_jdft_dir)

    copy_files=["fort.10", "pseudo.dat"]
    for file in copy_files:
        shutil.copy(os.path.join(jasopt_dir,file), os.path.join(vmc_jdft_dir, file))

    vmc_jdft_test_dir = os.path.join(vmc_jdft_dir, "vmc_jdft_test")
    os.makedirs(vmc_jdft_test_dir, exist_ok=True)
    os.chdir(vmc_jdft_test_dir)

    copy_files=["fort.10", "pseudo.dat"]
    for file in copy_files:
        shutil.copy(os.path.join(vmc_jdft_dir, file), os.path.join(vmc_jdft_test_dir, file))

    ## estimated necesary steps per optimization to achieve the target error bar.
    assert vmc_trial_steps > vmc_bin_block * vmc_warmupblocks
    vmc_genius_jdft=VMC_genius(
                     vmcsteps=vmc_trial_steps,
                     bin_block=vmc_bin_block,
                     warmupblocks=vmc_warmupblocks,
                     num_walkers=vmc_num_walkers,
                     twist_average=vmc_twist_average,
                     force_calc_flag=vmc_force_calc_flag
                     )

    vmc_genius_jdft.run_all()
    energy, error= vmc_genius_jdft.energy, vmc_genius_jdft.energy_error
    logger.info(f"The error bar of the vmc energy {error:.5f} Ha per mcmc step={(vmc_trial_steps-vmc_warmupblocks*vmc_bin_block)}")
    vmc_jdft_steps_estimated_proper=int((vmc_trial_steps-vmc_warmupblocks*vmc_bin_block) * (error / vmc_target_error_bar) ** 2)
    logger.info(f"The target error bar per optstep is {vmc_target_error_bar:.5f} Ha")
    logger.info(f"The estimated steps to achieve the target error bar is {vmc_jdft_steps_estimated_proper:d} steps")
    estimated_time_for_1_generation=vmc_genius_jdft.get_estimated_time_for_1_generation()
    estimated_time = estimated_time_for_1_generation * vmc_jdft_steps_estimated_proper

    os.chdir(vmc_jdft_dir)

    #! VMC (production)

    logger.info(f"Production run starts.")
    logger.info(f"Estimated time = {estimated_time:.0f} sec.")

    ## estimated necesary steps per optimization to achieve the target error bar.
    vmc_genius_jdft=VMC_genius(
                     vmcsteps=vmc_jdft_steps_estimated_proper,
                     bin_block=vmc_bin_block,
                     warmupblocks=vmc_warmupblocks,
                     num_walkers=vmc_num_walkers,
                     twist_average=vmc_twist_average,
                     force_calc_flag=vmc_force_calc_flag
                     )
    vmc_genius_jdft.run_all()
    energy, error= vmc_genius_jdft.energy, vmc_genius_jdft.energy_error
    logger.info(f"Final VMC-JDFT energy = {energy:.5f} +- {error:3f} Ha")

    with open(os.path.join(pkl_dir, vmc_jdft_pkl), "wb") as f:
        pickle.dump(vmc_genius_jdft, f)

    os.chdir(root_dir)

else:
    with open(os.path.join(pkl_dir, vmc_jdft_pkl), "rb") as f:
        vmc_genius_jdft=pickle.load(f)

#******************
#! LRDMC
#******************
lrdmc_jdft_dir=os.path.join(root_dir, "lrdmc_jdft")
if lrdmc_jdft_rerun or not os.path.isfile(os.path.join(pkl_dir, lrdmc_jdft_pkl)):
    logger.info(f"Start: LRDMC(JDFT) calculation")
    os.makedirs(lrdmc_jdft_dir, exist_ok=True)
    os.chdir(lrdmc_jdft_dir)
    copy_files=["fort.10", "pseudo.dat"]
    for file in copy_files:
        shutil.copy(os.path.join(vmc_jdft_dir,file), os.path.join(lrdmc_jdft_dir,file))

    lrdmc_jdft_test_dir = os.path.join(lrdmc_jdft_dir, "lrdmc_jdft_test")
    os.makedirs(lrdmc_jdft_test_dir, exist_ok=True)
    os.chdir(lrdmc_jdft_test_dir)

    copy_files=["fort.10", "pseudo.dat"]
    for file in copy_files:
        shutil.copy(os.path.join(lrdmc_jdft_dir, file), os.path.join(lrdmc_jdft_test_dir, file))

    ## estimated necesary steps per optimization to achieve the target error bar.
    assert lrdmc_trial_steps > lrdmc_bin_block * lrdmc_warmupblocks

    with open(os.path.join(pkl_dir, vmc_jdft_pkl), "rb") as f:
        vmc_genius_jdft=pickle.load(f)

    energy, error = vmc_genius_jdft.energy, vmc_genius_jdft.energy_error

    trial_etry=energy
    trial_alat=np.min(lrdmc_alat_list)

    lrdmc_genius_jdft=LRDMC_genius(
                     lrdmcsteps=lrdmc_trial_steps,
                     bin_block=lrdmc_bin_block,
                     correcting_factor=lrdmc_correcting_factor,
                     warmupblocks=lrdmc_warmupblocks,
                     num_walkers=lrdmc_num_walkers,
                     alat=trial_alat,
                     etry=trial_etry,
                     twist_average=lrdmc_twist_average,
                     force_calc_flag=lrdmc_force_calc_flag
                     )

    lrdmc_genius_jdft.run_all()
    energy, error= lrdmc_genius_jdft.energy, lrdmc_genius_jdft.energy_error
    logger.info(f"The error bar of the lrdmc energy (alat={trial_alat}) is {error:.5f} Ha per mcmc step={(lrdmc_trial_steps-lrdmc_warmupblocks*lrdmc_bin_block)}")
    logger.info(f"The target error bar per optstep is {lrdmc_target_error_bar:.5f} Ha")
    lrdmc_jdft_steps_estimated_proper = int((lrdmc_trial_steps-lrdmc_warmupblocks*lrdmc_bin_block) * (error / lrdmc_target_error_bar) ** 2)
    logger.info(f"The estimated steps to achieve the target error bar is {lrdmc_jdft_steps_estimated_proper:d} steps")
    estimated_time_for_1_generation=lrdmc_genius_jdft.get_estimated_time_for_1_generation()
    estimated_time = estimated_time_for_1_generation * lrdmc_jdft_steps_estimated_proper

    os.chdir(lrdmc_jdft_dir)

    #! LRDMC (production)
    logger.info(f"Production run starts.")
    logger.info(f"Estimated time for the minimum alat = {estimated_time * (np.min(lrdmc_alat_list)/np.max(lrdmc_alat_list))**2:.0f} sec.")

    ## estimated necesary steps per optimization to achieve the target error bar.
    lrdmc_ext_genius_jdft=LRDMC_ext_genius(
                     lrdmcsteps=lrdmc_jdft_steps_estimated_proper,
                     bin_block=lrdmc_bin_block,
                     correcting_factor=lrdmc_correcting_factor,
                     warmupblocks=lrdmc_warmupblocks,
                     num_walkers=lrdmc_num_walkers,
                     alat_list=lrdmc_alat_list,
                     etry=trial_etry,
                     twist_average=lrdmc_twist_average,
                     force_calc_flag=lrdmc_force_calc_flag
                     )
    lrdmc_ext_genius_jdft.run_all()
    energy, error= lrdmc_ext_genius_jdft.energy, lrdmc_ext_genius_jdft.energy_error
    logger.info(f"Final LRDMC-JDFT energy (a->0) = {energy:.5f} +- {error:3f} Ha")

    with open(os.path.join(pkl_dir, lrdmc_jdft_pkl), "wb") as f:
        pickle.dump(lrdmc_ext_genius_jdft, f)

    os.chdir(root_dir)

else:
    with open(os.path.join(pkl_dir, lrdmc_jdft_pkl), "rb") as f:
        lrdmc_ext_genius_jdft=pickle.load(f)

logger.info("LRDMC (WF=JDFT,a->0) workflow ends.")