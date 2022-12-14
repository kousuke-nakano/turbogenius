#!/usr/bin/env python
# coding: utf-8

"""

converter: TREXIO to TurboRVB WF (fort.10)

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

# module information pyscf_to_trexio
__author__ = "Kosuke Nakano"
__copyright__ = "Copyright 2022, Kosuke Nakano"
__version__ = "1.0.0"
__maintainer__ = "Kosuke Nakano"
__email__ = "kousuke_1123@icloud.com"
__date__ = "April 27, 2022"

# # TREXIO->TurboRVB conversion using pyturbo

# import python modules
import os,sys
import shutil
import argparse
import numpy as np
import glob

# logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)

# import pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pyturbo.structure import Structure, Cell
from pyturbo.pseudopotentials import Pseudopotentials
from pyturbo.basis_set import Det_Basis_sets, Jas_Basis_sets
from pyturbo.io_fort10 import IO_fort10
from pyturbo.makefort10 import Makefort10
from pyturbo.convertfort10mol import Convertfort10mol
from pyturbo.utils.utility import turbo_prim_orb_type_num, turbo_cont_orb_type_num, return_orbchr
from pyturbo.utils.utility import return_atomic_number

#import turbo-genius modules
from utils_workflows.env import turbo_genius_root, turbo_genius_tmp_dir
from trexio_wrapper import Trexio_wrapper_r

def trexio_to_turborvb_wf(trexio_file:str,
                          jas_basis_sets=Jas_Basis_sets(),
                          max_occ_conv:int=0,
                          mo_num_conv:int=-1,
                          only_mol:bool=True,
                          cleanup:bool=True
                          )->None:
    """
        Convert trexio file to TurboRVB WF file (fort.10)

        Args:
            trexio_file (str): TREXIO file name
            jas_basis_sets (Jas_basis_sets): Jastrow basis sets added to the TREXIO WF.
            max_occ_conv (int): maximum occ used for the conv, not used with mo_num
            mo_num_conv (int): num mo used for the conv, not used with max occ
            only_mol (bool): if True, only moleculer orbitals option = True in convertfort10mol
            cleanup (bool): clean up temporary files

    """
    #os.environ["DYLD_LIBRARY_PATH"]=os.environ["LIBRARY_PATH"]
    
    # prefix and file names
    logger.info(f"Input TREXIO file = {trexio_file}")
    trexio_r=Trexio_wrapper_r(trexio_file=trexio_file)

    # read electron num
    num_ele_up = trexio_r.num_ele_up
    num_ele_dn = trexio_r.num_ele_dn
    num_ele_total = num_ele_up + num_ele_dn
    
    # read structure info.
    nucleus_num_r = trexio_r.nucleus_num_r
    labels_r = trexio_r.labels_r
    charges_r = trexio_r.charges_r
    coords_r = trexio_r.coords_r
    total_charge =  np.sum(charges_r) - num_ele_total
    
    atomic_number_list=[return_atomic_number(Z) for Z in labels_r]
    atomic_number_unique=list(set(atomic_number_list))
    element_list=labels_r
    
    # check data
    logger.debug(f"nucleus_num from TREXIO file ---> {nucleus_num_r}")
    logger.debug(f"nucleus_label from TREXIO file ---> {labels_r}")
    logger.debug(f"nucleus_charge from TREXIO file ---> {charges_r}")
    logger.debug(f"nucleus_coord from TREXIO file ----> {coords_r}")
    
    # structure
    # Check if the given TREXIO file contains a crystal structure or an open system.
    # here check.... to be implemented
    pbc_flag=trexio_r.periodic
    
    if pbc_flag: # Crystal structure
        cell_a=trexio_r.cell_a
        cell_b=trexio_r.cell_b
        cell_c = trexio_r.cell_c
        k1,k2,k3 = trexio_r.k_point
        #phase_up=[+k1, +k2, +k3]
        #phase_dn=[-k1, -k2, -k3]
        # Thank you very much Michele, this opposite spin approach is obviously not general
        # because psi(k) != psi(-k) except for TRIM unless the crystal structure has the inversion symmetry.
        # so, we have not choice but to choosing the same phase for up and dn.
        phase_up=[+k1, +k2, +k3]
        phase_dn=[+k1, +k2, +k3]

        cell=Cell(
                  vec_a=cell_a,
                  vec_b=cell_b,
                  vec_c=cell_c,
                  )
        structure=Structure(
                        cell=cell,
                        atomic_numbers=atomic_number_list,
                        element_symbols=labels_r,
                        positions=coords_r
                        )
        complex_flag=trexio_r.complex_flag

    else: # open system
        structure=Structure(
                        atomic_numbers=atomic_number_list,
                        element_symbols=labels_r,
                        positions=coords_r
                        )
        phase_up=[0.0, 0.0, 0.0]
        phase_dn=[0.0, 0.0, 0.0]
        complex_flag = trexio_r.complex_flag

    logger.info(complex_flag)
    # view
    #from ase.visualize import view
    #atom=structure.get_ase_atom()
    #view(atom)
    
    # Reading basis sets info
    basis_type=trexio_r.basis_type
    basis_shell_num=trexio_r.basis_shell_num
    basis_shell_index=trexio_r.basis_shell_index
    basis_prim_num=trexio_r.basis_prim_num
    basis_nucleus_index=trexio_r.basis_nucleus_index
    basis_shell_ang_mom=trexio_r.basis_shell_ang_mom
    basis_shell_factor=trexio_r.basis_shell_factor
    basis_shell_index=trexio_r.basis_shell_index
    basis_exponent=trexio_r.basis_exponent
    basis_coefficient=trexio_r.basis_coefficient
    basis_prim_factor=trexio_r.basis_prim_factor
    
    logger.debug(len(basis_nucleus_index))
    logger.debug(len(basis_shell_index))
    logger.debug(basis_nucleus_index)
    logger.debug(basis_shell_index)
    
    # Pseudo potentials info
    ecp_max_ang_mom_plus_1=trexio_r.ecp_max_ang_mom_plus_1
    ecp_z_core=trexio_r.ecp_z_core
    ecp_num=trexio_r.ecp_num
    ecp_ang_mom=trexio_r.ecp_ang_mom
    ecp_nucleus_index=trexio_r.ecp_nucleus_index
    ecp_exponent=trexio_r.ecp_exponent
    ecp_coefficient=trexio_r.ecp_coefficient
    ecp_power=trexio_r.ecp_power
    
    # ao info
    ao_cartesian=trexio_r.ao_cartesian
    ao_num=trexio_r.ao_num
    ao_shell=trexio_r.ao_shell
    ao_normalization=trexio_r.ao_normalization
    
    logger.debug(len(ao_shell))
    logger.debug(ao_shell)
    
    # mo info
    mo_type=trexio_r.mo_type
    mo_num=trexio_r.mo_num
    mo_coefficient=trexio_r.mo_coefficient
    mo_occupation=trexio_r.mo_occupation
    if complex_flag:
        mo_coefficient_imag = trexio_r.mo_coefficient_imag
        mo_coefficient = [[complex(i,j) for i,j in zip(mo_real, mo_imag)] for mo_real, mo_imag in zip(mo_coefficient, mo_coefficient_imag)]
    
    logger.debug(ecp_max_ang_mom_plus_1)
    logger.debug(ecp_z_core)
    logger.debug(ecp_num)
    logger.debug(ecp_ang_mom)
    logger.debug(ecp_nucleus_index)
    logger.debug(ecp_exponent)
    logger.debug(ecp_coefficient)
    logger.debug(ecp_power)
    
    # basis sets
    shell_ang_mom_turbo_notation = []
    
    for i, ang_mom in enumerate(basis_shell_ang_mom):
        if list(basis_shell_index).count(i) > 1: # contracted orbitals
            shell_ang_mom_turbo_notation.append(turbo_cont_orb_type_num(return_orbchr(ang_mom)))
        else: # uncontracted orbitals
            shell_ang_mom_turbo_notation.append(turbo_prim_orb_type_num(return_orbchr(ang_mom)))
            
    # if complex_flag is true, we need a dummy coeff_imag!!
    if complex_flag:
        basis_coefficient_imag=[0.0] * len(basis_coefficient)
    else:
        basis_coefficient_imag=[]
    det_basis_sets = Det_Basis_sets(
                    nucleus_index=basis_nucleus_index,
                    shell_ang_mom=basis_shell_ang_mom,
                    shell_ang_mom_turbo_notation=shell_ang_mom_turbo_notation,
                    shell_factor=basis_shell_factor,
                    shell_index=basis_shell_index,
                    exponent=basis_exponent,
                    coefficient=basis_coefficient,
                    coefficient_imag=basis_coefficient_imag,
                    prim_factor=basis_prim_factor,
                    )
    
    # Pseudopotentials
    cutoff = [0.0] * len(ecp_z_core)
    pseudopotentials=Pseudopotentials(
                    max_ang_mom_plus_1=ecp_max_ang_mom_plus_1,
                    z_core=ecp_z_core,
                    cutoff=cutoff,
                    nucleus_index=ecp_nucleus_index,
                    element_list=element_list,
                    ang_mom=ecp_ang_mom,
                    exponent=ecp_exponent,
                    coefficient=ecp_coefficient,
                    power=ecp_power
                    )
    pseudopotentials.set_cutoffs()
    logger.debug(pseudopotentials.cutoff)
    pseudopotentials.write_pseudopotential_turborvb_file()

    # set jastrow_type:
    if jas_basis_sets.shell_num==0:
        jastrow_type = 0
    else:
        jastrow_type = -6 # the standard choice for PP calc.

    # makefort10
    namelist = Makefort10.read_default_namelist(structure=structure, 
                                                jastrow_type=jastrow_type,
                                                neldiff = num_ele_up - num_ele_dn,
                                                nel = num_ele_up + num_ele_dn)

    # check input related to phases
    if not len(phase_up) == 3:
        logger.error(f"len(phase_up) should be 3.")
        raise ValueError
    if not len(phase_dn) == 3:
        logger.error(f"len(phase_up) should be 3.")
        raise ValueError

    # check if opposite or same
    if np.array(phase_up) == np.array(phase_dn):
        logger.info(f"phase up == phase dn or at gamma points")
        logger.warning(f"forcesymm is true even for gamma points. The option to deactivate this will be implemented in the future.")
        logger.warning(f"forcesymm option will be activated.")
        namelist.set_parameter(parameter="forcesymm", value='.true.', namelist="&symmetries")

    elif np.array(phase_up) == -1 * np.array(phase_dn):
        logger.info(f"phase up == -1 * phase dn")
        logger.warning(f"forcesymm option will be deactivated.")
        namelist.set_parameter(parameter="forcesymm", value='.false.', namelist="&symmetries")

    else:
        logger.error(f"phase up is not equal to +1 * phase up and -1 * phase dn")
        raise ValueError

    # phase-up and phase-dn
    namelist.set_parameter(parameter="phase(1)", value=phase_up[0], namelist="&system")
    namelist.set_parameter(parameter="phase(2)", value=phase_up[1], namelist="&system")
    namelist.set_parameter(parameter="phase(3)", value=phase_up[2], namelist="&system")
    namelist.set_parameter(parameter="phasedo(1)", value=phase_dn[0], namelist="&system")
    namelist.set_parameter(parameter="phasedo(2)", value=phase_dn[1], namelist="&system")
    namelist.set_parameter(parameter="phasedo(3)", value=phase_dn[2], namelist="&system")

    if complex_flag:
        namelist.set_parameter(parameter="complexfort10", value='.true.', namelist="&system")
        #logger.error("complex is not implemented in the converter.")
    else:
        namelist.set_parameter(parameter="complexfort10", value='.false.', namelist="&system")

    makefort10 = Makefort10(
        structure=structure,
        det_basis_sets=det_basis_sets,
        jas_basis_sets=jas_basis_sets,
        pseudo_potentials=pseudopotentials,
        namelist=namelist
    )
    makefort10.generate_input('makefort10.input')
    makefort10.run()
    
    # rename fort.10
    shutil.move("fort.10_new", "fort.10_in")

    #here, specify how many MOs are used for the conversion!!
    if max_occ_conv == 0 and mo_num_conv == -1: # default
        logger.info(f"All the MOs ={mo_num} are used for the conversion.")
        mo_num_use=mo_num
        mo_coefficient=mo_coefficient
        mo_occ=mo_occupation
    elif max_occ_conv != 0 and mo_num_conv != -1:
        logger.error("max_occ_conv and mo_num_conv options cannot be used at the same time.")
        raise ValueError
    else:
        if max_occ_conv != 0:
            mo_used_index = []
            logger.info(f"mo_occ_thr={max_occ_conv}")
            logger.info(f"mo_occupation")
            logger.info(mo_occupation)
            for i in range(len(mo_occupation)):
                mo_used_index.append(i)
                if i==len(mo_occupation)-1:
                    break
                if mo_occupation[i] >= max_occ_conv and mo_occupation[i + 1] < max_occ_conv:
                    logger.info(f"mo_occ < {max_occ_conv} is 1-{i + 1}")
                    break
            mo_num_use = len(mo_used_index)
            mo_coefficient = [mo_coefficient[ppp] for ppp in mo_used_index]
            mo_occ = [mo_occupation[ppp] for ppp in mo_used_index]

        else: # mo_num_conv != -1:
            mo_num_use=mo_num_conv
            mo_coefficient=mo_coefficient[0:mo_num_use]
            mo_occ = mo_occupation[0:mo_num_use]

    # in_convertfort10mol class
    convertfort10mol = Convertfort10mol.parse_from_default_namelist()
    if only_mol:
        convertfort10mol.set_parameter("only_molecular", '.true.', "&control")
    else:
        convertfort10mol.set_parameter("only_molecular", '.false.', "&control")
    convertfort10mol.set_parameter("nmol", mo_num_use - (num_ele_up - num_ele_dn), "&molec_info")
    convertfort10mol.set_parameter("nmolmax", num_ele_dn, "&molec_info")
    convertfort10mol.generate_input()
    convertfort10mol.run()
    
    # rename and remove
    shutil.copy("fort.10_new", "fort.10")
    
    # and read and write MOs to fort.10_new -> fort.10_new_mo
    io_fort10=IO_fort10(f"fort.10")
    
    # check cartesian / spherical flag:
    if ao_cartesian!=0:
        logger.debug("basis set is represent with cartesian")
    else:
        logger.debug("basis set is represent with spherical")
    
    mo_coefficient_turbo=[]
    mo_exponent_turbo=[]
    num=mo_num_use
    #num=1
    for mo_i in range(num):
        logger.debug(f"============{mo_i+1}-th MO================")
        mo_coefficient_list=[]
        mo_exponent_list=[]
        mo_nucleus_index_list=[]
        mo_l_list=[]
        mo_m_list=[]
        #logger.debug(f"i={i}")
        #logger.debug(mo_coefficient[i])
        #logger.debug(ao_shell)
        
        for ao_i, shell_index in enumerate(ao_shell):
            logger.debug(f"--ao_i={ao_i}, shell_index={shell_index}--")
            # for real spherical harmonic notations, 
            # see [https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics]
            
            # initialization
            if ao_i==0:
                mo_coefficient_list_for_reordering=[]
                mo_exponent_list_for_reordering=[]
                mo_nucleus_index_list_for_reordering=[]
                current_ang_mom=-1
    
            # read ang_mom (i.e., angular momentum of the shell)
            ang_mom=basis_shell_ang_mom[shell_index]
            previous_ang_mom=current_ang_mom
            current_ang_mom=ang_mom
            
            # set multiplicity
            multiplicity = 2 * ang_mom + 1
            logger.debug(f"multiplicity = {multiplicity}")
            
            # check if the buffer is null, when ang_mom changes
            if previous_ang_mom != current_ang_mom:
                assert len(mo_coefficient_list_for_reordering) == 0
                assert len(mo_exponent_list_for_reordering) == 0
            
            if current_ang_mom==0: # s shell
                logger.debug("s shell/no permutation is needed.")
                logger.debug("(trexio  notation): s(m=0)")
                logger.debug("(makefun notation): s(m=0)")
                reorder_index=[0]
                reorder_m_list=[0]
                reorder_l_list=[0] * 1
            
            elif current_ang_mom==1: # p shell
                
                if ao_cartesian!=0:
                    logger.debug("p shell/no permutation is needed.")
                    logger.debug("(trexio  notation): px(m=+1), py(m=-1), pz(m=0)")
                    logger.debug("(makefun notation): px(m=+1), py(m=-1), pz(m=0)")
                    reorder_index=[0,1,2]
                    reorder_m_list = [+1,-1,0]
                    reorder_l_list = [1] * 3
                else:
                    logger.debug("p shell/permutation is needed.")
                    logger.debug("(trexio  notation): pz(m=0), px(m=+1), py(m=-1)")
                    logger.debug("(makefun notation): px(m=+1), py(m=-1), pz(m=0)")
                    reorder_index=[1,2,0]
                    reorder_m_list = [+1,-1,0]
                    reorder_l_list = [1] * 3
    
            elif current_ang_mom==2: # d shell
                
                if ao_cartesian!=0:
                    logger.debug("Cartesian notation for d shell is not implemented yet! Sorry.")
                    raise NotImplementedError
                else:
                    logger.debug("d shell/permutation is needed.")
                    logger.debug("(trexio  notation): dz2(m=0), dzx=(m=+1), dyz(m=-1), dx2-y2(m=+2), dxy(m=-2)")
                    logger.debug("(makefun notation): dz2(m=0), dx2-y2(m=+2), dxy(m=-2), dyz(m=-1), dzx=(m=+1)")
                    reorder_index=[0,3,4,2,1]
                    reorder_m_list = [0,+2,-2,-1,+1]
                    reorder_l_list = [2] * 5
                
            elif current_ang_mom==3: # f shell
                
                if ao_cartesian!=0:
                    logger.debug("Cartesian notation for f shell is not implemented yet! Sorry.")
                    raise NotImplementedError
                
                else:
                    logger.debug("f shell/no permutation is needed.")
                    logger.debug("(trexio  notation): f3,0(l=0), f3,+1(l=+1), f3,-1(l=-1), f3,+2(l=+2), f3,-2(l=-2), f3,+3(l=+3), f3,-3(l=-3)")
                    logger.debug("(makefun notation): f3,0(l=0), f3,+1(l=+1), f3,-1(l=-1), f3,+2(l=+2), f3,-2(l=-2), f3,+3(l=+3), f3,-3(l=-3)")
                    reorder_index=[0,1,2,3,4,5,6]
                    reorder_m_list = [0, +1, -1, +2, -2, +3, -3]
                    reorder_l_list = [3] * 7

            elif current_ang_mom == 4:  # g shell

                if ao_cartesian != 0:
                    logger.debug("Cartesian notation for g shell is not implemented yet! Sorry.")
                    raise NotImplemented

                else:
                    logger.debug("g shell/no permutation is needed.")
                    logger.debug(
                        "(trexio  notation): g4,0(l=0), g4,+1(l=+1), g4,-1(l=-1), g4,+2(l=+2), g4,-2(l=-2), g4,+3(l=+3), g4,-3(l=-3), g4,+4(l=+4), g4,-4(l=-4)")
                    logger.debug(
                        "(makefun notation): g4,0(l=0), g4,+1(l=+1), g4,-1(l=-1), g4,+2(l=+2), g4,-2(l=-2), g4,+3(l=+3), g4,-3(l=-3), g4,+4(l=+4), g4,-4(l=-4)")
                    reorder_index = [0, 1, 2, 3, 4, 5, 6, 7, 8]
                    reorder_m_list = [0, +1, -1, +2, -2, +3, -3, +4, -4]
                    reorder_l_list = [4] * 9

            elif current_ang_mom == 5:  # h shell

                if ao_cartesian != 0:
                    logger.debug("Cartesian notation for h shell is not implemented yet! Sorry.")
                    raise NotImplemented

                else:
                    logger.debug("h shell/no permutation is needed.")
                    logger.debug(
                        "(trexio  notation): h5,0(l=0), h5,+1(l=+1), h5,-1(l=-1), h5,+2(l=+2), h5,-2(l=-2), h5,+3(l=+3), h5,-3(l=-3), h5,+4(l=+4), h5,-4(l=-4), h5,+5(l=+5), h5,-5(l=-5)")
                    logger.debug(
                        "(makefun notation): h5,0(l=0), h5,+1(l=+1), h5,-1(l=-1), h5,+2(l=+2), h5,-2(l=-2), h5,+3(l=+3), h5,-3(l=-3), h5,+4(l=+4), h5,-4(l=-4), h5,+5(l=+5), h5,-5(l=-5)")
                    reorder_index = [0, 1, 2, 3, 4, 5, 6, 7, 8 ,9 ,10]
                    reorder_m_list = [0, +1, -1, +2, -2, +3, -3, +4, -4, +5, -5]
                    reorder_l_list = [5] * 11
            
            elif current_ang_mom == 6:  # i shell

                if ao_cartesian != 0:
                    logger.debug("Cartesian notation for i shell is not implemented yet! Sorry.")
                    raise NotImplemented

                else:
                    logger.debug("i shell/no permutation is needed.")
                    logger.debug(
                        "(trexio  notation): i6,0(l=0), i6,+1(l=+1), i6,-1(l=-1), i6,+2(l=+2), i6,-2(l=-2), i6,+3(l=+3), i6,-3(l=-3), i6,+4(l=+4), i6,-4(l=-4), i6,+5(l=+5), i6,-5(l=-5), i6,+6(l=+6), i6,-6(l=-6)")
                    logger.debug(
                        "(makefun notation): i6,0(l=0), i6,+1(l=+1), i6,-1(l=-1), i6,+2(l=+2), i6,-2(l=-2), i6,+3(l=+3), i6,-3(l=-3), i6,+4(l=+4), i6,-4(l=-4), i6,+5(l=+5), i6,-5(l=-5), i6,+6(l=+6), i6,-6(l=-6)")
                    reorder_index = [0, 1, 2, 3, 4, 5, 6, 7, 8 ,9 ,10, 11, 12]
                    reorder_m_list = [0, +1, -1, +2, -2, +3, -3, +4, -4, +5, -5, +6, -6]
                    reorder_l_list = [6] * 13

            else:
                logger.error(f" Angular momentum={ang_mom} is not implemented in TurboRVB!!!")
                raise NotImplementedError
                
            basis_index_list = [n for n, v in enumerate(basis_shell_index) if v == shell_index]
            
            mo_coeff_local_buffer=[]
            mo_exponent_local_buffer = []
            mo_nucleus_index_local_buffer = []
    
            logger.debug(basis_index_list)
            logger.debug(basis_nucleus_index)
            for basis_index in basis_index_list:
                mo_coeff_local_buffer.append(mo_coefficient[mo_i][ao_i] * basis_coefficient[basis_index])
                mo_exponent_local_buffer.append(basis_exponent[basis_index])
                mo_nucleus_index_local_buffer.append(basis_nucleus_index[shell_index])
    
            mo_coefficient_list_for_reordering.append(mo_coeff_local_buffer)
            mo_exponent_list_for_reordering.append(mo_exponent_local_buffer)
            mo_nucleus_index_list_for_reordering.append(mo_nucleus_index_local_buffer)
    
            # write MOs!!
            if len(mo_coefficient_list_for_reordering) == multiplicity:
                logger.debug("--write MOs!!--")
                mo_coeff_list_reordered = [mo_coefficient_list_for_reordering[i] for i in reorder_index]
                mo_exp_list_reordered = [mo_exponent_list_for_reordering[i] for i in reorder_index]
                mo_nucleus_index_list_reordered = [mo_nucleus_index_list_for_reordering[i] for i in reorder_index]
                mo_m_list_reordered = reorder_m_list * len(basis_index_list)
                mo_l_list_reordered = reorder_l_list * len(basis_index_list)
    
                # reordered again is needed here for contracted shell!
                # because order is not px ... py ... pz but px,py,pz,px,py,pz ...
                for ii in range(len(mo_coeff_list_reordered[0])):
                    for mo_coeff_reordered in mo_coeff_list_reordered:
                        mo_coefficient_list.append(mo_coeff_reordered[ii])
                for ii in range(len(mo_exp_list_reordered[0])):
                    for mo_exp_reordered in mo_exp_list_reordered:
                        mo_exponent_list.append(mo_exp_reordered[ii])
                for ii in range(len(mo_nucleus_index_list_reordered[0])):
                    for mo_nucleus_index_reordered in mo_nucleus_index_list_reordered:
                        mo_nucleus_index_list.append(mo_nucleus_index_reordered[ii])
                """
                for ii in range(len(mo_m_list_reordered[0])):
                    for mo_m_reordered in mo_m_list_reordered:
                        mo_m_list.append(mo_m_reordered[ii])
                for ii in range(len(mo_angmom_list_reordered[0])):
                    for mo_angmom_reordered in mo_angmom_list_reordered:
                        mo_l_list.append(mo_angmom_reordered[ii])
                """
                mo_m_list += mo_m_list_reordered
                mo_l_list += mo_l_list_reordered
    
                # reset buffer after writing MOs
                mo_coefficient_list_for_reordering=[]
                mo_exponent_list_for_reordering=[]
                mo_nucleus_index_list_for_reordering=[]
    
        logger.debug(len(mo_exponent_list))
        logger.debug(len(mo_coefficient_list))
        logger.debug(len(mo_nucleus_index_list))
        logger.debug(len(mo_m_list))
        logger.debug(len(mo_l_list))
    
        assert len(mo_coefficient_list) == len(mo_exponent_list)
        assert len(mo_coefficient_list) == len(mo_nucleus_index_list)
        assert len(mo_coefficient_list) == len(mo_m_list)
        assert len(mo_coefficient_list) == len(mo_l_list)
    
        logger.debug('========================================')
        logger.debug('End buffering MOs')
        logger.debug('========================================')
    
        # remove duplicated exponents!
        # Note: TurboRVB internally removes duplicated exponents.
        logger.debug(mo_nucleus_index_list)
        logger.debug(len(mo_nucleus_index_list))
    
        ref_exponent_index_list=[]
        removed_exponent_index_list=[]
        for unique_mo_nucleus_index in set(mo_nucleus_index_list):
            logger.debug(f"unique_mo_nucleus_index = {unique_mo_nucleus_index}.")
            mo_nucleus_index=[i for i, x in enumerate(mo_nucleus_index_list) if x == unique_mo_nucleus_index]
            logger.debug(mo_nucleus_index)
            mo_exponent_list_nucleus=[mo_exponent_list[i] for i in mo_nucleus_index]
            mo_m_list_nucleus=[mo_m_list[i] for i in mo_nucleus_index]
            mo_l_list_nucleus=[mo_l_list[i] for i in mo_nucleus_index]
            logger.debug(mo_exponent_list_nucleus)
    
            unique_exp_list=list(set(mo_exponent_list_nucleus))
            logger.debug(mo_exponent_list_nucleus)
            logger.debug(unique_exp_list)
    
            logger.debug(len(mo_exponent_list_nucleus))
            logger.debug(len(unique_exp_list))
            if len(unique_exp_list) != len(mo_exponent_list_nucleus):
                logger.debug(f"Duplicated exponents exist for nucleus index = {unique_mo_nucleus_index}.")
                for unique_exp in unique_exp_list:
                    unique_exponent_index=[i for i, x in enumerate(mo_exponent_list_nucleus) if x == unique_exp]
                    mo_e_list_u = [mo_exponent_list_nucleus[i] for i in unique_exponent_index]
                    assert len(set(mo_e_list_u)) == 1
                    mo_l_list_u = [mo_l_list_nucleus[i] for i in unique_exponent_index]
                    mo_m_list_u = [mo_m_list_nucleus[i] for i in unique_exponent_index]
    
                    while len(unique_exponent_index) > 0 and len(mo_l_list_u) > 0 and len(mo_m_list_u) > 0:
                        ref_flag = False
                        removed_uu_list = []
                        removed_exponent_index_list_ = []
                        logger.debug(unique_exponent_index)
                        logger.debug(mo_l_list_nucleus)
                        logger.debug(mo_m_list_nucleus)
                        for uu, (e_i, mo_l, mo_m) in enumerate(zip(unique_exponent_index, mo_l_list_u, mo_m_list_u)):
                            if not ref_flag:
                                ref_e_i=e_i
                                mo_l_ref=mo_l
                                mo_m_ref=mo_m
                                removed_uu_list.append(uu)
                                ref_flag=True
                            else:
                                if mo_l_ref == mo_l and mo_m_ref == mo_m:
                                    removed_uu_list.append(uu)
                                    removed_exponent_index_list_.append(e_i)
    
                        logger.debug(removed_uu_list)
                        for kk in reversed(sorted(removed_uu_list)):
                            unique_exponent_index.pop(kk)
                            mo_l_list_u.pop(kk)
                            mo_m_list_u.pop(kk)
    
                        logger.debug(ref_e_i)
                        logger.debug(removed_exponent_index_list_)
                        logger.debug(mo_nucleus_index[ref_e_i])
                        logger.debug([mo_nucleus_index[pp] for pp in removed_exponent_index_list_])
                        ref_exponent_index_list.append(mo_nucleus_index[ref_e_i])
                        removed_exponent_index_list.append([mo_nucleus_index[pp] for pp in removed_exponent_index_list_])
    
            else:
                logger.debug(f"No duplicated exponent is found for nucleus index ={unique_mo_nucleus_index}")
    
        # Finally, here, duplicated exponents are removed.
        # because of using pop here.
        # if pop is used inside the above for loop, indices change during the for loop.
        pop_index_list = []
        for ref_exponent_index, removed_exponent_index in zip(ref_exponent_index_list, removed_exponent_index_list):
            for r_index in reversed(sorted(removed_exponent_index)):
                mo_coefficient_list[ref_exponent_index] += mo_coefficient_list[r_index]
                pop_index_list.append(r_index)
        logger.debug(sorted(pop_index_list))
        logger.debug(len(sorted(pop_index_list)))
        for pop_index in reversed(sorted(pop_index_list)):
            mo_coefficient_list.pop(pop_index)
    
        mo_coefficient_turbo.append(mo_coefficient_list)
        mo_exponent_turbo.append(mo_exponent_list)

    # fort.10 MO replace
    logger.debug(len(io_fort10.f10detbasissets.mo_coefficient[0]))
    logger.debug(len(mo_coefficient_turbo[0]))

    # molecular orbital swapped, spin polarized cases.
    mo_coefficient_turbo_unpaired=[]
    pop_lst=[num_ele_dn+nel_diff for nel_diff in range(0, num_ele_up - num_ele_dn)]
    for p in reversed(pop_lst):
        logger.debug("molecular orbital swapped")
        mo_coefficient_turbo_unpaired.append(mo_coefficient_turbo.pop(p))
    for m in reversed(mo_coefficient_turbo_unpaired):
        mo_coefficient_turbo.append(m)

    # fort.10 MO replace
    logger.info("Writing obtained MOs to fort.10....(It might take a while).")
    if not complex_flag:
        io_fort10.f10detbasissets.mo_coefficient = mo_coefficient_turbo
    else:
        # complex case,
        # mo_coefficient_turbo -> mo_coefficient_turbo_real, mo_coefficient_turbo_imag
        mo_coefficient_turbo_real = []
        mo_coefficient_turbo_imag   = []
        for mo__ in mo_coefficient_turbo:
            mo_real_b = []
            mo_imag_b = []
            for coeff in mo__:
                mo_real_b.append(coeff.real)
                mo_imag_b.append(coeff.imag)

            mo_real_b_up=list(np.array(mo_real_b) * +1)
            mo_imag_b_up=list(np.array(mo_imag_b) * +1) # up phase is +
            # these commented lines are wrong!! In general, the wf does not symmetric with respect to the time reversal except for TRIM points.
            #mo_real_b_dn=list(np.array(mo_real_b) * +1) # real part is the same as the up spin.
            #mo_imag_b_dn=list(np.array(mo_imag_b) * -1) # dn phase is -. because the opposite phase is attached in turbo with option double k-grid=.true.
            mo_real_b_dn = list(np.array(mo_real_b) * +1)
            mo_imag_b_dn = list(np.array(mo_imag_b) * +1)

            # turbo seems to treat up and dn MOs individually for complex cases!!
            # trexio currently supports only restricted open or closed shells, => This is a TODO
            # so, here just duplicate the MOs. => is it ok??
            # Yes. the order seems correct, i.e, up, dn, up, dn....
            # but, let me see,,, for spin-polarized cases??
            if num_ele_up != num_ele_dn:
                logger.error("spin-polarized case for a complex system is not implemented yet.")
                raise NotImplementedError

            mo_coefficient_turbo_real.append(mo_real_b_up) # up
            mo_coefficient_turbo_imag.append(mo_imag_b_up) # up
            mo_coefficient_turbo_real.append(mo_real_b_dn) # dn
            mo_coefficient_turbo_imag.append(mo_imag_b_dn) # dn

        logger.debug(mo_coefficient_turbo)
        logger.debug(len(io_fort10.f10detbasissets.mo_coefficient))
        logger.debug(len(mo_coefficient_turbo_real))
        logger.debug(len(io_fort10.f10detbasissets.mo_coefficient_imag))
        logger.debug(len(mo_coefficient_turbo_imag))
        io_fort10.f10detbasissets.mo_coefficient = mo_coefficient_turbo_real
        io_fort10.f10detbasissets.mo_coefficient_imag = mo_coefficient_turbo_imag
    
    # clean up
    if cleanup:
        logger.info("Cleaning up temporary files.")
        files=["fort.10_in","fort.10_new","out_make","out_mol","symmetries.dat"]
        for file_ in files:
            if os.path.isfile(file_): os.remove(file_)

    logger.debug("END of conversion")

def main():
    # parser.add_argument
    from database_setup import all_electron_basis_set_list, ccECP_basis_set_list, BFD_basis_set_list
    from utils_workflows.utility import prompt
    from database_setup import database_setup

    parser = argparse.ArgumentParser(description='This program is a python-based script for converting a TREXIO file to a TurboRVB Wavefunction file')
    parser.add_argument('trexio_file', help='Name of TREXIO file')
    parser.add_argument("-jasbasis", "--jas_basis_sets",
                  help=f'Specify a basis set for the Jastrow part {all_electron_basis_set_list}',
                  default=None,
                  choices=all_electron_basis_set_list,
                  )
    parser.add_argument("-jascont", "--jas_contracted_flag",
                  help='Contraction flag for the jastrow part',
                  default=False,
                  action="store_true"
                  )
    parser.add_argument("-jascutbasis", "--jas_cut_basis_option",
                  help='Cutting the jastrow basis set according to a default criteria',
                  default=False,
                  action="store_true"
                  )
    parser.add_argument("-twist", "--twist_average",
                        help="flag for twist average",
                        default=False,
                        action="store_true"
                        )
    parser.add_argument('-log', '--loglevel', help='logger setlevel', default="INFO", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])
    parser.add_argument('-c', '--cleanup', help='clean up temporary files', default=True, action="store_true")
    args = parser.parse_args()

    logger = getLogger('Turbo-Genius').getChild(__name__)
    logger.setLevel(args.loglevel)
    stream_handler = StreamHandler()
    stream_handler.setLevel(args.loglevel)
    if args.loglevel in {'DEBUG'}:
        handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    else:
        handler_format = Formatter('%(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    logger.info("    ")
    logger.info("Program   : TREXIO-TurboRVB converter")
    logger.info("Author    : 2022 Kosuke Nakano (SISSA/JAIST)")
    logger.info("E-mail    : kousuke_1123@icloud.com")
    logger.info("Website   : https://www.kosuke-nakano-research.info")
    logger.info("License   : MIT lisence")
    logger.info("    ")
    logger.info("Please cite the following paper when publishing your work(s).")
    logger.info("    ")
    logger.info("TurboRVB:")
    logger.info("    K. Nakano et al. J. Chem. Phys., 152, 204121 (2020)")
    logger.info("Turbo-Genius:")
    logger.info("    K. Nakano et al. in preparation.")
    logger.info("TREXIO:")
    logger.info("    TREX project, in preparation.")
    logger.info("    ")

    # twist average setting
    if args.twist_average:
        with open(os.path.join(os.getcwd(), "kp_info.dat"), "r") as f:
            lines = f.readlines()
        k_num = len(lines)
    else:
        k_num = 1

    for num in range(k_num):
        trexio_r = Trexio_wrapper_r(trexio_file=args.trexio_file)
        element_symbols = trexio_r.labels_r

        # jastrow setting
        if args.jas_basis_sets is not None:
            database_setup(database="BSE")
            jas_basis_files = [];
            jas_basis_choice = {}

            def database_founder(data_sets_list, element, data_choice, prefix="basis_set"):
                if len(data_sets_list) == 0:
                    logger.error(f"The chosen {prefix} is not found in the database!!")
                    raise NotImplementedError
                elif len(data_sets_list) == 1:
                    data_set_found = data_sets_list[0]
                    logger.info(f"The chosen {prefix} is found, {os.path.basename(data_set_found)}")
                    return data_set_found, data_choice
                else:  # >= 2
                    if element not in data_choice.keys():
                        logger.info(f"More than two {prefix}s are found!")

                        def checker(choice):
                            try:
                                if int(choice) in range(len(data_sets_list)):
                                    return True
                                else:
                                    return False
                            except ValueError:
                                return False

                        b_list_shown = [f'{i}:{os.path.basename(d)}' for i, d in enumerate(data_sets_list)]
                        b_index = int(prompt(f"Choose one of them, 0,1,.. from {b_list_shown}:", checker=checker))
                        data_set_found = data_sets_list[b_index]
                        data_choice[element] = data_set_found
                        logger.info(f"The chosen {prefix} is {os.path.basename(data_set_found)}")
                        return data_set_found, data_choice

                    else:
                        data_set_found = data_choice[element]
                        logger.info(f"The chosen {prefix} is found, {os.path.basename(data_set_found)}")
                        return data_set_found, data_choice

            # jas. basis set
            for element in element_symbols:
                jas_basis_sets_list = glob.glob(
                    os.path.join(turbo_genius_tmp_dir, "basis_set", "BSE", f"{element}_{args.jas_basis_sets}*.basis"))
                logger.debug(jas_basis_sets_list)
                jas_basis_chosen, jas_basis_choice = database_founder(data_sets_list=jas_basis_sets_list,
                                                                      element=element, data_choice=jas_basis_choice,
                                                                      prefix="basis_set")
                jas_basis_files.append(jas_basis_chosen)
            jas_basis_sets = Jas_Basis_sets.parse_basis_sets_from_gamess_format_files(files=jas_basis_files)

            if not args.jas_contracted_flag:
                jas_basis_sets.contracted_to_uncontracted()

            if args.jas_cut_basis_option:
                # cut basis, jas_basis, according to max criteria, exponents > max (det part)
                for nuc, element in enumerate(element_symbols):
                    #thr_exp = 8 * return_atomic_number(element) ** 2
                    thr_exp = 4 * return_atomic_number(element) # not 8*Z**2 but 4*Z
                    jas_basis_sets.cut_orbitals(thr_exp=thr_exp, nucleus_index=nuc, method="larger")
                    thr_angmom = jas_basis_sets.get_largest_angmom(nucleus_index=nuc)
                    jas_basis_sets.cut_orbitals(thr_angmom=thr_angmom, nucleus_index=nuc, method="larger-angmom")

        # jastrow is None
        else:
            jas_basis_sets = Jas_Basis_sets()

        ## trexio -> turborvb_wf
        # conversion
        trexio_to_turborvb_wf(trexio_file=args.trexio_file, cleanup=args.cleanup, max_occ_conv=1.0e-3,
                              jas_basis_sets=jas_basis_sets)

        if args.twist_average:
            turborvb_scratch_dir = os.path.join(os.getcwd(), "turborvb.scratch")
            os.makedirs(turborvb_scratch_dir, exist_ok=True)
            shutil.move(os.path.join(trexio_dir, "fort.10"), os.path.join(turborvb_scratch_dir, "fort.10_{:0>6}".format(num)))

    if args.twist_average:
        shutil.copy(os.path.join(os.getcwd(), "fort.10_{:0>6}".format(0)), os.path.join(os.getcwd(), "fort.10"))


if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("INFO")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    main()

    """
    # moved to examples
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))

    os.chdir(os.path.join("../../", "tests", "trexio_to_turborvb"))
    trexio_file="diamond_trexio_q.hdf5"

    trexio_to_turborvb_wf(trexio_file=trexio_file,
                          jas_basis_sets=Jas_Basis_sets(),
                          max_occ_conv=1.0e-4,  # maximum occ used for the conv, not used with mo_num
                          mo_num_conv=-1,  # num mo used for the conv, not used with max occ
                          cleanup=True)
    """