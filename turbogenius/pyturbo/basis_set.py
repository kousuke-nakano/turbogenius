#!python -u
# -*- coding: utf-8 -*-

#python modules
import os, sys
import re
import numpy as np

#logger set
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo').getChild(__name__)
stream_handler = StreamHandler()

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils.utility import return_ang_mom, turbo_prim_orb_type_num, turbo_cont_orb_type_num, return_orbchr
from utils.env import pyturbo_root
from utils.downloader import BSE, BFD, ccECP

# we should consider also complex cases... gosh.
class Basis_sets:

    def __init__(self,
                nucleus_index=[],
                shell_ang_mom=[],
                shell_ang_mom_turbo_notation=[],
                shell_factor=[],
                shell_index=[],
                exponent=[],
                coefficient=[],
                coefficient_imag=[],
                prim_factor=[],
    ):

        self.nucleus_index = nucleus_index
        self.shell_ang_mom = shell_ang_mom
        self.shell_ang_mom_turbo_notation = shell_ang_mom_turbo_notation
        self.shell_factor = shell_factor
        self.shell_index = shell_index
        self.exponent = exponent
        self.coefficient = coefficient
        self.coefficient_imag = coefficient_imag
        self.prim_factor = prim_factor

        logger.debug(self.nucleus_index)
        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_ang_mom_turbo_notation)
        logger.debug(self.shell_factor)
        logger.debug(self.shell_index)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.coefficient_imag)
        logger.debug(self.prim_factor)

        # assertion!!
        assert len(self.exponent) == len(self.coefficient)
        assert len(self.exponent) == len(self.prim_factor)
        assert len(self.exponent) == len(self.shell_index)

        assert len(self.nucleus_index) == len(self.shell_ang_mom)
        assert len(self.nucleus_index) == len(self.shell_ang_mom_turbo_notation)
        assert len(self.nucleus_index) == len(self.shell_factor)

        assert len(self.coefficient) == len(self.coefficient_imag) or len(self.coefficient_imag) == 0

    @property
    def prim_num(self):
        return len(self.exponent)
    @property
    def shell_num(self):
        return len(self.shell_ang_mom)
    @property
    def nuclei_num(self):
        return len(set(self.nucleus_index))
    @property
    def complex_flag(self):
        if len(self.coefficient_imag) == 0:
            return False
        else:
            return True

    def cut_exponents(self, thr_exp, nucleus_index=None, method="larger"):
        # method = larger, smaller, equil

        logger.debug("--Before cut--")
        logger.debug(self.nucleus_index)
        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_ang_mom_turbo_notation)
        logger.debug(self.shell_factor)
        logger.debug(self.shell_index)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.coefficient_imag)
        logger.debug(self.prim_factor)

        if nucleus_index is None:
            nucleus_index=list(set(self.nucleus_index))
        elif type(nucleus_index) == int:
            nucleus_index = [nucleus_index]
        elif type(nucleus_index) == list:
            pass
        else:
            logger.error("type(nucleus_index) seems wrong.")
            raise NotImplementerErorr

        cut_prim_index=[]
        for nuc in nucleus_index:
            shell_index_list = [i for i, x in enumerate(self.nucleus_index) if x == nuc]
            logger.debug(shell_index_list)
            for shell in shell_index_list:
                prim_index_list = [i for i, x in enumerate(self.shell_index) if x == shell]
                logger.debug(prim_index_list)
                for prim_index in prim_index_list:
                    if method=="larger":
                        if self.exponent[prim_index] >= thr_exp:
                            cut_prim_index.append(prim_index)
                    elif method=="smaller":
                        if self.exponent[prim_index] <= thr_exp:
                            cut_prim_index.append(prim_index)
                    elif method=="equil":
                        if self.exponent[prim_index] == thr_exp:
                            cut_prim_index.append(prim_index)
                    else:
                        logger.error(f"Not implemented method={method}")
                        raise NotImplementerErorr

        for prim_index in reversed(cut_prim_index):
            self.remove_primitive_basis(prim_index=prim_index)

        logger.debug("--After cut--")
        logger.debug(self.nucleus_index)
        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_ang_mom_turbo_notation)
        logger.debug(self.shell_factor)
        logger.debug(self.shell_index)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.coefficient_imag)
        logger.debug(self.prim_factor)

    def remove_primitive_basis(self, prim_index):
        logger.debug(f"prim_index={prim_index} will be removed.")

        logger.debug("--Before remove basis--")
        logger.debug(self.nucleus_index)
        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_ang_mom_turbo_notation)
        logger.debug(self.shell_factor)
        logger.debug(self.shell_index)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.coefficient_imag)
        logger.debug(self.prim_factor)

        #self.nucleus_index = nucleus_index
        #self.shell_ang_mom = shell_ang_mom
        #self.shell_ang_mom_turbo_notation = shell_ang_mom_turbo_notation
        #self.shell_factor = shell_factor
        #self.shell_index = shell_index

        self.exponent.pop(prim_index)
        self.coefficient.pop(prim_index)
        if self.complex_flag: self.coefficient_imag.pop(prim_index)
        self.prim_factor.pop(prim_index)

        # whether or not the corresponding shell should be removed.
        shell_index=self.shell_index.pop(prim_index)

        if shell_index in self.shell_index:
            logger.debug("Removed a primitive basis in a contracted shell")
        else:
            logger.debug("Removed a primitive shell")

            self.nucleus_index.pop(shell_index)
            self.shell_ang_mom.pop(shell_index)
            self.shell_ang_mom_turbo_notation.pop(shell_index)
            self.shell_factor.pop(shell_index)

            for i in range(len(self.shell_index[prim_index:])):
                self.shell_index[prim_index+i] -= 1

        logger.debug("--After remove basis--")
        logger.debug(self.nucleus_index)
        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_ang_mom_turbo_notation)
        logger.debug(self.shell_factor)
        logger.debug(self.shell_index)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.coefficient_imag)
        logger.debug(self.prim_factor)

    def contracted_to_uncontracted(self):
        logger.info("Conversion of basis sets, contracted -> uncontracted")
        logger.debug("--Before conversion--")
        logger.debug(self.nucleus_index)
        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_ang_mom_turbo_notation)
        logger.debug(self.shell_factor)
        logger.debug(self.shell_index)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.coefficient_imag)
        logger.debug(self.prim_factor)

        assert len(self.exponent) == len(self.coefficient)
        nucleus_index=[]
        shell_ang_mom=[]
        shell_ang_mom_turbo_notation=[]
        shell_factor=[]
        shell_index=[]
        exponent=[]
        coefficient=[]
        coefficient_imag=[]
        prim_factor=[]

        for nuc in range(len(set(self.nucleus_index))):
            exponent_nuc=[]
            shell_ang_mom_nuc=[]
            shell_index_list=[i for i, x in enumerate(self.nucleus_index) if x == nuc]

            for i_shell in shell_index_list:
                prim_index_list = [i for i, x in enumerate(self.shell_index) if x == i_shell]
                exponent_n=[self.exponent[i] for i in prim_index_list]
                coefficient_n=[self.coefficient[i] for i in prim_index_list]
                if self.complex_flag:
                    coefficient_imag_n=[self.coefficient_imag[i] for i in prim_index_list]
                prim_factor_n=[self.prim_factor[i] for i in prim_index_list]

                shell_ang_mom_n = self.shell_ang_mom[i_shell]
                shell_ang_mom_turbo_notation_n = self.shell_ang_mom_turbo_notation[i_shell]
                shell_factor_n = self.shell_factor[i_shell]

                assert len(exponent_n) == len(coefficient_n)
                prim_num=len(exponent_n)
                for p in range(prim_num):
                    if exponent_n[p] not in exponent_nuc:
                        store=True
                    else: # exponent_n[p] in exponent:
                        e_index_list=[i for i, x in enumerate(exponent_nuc) if x == exponent_n[p]]
                        if all([shell_ang_mom_n != shell_ang_mom_nuc[e] for e in e_index_list]):
                            store=True
                        else:
                            store=False

                    if store:
                        exponent_nuc.append(exponent_n[p]) # prim
                        shell_ang_mom_nuc.append(shell_ang_mom_n) #shell

                        #prim
                        nucleus_index.append(nuc)
                        exponent.append(exponent_n[p])
                        coefficient.append(1.0)
                        if self.complex_flag:
                            coefficient_imag.append(0.0)
                        prim_factor.append(prim_factor_n[p])
                        shell_index.append(len(shell_index))

                        #shell
                        shell_factor.append(None)
                        shell_ang_mom.append(shell_ang_mom_n)
                        shell_ang_mom_turbo_notation.append(
                            turbo_prim_orb_type_num(return_orbchr(shell_ang_mom_n))
                        )

        self.nucleus_index=nucleus_index
        self.exponent=exponent
        self.coefficient=coefficient
        self.coefficient_imag = coefficient_imag
        self.prim_factor=prim_factor
        self.shell_index=shell_index
        self.shell_factor=shell_factor
        self.shell_ang_mom=shell_ang_mom
        self.shell_ang_mom_turbo_notation=shell_ang_mom_turbo_notation

        logger.debug("--After conversion--")
        logger.debug(self.nucleus_index)
        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_ang_mom_turbo_notation)
        logger.debug(self.shell_factor)
        logger.debug(self.shell_index)
        logger.debug(self.exponent)
        logger.debug(self.coefficient)
        logger.debug(self.coefficient_imag)
        logger.debug(self.prim_factor)

    @classmethod
    def parse_basis_sets_from_gamess_format_files(cls, files=[]):
        texts=[]
        for file in files:
            with open(file, 'r') as f:
                text_buffer=f.readlines()
            text = "".join(text_buffer)
            texts.append(text)

        return cls.parse_basis_sets_from_texts(texts=texts, format="gamess")

    @classmethod
    def parse_basis_sets_from_eCEPP_format_files(cls, files=[]):
        texts=[]
        for file in files:
            with open(file, 'r') as f:
                text_buffer=f.readlines()
            text = "".join(text_buffer)
            texts.append(text)

        return cls.parse_basis_sets_from_texts(texts=texts, format="eCEPP")

    @classmethod
    def parse_basis_sets_from_texts(cls, texts=[], format="gamess"):
        nucleus_index = []
        shell_ang_mom = []
        shell_ang_mom_turbo_notation = []
        shell_factor = []
        shell_index = []
        exponent = []
        coefficient = []
        coefficient_imag = []
        prim_factor = []

        shell_num=0
        prim_num=0

        for nuc_i, text in enumerate(texts):

            if format=="gamess":
                basis_set=Basis_set.parse_basis_set_info_from_gamess_format_text(text)
            elif format=="eCEPP":
                basis_set=Basis_set.parse_basis_set_info_from_eCEPP_format_text(text)
            else:
                logger.error(f"format = {format} is not implemented")
                raise NotImplementedError
            # storing
            nucleus_index+=[nuc_i]*basis_set.shell_num
            shell_ang_mom+=basis_set.shell_ang_mom
            shell_factor+=[1.0]*basis_set.shell_num
            shell_index += list(np.array(basis_set.shell_index) + shell_num)
            exponent+=basis_set.exponent_list
            coefficient+=basis_set.coefficient_list
            prim_factor+=[1.0]*basis_set.prim_num # for the time being. / this can be computed analytically.

            shell_num+=basis_set.shell_num
            prim_num+=basis_set.prim_num

        for i, ang_mom in enumerate(shell_ang_mom):
            if shell_index.count(i) > 1:
                contraction = True
            else:
                contraction = False

            if contraction:
                shell_ang_mom_turbo_notation.append(turbo_cont_orb_type_num(orb_type_chr=return_orbchr(ang_mom=ang_mom)))
            else:
                shell_ang_mom_turbo_notation.append(turbo_prim_orb_type_num(orb_type_chr=return_orbchr(ang_mom=ang_mom)))

        return cls(
            nucleus_index = nucleus_index,
            shell_ang_mom = shell_ang_mom,
            shell_ang_mom_turbo_notation = shell_ang_mom_turbo_notation,
            shell_factor = shell_factor,
            shell_index = shell_index,
            exponent = exponent,
            coefficient = coefficient,
            coefficient_imag = coefficient_imag,
            prim_factor = prim_factor
        )

class Basis_set:

    def __init__(self,
                 shell_ang_mom = [],
                 shell_index = [],
                 exponent_list = [],
                 coefficient_list = [],
                 coefficient_imag_list = [],
                 ):

        self.shell_ang_mom = shell_ang_mom
        self.shell_index = shell_index
        self.exponent_list = exponent_list
        self.coefficient_list = coefficient_list
        self.coefficient_imag_list = coefficient_imag_list

        logger.debug(self.shell_ang_mom)
        logger.debug(self.shell_index)
        logger.debug(self.exponent_list)
        logger.debug(self.coefficient_list)
        logger.debug(self.coefficient_imag_list)

        # assertion!
        assert len(self.exponent_list) == len(self.coefficient_list)
        assert len(self.exponent_list) == len(self.shell_index)
        assert len(self.coefficient_list) == len(self.coefficient_imag_list) or len(self.coefficient_imag_list) == 0


    @property
    def prim_num(self):
        return len(self.exponent_list)
    @property
    def shell_num(self):
        return len(self.shell_ang_mom)
    @property
    def complex_flag(self):
        if len(self.coefficient_imag) == 0:
            return False
        else:
            return True

    @classmethod
    def parse_basis_set_info_from_gamess_format_file(cls, file):
        with open(file, 'r') as f:
            text = f.readlines()
        text="".join(text)
        return cls.parse_basis_set_info_from_gamess_format_text(text=text)

    @classmethod
    def parse_basis_set_info_from_eCEPP_format_file(cls, file):
        with open(file, 'r') as f:
            text = f.readlines()
        text="".join(text)
        return cls.parse_basis_set_info_from_eCEPP_format_text(text=text)

    @classmethod
    def parse_basis_set_info_from_eCEPP_format_text(cls, text):
        shell_ang_mom=[]
        shell_index=[]
        exponent_list=[]
        coefficient_list=[]
        coefficient_imag_list = []

        lines = text.split("\n")
        lines = [line for line in lines if re.match(r'^\s*\S+.*', line)] # remove blank lines

        logger.info(lines)
        for line in lines:
            if re.match(r"^!.*", line):
                m=re.match(r"^!([a-z]).*", line)
                l_str=m.group(1)
                logger.debug(f"ang_mom = {l_str}")
                ang_mom=return_ang_mom(orb_typ_chr=l_str)
            else:
                s_line=(line.rstrip(";").split(","))
                fl=s_line[0]
                if fl==l_str:
                    logger.debug("Exponent line")
                    r_exponent_list=list(map(float,s_line[2:]))
                    logger.debug(r_exponent_list)
                    for r_exponent in r_exponent_list:
                        exponent_list.append(r_exponent)
                elif fl=="c":
                    logger.debug("Contraction coeff. line")
                    r_coefficient_list = list(map(float, s_line[2:]))
                    logger.debug(r_coefficient_list)
                    r_shell_index = len(shell_ang_mom)
                    shell_ang_mom.append(ang_mom)
                    for r_coefficient in r_coefficient_list:
                        shell_index.append(r_shell_index)
                        coefficient_list.append(r_coefficient)

        return cls(shell_ang_mom = shell_ang_mom,
                   shell_index = shell_index,
                   exponent_list = exponent_list,
                   coefficient_list = coefficient_list,
                   coefficient_imag_list = coefficient_imag_list
                   )

    @classmethod
    def parse_basis_set_info_from_gamess_format_text(cls, text):
        shell_ang_mom=[]
        shell_index=[]
        exponent_list=[]
        coefficient_list=[]
        coefficient_imag_list=[]

        #index_num = [n for n, v in enumerate(text) if v == "****\n"]
        #orbital = text[index_num[0] + 2:index_num[1]]
        lines=text.split("\n")
        lines = [line for line in lines if re.match(r'^\s*\S+.*', line)] # remove blank lines

        # reading the orbitals
        i = 0
        shell_num=0
        while i < len(lines):
            line = lines[i].split()
            orb_type = str(line[0]).lower()
            #logger.info(orb_type)
            num_contraction = int(line[1])

            # uncontracted basis
            if num_contraction == 1:
                orb_type_chr = str(orb_type)
                ang_mom = int(return_ang_mom(orb_type_chr))

                if len(lines[i + 1].split()) == 3:
                    # gamess format
                    stride = 1
                elif len(lines[i + 1].split()) == 2:
                    # gaussian format
                    stride = 0
                else:
                    raise ValueError("The pseudo potential format seems wrong")

                exponent = float(lines[i + 1].split()[stride+0].replace("D", "E"))
                coefficient = float(lines[i + 1].split()[stride+1].replace("D", "E"))

                shell_ang_mom.append(ang_mom)
                shell_index.append(shell_num)
                exponent_list.append(exponent)
                coefficient_list.append(coefficient)

                i = i + 1 + 1
                shell_num+=1

            # contracted basis
            else:
                exponent_list_ = []
                coefficient_list_ = []

                orb_type_chr = str(orb_type)
                ang_mom = int(return_ang_mom(orb_type_chr))

                shell_ang_mom.append(ang_mom)

                for j in range(1, num_contraction + 1):

                    if len(lines[i + j].split()) == 3:
                        # gamess format
                        stride = 1
                    elif len(lines[i + j].split()) == 2:
                        # gaussian format
                        stride = 0
                    else:
                        raise ValueError("The pseudo potential format seems wrong")

                    exponent_list_.append(float(float(lines[i + j].split()[stride+0].replace("D", "E"))))
                    coefficient_list_.append(float(float(lines[i + j].split()[stride+1].replace("D", "E"))))
                    shell_index.append(shell_num)

                exponent_list+=exponent_list_
                coefficient_list+=coefficient_list_

                i = i + j + 1
                shell_num += 1

        return cls(shell_ang_mom = shell_ang_mom,
                   shell_index = shell_index,
                   exponent_list = exponent_list,
                   coefficient_list = coefficient_list,
                   coefficient_imag_list = coefficient_imag_list
                   )

    def to_text_gamess_format(self):
        text=""
        for yy, ang_mom in enumerate(self.shell_ang_mom):
            t_shell_index=[i for i,x in enumerate(self.shell_index) if x==yy]
            line=f"{return_orbchr(ang_mom=ang_mom).upper()}   {len(t_shell_index)}\n"
            text+=line

            t_exponent_list = [self.exponent_list[i] for i in t_shell_index]
            t_coefficient_list = [self.coefficient_list[i] for i in t_shell_index]
            assert len(t_exponent_list) == len(t_coefficient_list)

            for kk, (exponent, coefficient) in enumerate(zip(t_exponent_list, t_coefficient_list)):
                line = f" {kk+1:d}    {exponent:.8f}   {coefficient:8f}\n"
                text += line

        return text

    #todo! this class should hold element attribute.
    def to_text_nwchem_format(self, element):
        text=""
        for yy, ang_mom in enumerate(self.shell_ang_mom):
            t_shell_index=[i for i,x in enumerate(self.shell_index) if x==yy]
            line=f"{element} {return_orbchr(ang_mom=ang_mom).upper()}\n"
            text+=line

            t_exponent_list = [self.exponent_list[i] for i in t_shell_index]
            t_coefficient_list = [self.coefficient_list[i] for i in t_shell_index]
            assert len(t_exponent_list) == len(t_coefficient_list)

            for kk, (exponent, coefficient) in enumerate(zip(t_exponent_list, t_coefficient_list)):
                line = f"      {exponent:.8f}   {coefficient:8f}\n"
                text += line

        return text

class Det_Basis_sets(Basis_sets):
    def __init__(self,
                 nucleus_index=[],
                 shell_ang_mom=[],
                 shell_ang_mom_turbo_notation=[],
                 shell_factor=[],
                 shell_index=[],
                 exponent=[],
                 coefficient=[],
                 coefficient_imag=[],
                 prim_factor=[],
                 hyb_nucleus_index=[],
                 hyb_param_num=[],
                 hyb_shell_ang_mom=[],
                 hyb_shell_ang_mom_turbo_notation=[],
                 hyb_prim_index=[],
                 hyb_coefficient=[],
                 hyb_coefficient_imag=[],
                 number_of_additional_hybrid_orbitals=[]
                 ):
        super().__init__(nucleus_index=nucleus_index,
                         shell_ang_mom=shell_ang_mom,
                         shell_ang_mom_turbo_notation=shell_ang_mom_turbo_notation,
                         shell_factor=shell_factor,
                         shell_index=shell_index,
                         exponent=exponent,
                         coefficient=coefficient,
                         coefficient_imag=coefficient_imag,
                         prim_factor=prim_factor,
                         )
        self.hyb_nucleus_index = hyb_nucleus_index
        self.hyb_param_num = hyb_param_num
        self.hyb_shell_ang_mom = hyb_shell_ang_mom
        self.hyb_shell_ang_mom_turbo_notation = hyb_shell_ang_mom_turbo_notation
        self.hyb_prim_index = hyb_prim_index
        self.hyb_coefficient = hyb_coefficient
        self.hyb_coefficient_imag = hyb_coefficient_imag
        self.number_of_additional_hybrid_orbitals = number_of_additional_hybrid_orbitals

class Jas_Basis_sets(Basis_sets):
    def __init__(self,
                 nucleus_index=[],
                 shell_ang_mom=[],
                 shell_ang_mom_turbo_notation=[],
                 shell_factor=[],
                 shell_index=[],
                 exponent=[],
                 coefficient=[],
                 coefficient_imag=[],
                 prim_factor=[],
                 number_of_additional_hybrid_orbitals=[]
                 ):
        super().__init__(nucleus_index=nucleus_index,
                         shell_ang_mom=shell_ang_mom,
                         shell_ang_mom_turbo_notation=shell_ang_mom_turbo_notation,
                         shell_factor=shell_factor,
                         shell_index=shell_index,
                         exponent=exponent,
                         coefficient=coefficient,
                         coefficient_imag=coefficient_imag,
                         prim_factor=prim_factor,
                         )

        self.number_of_additional_hybrid_orbitals = number_of_additional_hybrid_orbitals
"""
class Hybrid_orbital:
    ...
class Atomic_orbital:
    ...
class Molecular_orbital:
    ...
"""

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    # moved to examples
    # eccep conveter -> gamess format
    basis_sets_test_dir = os.path.join(pyturbo_root, "tests", "basis_sets")
    os.chdir(basis_sets_test_dir)
    file_name = "aug-cc-pV5Z-eCEPP.dat_O"
    basis_sets = Basis_sets.parse_basis_sets_from_eCEPP_format_files(files=[file_name, file_name])
    basis_sets.contracted_to_uncontracted()
    #print(basis_set.to_text_gamess_format())
    #print(basis_set.to_text_nwchem_format(element="O"))
    #basis_sets.remove_primitive_basis(prim_index=10)
    #basis_sets.cut_exponents(thr_exp=20, nucleus_index=[0])
