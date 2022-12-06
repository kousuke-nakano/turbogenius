#!python -u
# -*- coding: utf-8 -*-

"""

pyturbo: utilities

Todo:
    * docstrings are not completed.

"""

from __future__ import print_function

# python modules
import os, sys
import re
import shutil
import platform
import subprocess
import linecache
import numpy as np
from scipy.io import FortranFile

# python special module
from pymatgen.core.periodic_table import Element, ElementBase

# set logger
from logging import config, getLogger, StreamHandler, Formatter

logger = getLogger('pyturbo').getChild(__name__)

# turbogenius module
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from env import pyturbo_root

def get_linenum_fort12(fort12="fort.12"):
    # check column length of fort.12
    f = FortranFile(fort12, 'r')
    a = f.read_reals(dtype='float64')
    column_length = len(a)
    f.close()
    head = ("head", "<i")
    tail = ("tail", "<i")
    dt = np.dtype([head, ("a", "<{}d".format(column_length)), tail])
    fd = open(fort12, "r")
    fort12_b = np.fromfile(fd, dtype=dt, count=-1)
    data_length = len(fort12_b)
    fd.close()
    return data_length

def return_element_symbol(atomic_number):
    atomic_number = int(float(atomic_number))
    return str(ElementBase.from_Z(atomic_number))
def return_atomic_number(element):
    element = str(element)
    return Element(element).number
def remove_file(file):
    if os.path.isfile(file): os.remove(file)
def copy_file(from_file,to_file):
    try:
        shutil.copy(from_file, to_file)
    except shutil.SameFileError:
        pass
def file_check(file):
    if not os.path.exists(file):
        raise FileNotFoundError(f"{file} is not found.")
def file_check_flag(file):
    if os.path.exists(file):
        return True
    else:
        return False
def get_line_from_file(file, line_no):
    with open(file, 'r') as f:
        data = f.readlines()
    return data[line_no]
def prompt(text, checker):
    """Loop input() *FOREVER* while input is invalid."""
    while True:
        output = input(text)
        if checker(output):
            return output
def get_str_variable_type_auto(variable):
    logger.debug(f"variable={variable}")
    logger.debug(f"isdecimal={variable.isdecimal()}")

    if not variable.replace('-', '').isdecimal():
        try:
            # Convert d to e in scientific notation. python supports e only
            var = variable.replace('d', 'e').replace('D', 'e')
            result = float(var)
            # logger.debug("This is float")
            return result
        except ValueError:
            # logger.debug("This is str")
            return str(variable).strip("\'")
    else:
        # logger.debug("This is int")
        return (int(variable))
def turbo_prim_orb_type_num(orb_type_chr):
    if orb_type_chr == "s":
        return (16)
    elif orb_type_chr == "p":
        return (36)
    elif orb_type_chr == "d":
        return (37)
        # 68 -> 37 !! 12/May/2022 because both are the same in makefun.f90 (see. case(37,68)),
        # but 500 is usually associated with 37 (see. ioptorbcontr.f90)
    elif orb_type_chr == "f":
        return (48)
    elif orb_type_chr == "g":
        return (51) #51 or 88?
    elif orb_type_chr == "h":
        return (72)
    elif orb_type_chr == "i":
        return (73)
    else:
        logger.error(f"orb_type_chr={orb_type_chr} is not implemented.")
        raise NotImplementedError(f"Not implemented orb_type_chr={orb_type_chr}")
def turbo_cont_orb_type_num(orb_type_chr):
    if orb_type_chr == "s":
        return (300)
    elif orb_type_chr == "p":
        return (400)
    elif orb_type_chr == "d":
        return (500)
    elif orb_type_chr == "f":
        return (600)
    elif orb_type_chr == "g":
        return (700)
    elif orb_type_chr == "h":
        return (800)
    elif orb_type_chr == "i":
        return (900)
    else:
        logger.error(f"orb_type_chr={orb_type_chr} is not implemented.")
        raise NotImplementedError(f"Not implemented orb_type_chr={orb_type_chr}")
def return_ang_mom(orb_typ_chr):
    if orb_typ_chr in {"s"}:
        return (0)
    elif orb_typ_chr in {"p"}:
        return (1)
    elif orb_typ_chr in {"d"}:
        return (2)
    elif orb_typ_chr in {"f"}:
        return (3)
    elif orb_typ_chr in {"g"}:
        return (4)
    elif orb_typ_chr in {"h"}:
        return (5)
    elif orb_typ_chr in {"i"}:
        return (6)
    else:
        logger.error(f"orb_type={orb_typ_chr} is not implemented.")
        raise NotImplementedEroor(f"orb_type={orb_typ_chr} is not implemented.")
def return_orbchr(ang_mom):
    if ang_mom == 0:
        return ("s")
    elif ang_mom == 1:
        return ("p")
    elif ang_mom == 2:
        return ("d")
    elif ang_mom == 3:
        return ("f")
    elif ang_mom == 4:
        return ("g")
    elif ang_mom == 5:
        return ("h")
    elif ang_mom == 6:
        return ("i")
    else:
        raise NotImplementedError(f"angmom={ang_mom} is not implemented.")
def return_orb_type_chr(num_orb_type):
    if num_orb_type in {16, 300}:
        return ("s")
    elif num_orb_type in {36, 400}:
        return ("p")
    elif num_orb_type in {37, 68, 500}:
        return ("d")
    elif num_orb_type in {48, 600}:
        return ("f")
    elif num_orb_type in {51, 700}:
        return ("g")
    elif num_orb_type in {72, 800}:
        return ("h")
    elif num_orb_type in {73, 900}:
        return ("i")
    elif num_orb_type in {900000}:
        return ("hyb")
    elif num_orb_type in {1000000}:
        return ("mol")
    elif num_orb_type in {200}:
        return ("jas_const")
    else:
        raise NotImplementedError(f"Not suppported orb_type={num_orb_type}")
def return_contraction_flag(orb_type_num):
    if orb_type_num in {16, 100, 36, 37, 68, 48, 51, 72, 73}:
        contraction = False
    elif orb_type_num in {300, 400, 500, 600, 700, 800, 900}:
        contraction = True
    else:
        raise NotImplementedError(f"Not suported orb_type_num={orb_type_num}")
    return contraction
def return_num_twobody_and_flag_onebody(jastrow_type):
    if jastrow_type in {-15, -22, -26}:
        num_twobody = 1
        flag_onebody = True
    elif jastrow_type in {-27}:
        num_twobody = 2
        flag_onebody = True
    elif jastrow_type in {-5, -6}:
        num_twobody = 1
        flag_onebody = False
    elif jastrow_type in {0}:
        num_twobody = 0
        flag_onebody = False
    else:
        raise NotImplementedError(f"Jastrow = {jastrow_type} is not implemented.")

    return num_twobody, flag_onebody

# workarounds!! They work well but not pythonic. to be reafactored?
def pygrep_lineno(file, keyword):
    cmd=f"grep \'{keyword}\' -n {file} | cut -d ':' -f 1"
    sys_env = os.environ.copy()
    try:
        lineno=int(subprocess.check_output(cmd, shell=True, env=sys_env)) - 1
    except ValueError:
        cmd=f"grep -c '' {file}"
        lineno=int(subprocess.check_output(cmd, shell=True, env=sys_env))

    return lineno
def pysed_replace(file, value, lineno, index, inplace=False):
    if platform.system() == "Darwin":
        if shutil.which('gsed') is None:
            logger.error("The BSD sed on MacOS is not supported.")
            logger.error("Pls. install gsed via homebrew, i.e., brew install gnu-sed")
            raise NotImplementedError
        else:
            sed = "gsed"
            awk = "awk"
    else:
        sed = "sed"
        awk = "awk"

    sys_env = os.environ.copy()
    cmds=[]

    lineno+=1

    if re.match(r".*bash.*", sys_env["SHELL"]) or re.match(r".*zsh.*", sys_env["SHELL"]):
        cmds = [
            f"line=`{sed} -n {lineno}p {file}`",
            f"mod_line=`echo $line |  {awk} '{{FS=\" \";OFS=\" \"}}{{${index+1}={value}}}1'`",
        ]
    elif re.match(r".*csh.*", sys_env["SHELL"]) or re.match(r".*tsch.*", sys_env["SHELL"]):
        cmds = [
            f"line=`{sed} -n {lineno}p {file}`",
            f"mod_line=`echo $line |  {awk} '{{FS=\" \";OFS=\" \"}}{{${index+1}={value}}}1'`",
        ]

    if inplace:
        cmds += [f"cp {file} {file}_bak"]

    """
    # test case
    cmds += [
        f"{sed} \"{lineno}d\" {file} > {file}_", # test
        f"{sed} \"{lineno}i \\ $mod_line\" {file}_ > {file}_io", # test
        f"rm {file}_"
    ]

    """
    #"""
    # without IOs
    cmds += [
        f"{sed} -i \"{lineno}d\" {file}",
        f"{sed} -i \"{lineno}i \\ $mod_line\" {file}"
    ]
    #"""

    cmd="; ".join(cmds)
    logger.debug(cmd)
    subprocess.check_call(cmd, shell=True, env=sys_env)

def pysed_replace_lines(file, lineno_list, value_list, index_list, inplace=False, cmd_chunk_num=100):

    if platform.system() == "Darwin":
        if shutil.which('gsed') is None:
            logger.error("The BSD sed on MacOS is not supported.")
            logger.error("Pls. install gsed via homebrew, i.e., brew install gnu-sed")
            raise NotImplementedError
        else:
            sed = "gsed"
            awk = "awk"
    else:
        sed = "sed"
        awk = "awk"

    sys_env = os.environ.copy()

    assert len(lineno_list) == len(value_list)
    assert len(value_list) == len(index_list)

    if inplace:
        cmd = f"cp {file} {file}_bak"
        subprocess.check_call(cmd, shell=True, env=sys_env)

    counter = 0
    cmds = []
    for value_list_l, index_list_l, lineno in zip(value_list, index_list, lineno_list):
        lineno += 1 # because, grep starts from 1, but python starts from 0
        if re.match(r".*bash.*", sys_env["SHELL"]) or re.match(r".*zsh.*", sys_env["SHELL"]):
            cmds += [
                f"line=`{sed} -n {lineno}p {file}`"
            ]
            for value, index in zip(value_list_l, index_list_l):
                cmds += [
                    f"line=`echo $line |  {awk} '{{FS=\" \";OFS=\" \"}}{{${index+1}={value}}}1'`"
                ]
                # index+1 # because, grep starts from 1, but python starts from 0
        elif re.match(r".*csh.*", sys_env["SHELL"]) or re.match(r".*tsch.*", sys_env["SHELL"]):
            cmds += [
                f"line=`{sed} -n {lineno}p {file}`"
            ]
            for value, index in zip(value_list_l, index_list_l):
                cmds += [
                    f"line=`echo $line |  {awk} '{{FS=\" \";OFS=\" \"}}{{${index+1}={value}}}1'`"
                ]
                # index+1 # because, grep starts from 1, but python starts from 0

        # without IOs
        cmds += [
            f"{sed} -i \"{lineno}d\" {file}",
            f"{sed} -i \"{lineno}i \\ $line\" {file}"
        ]

        counter += 1
        if counter >= cmd_chunk_num:
            cmd="; ".join(cmds)
            logger.debug(cmd)
            subprocess.check_call(cmd, shell=True, env=sys_env)
            counter=0
            cmds=[]

    if len(cmds) !=0:
        cmd = "; ".join(cmds)
        logger.debug(cmd)
        subprocess.check_call(cmd, shell=True, env=sys_env)
        counter = 0
        cmds = []

def pygetline(filename, lineno, clearcache=True): # clearchache should be true!! as a default. # reasons for bugs.
    #logger.debug("get line!")
    line=linecache.getline(filename=filename, lineno=lineno+1)
    if clearcache:
        linecache.clearcache()
    return line

def remove_new_parameter_lines_in_fort10():

    output_buffer = [
            '#!/bin/bash',
            "line_num_key=`grep 'new parameters' -n fort.10 | cut -d ':' -f 1 | head -n 1`",
            'if [ -n "$line_num_key" ]; then',
            "line_num=`expr ${line_num_key} - 1`",
            "head -n ${line_num} fort.10 > fort.10_",
            "mv fort.10_ fort.10",
            "fi"
        ]
    output_buffer = '\n'.join(output_buffer)
    with open("./run_local.sh", "w") as f:
        f.writelines(output_buffer)

    sys_env = os.environ.copy()
    cmds=["chmod +x ./run_local.sh", "./run_local.sh", "./run_local.sh"]
    for cmd in cmds:
        subprocess.check_call(cmd, shell=True, env=sys_env)

if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    os.chdir(os.path.join(pyturbo_root, "tests", "fort10"))
    #pysed_replace(file="fort.10", value=5, lineno=2, index=1, inplace=True)
    #value=Value(value=15, lineno=1, index=0, file="fort.10")
    #print(value.v)
    #value.replace(38)
    #print(pygrep_lineno(file="fort.10", keyword="Ion coordinates"))
