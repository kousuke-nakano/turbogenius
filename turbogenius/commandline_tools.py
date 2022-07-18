#!python
# -*- coding: utf-8 -*-

"""

TurboGenius command-line tools

"""

#python modules
import os, sys
import shutil
import click

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pyturbo.io_fort10 import IO_fort10

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)

from turbo_genius_cli import cli, decorate_grpost, header, OptionEatAll

@cli.command(short_help = "visualize fort.10 structure by ASE")
@header
def view(
        operation:bool,
        log_level:str,
)->None:
    """
        Visualize a molecule or crystal structure written in fort.10

    """
    io_fort10=IO_fort10(fort10="fort.10")
    structure=io_fort10.f10structure.structure
    structure.view()

@cli.command(short_help = "convert fort.10 to a structure file by ASE-read")
@click.option("-s", "structure",
              help= 'Specify structure file name (e.g., xxxx.xsf)',
              default = "ase.xsf",
              type = str)
@header
def writestr(
        operation:bool,
        log_level:str,
        structure:str
):
    """
        Write a structure file (i.e., convert fort.10 to a structure file, e.g., xyz)

        Args:
            structure (str): structure file name. all formats supported by ASE are acceptable.

    """
    io_fort10=IO_fort10(fort10="fort.10")
    structure_=io_fort10.f10structure.structure
    structure_.write(structure)

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    from utils_workflows.env import turbo_genius_root
    os.chdir(os.path.join(turbo_genius_root, "tests", "copyjas"))

    # moved to examples
    copy_jastrow()
    