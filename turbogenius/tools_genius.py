#!python
# -*- coding: utf-8 -*-

"""

Useful tools for turbogenius

Todo:
    * refactoring assert sentences. The assert should not be used for any on-the-fly check.

"""

#python modules
import os, sys
import shutil

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from pyturbo.io_fort10 import IO_fort10
from pyturbo.utils.utility import file_check, remove_file, copy_file
from pyturbo.utils.env import turbo_copyjas_command
from pyturbo.utils.execute import run

#Logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)

def copy_jastrow(fort10_to="fort.10", fort10_from="fort.10_new", args=""):
    """
        Copy Jastrow factors

        Args:
            fort10_to (str): fort.10 to which jastrow factor is copied.
            fort10_from (str): fort.10 form which jastrow factor is copied.
            args (str): option string for copyjas.x
    """
    file_check(fort10_from)
    current_dir=os.getcwd()
    copy_file(fort10_to, os.path.join(current_dir, "fort.10"))
    copy_file(fort10_from, os.path.join(current_dir, "fort.10_new"))
    run(turbo_copyjas_command + " " + args, input_name=None, output_name="out_copyjas")
    #remove_file(os.path.join(current_dir, "fort.10_new"))

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
    