#!python
# -*- coding: utf-8 -*-

"""

Useful tools for turbogenius

Todo:
    * refactoring assert sentences.
    * The assert should not be used for any on-the-fly check.

"""

# python modules
import os
import shutil

# Logger
from logging import getLogger, StreamHandler, Formatter

# turboworkflows modules
from turbogenius.pyturbo.utils.utility import (
    file_check,
    copy_file,
)
from turbogenius.pyturbo.utils.env import turbo_copyjas_command
from turbogenius.pyturbo.utils.execute import run

logger = getLogger("Turbo-Genius").getChild(__name__)


def copy_jastrow(
    fort10_to: str = "fort.10",
    fort10_from: str = "fort.10_new",
    twist_flag: bool = False,
):
    """
    Copy Jastrow factors

    Args:
        fort10_to (str): fort.10 to which jastrow factor is copied.
        fort10_from (str): fort.10 form which jastrow factor is copied.
        twist_flag (bool): twist average or not.
    """
    current_dir = os.getcwd()
    file_check(fort10_from)
    file_check(fort10_to)
    copy_file(fort10_to, os.path.join(current_dir, "fort.10"))
    copy_file(fort10_from, os.path.join(current_dir, "fort.10_new"))

    run(
        turbo_copyjas_command + " ", input_name=None, output_name="out_copyjas"
    )

    if twist_flag:
        copy_jastrow_twist()


def copy_jastrow_twist():
    """
    Copy Jastrow factors with twist in turborvb.scratch
    """

    logger.info(
        "Additional commands are needed for averaging Jas. mat. with k average"
    )
    logger.info("cp fort.10 turborvb.scratch/fort.10;")
    logger.info("cp fort.10 turborvb.scratch/fort.10_new;")
    logger.info("cp kp_info.dat turborvb.scratch/kp_info.dat;")
    logger.info("cp parminimized.d turborvb.scratch/parminimized.d;")
    logger.info("cd turborvb.scratch/")
    logger.info("copyjas.x kpoints")

    # file check
    file_check("kp_info.dat")
    file_check("parminimized.d")
    file_check("turborvb.scratch")

    # copyjastrows with twist
    shutil.copy("fort.10", os.path.join("turborvb.scratch", "fort.10"))
    shutil.copy("fort.10", os.path.join("turborvb.scratch", "fort.10_new"))
    shutil.copy("kp_info.dat", os.path.join("turborvb.scratch", "kp_info.dat"))
    shutil.copy(
        "parminimized.d", os.path.join("turborvb.scratch", "parminimized.d")
    )

    # run copyjas
    current_dir = os.getcwd()
    os.chdir("turborvb.scratch")
    run(
        turbo_copyjas_command + " " + "kpoints",
        input_name=None,
        output_name="out_copyjas",
    )
    os.chdir(current_dir)


if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    from utils_workflows.env import turbo_genius_root

    os.chdir(os.path.join(turbo_genius_root, "tests", "copyjas"))

    # moved to examples
    copy_jastrow()
