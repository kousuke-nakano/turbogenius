#!python
# -*- coding: utf-8 -*-


"""

Database setup related classes and methods


"""

# python modules
import os

# Logger
from logging import getLogger, StreamHandler, Formatter

# turbogenius modules
from turbogenius.pyturbo.utils.downloader import ccECP, BSE, BFD
from turbogenius.utils_workflows.env import (
    turbo_genius_root,
    turbo_genius_tmp_dir,
)

logger = getLogger("Turbo-Genius").getChild(__name__)

all_electron_basis_set_list = BSE.list_of_basis_all
ccECP_basis_set_list = ccECP.list_of_basis_all
BFD_basis_set_list = BFD.list_of_basis_all
ecp_list = ["BFD", "ccECP"]


def database_setup(
    database: str = "", sleep_time: float = 1.5, force: bool = False
) -> None:
    """
    Downloading basis set and pseudo potential database from the Internet.

    Args:
        database (str): name of database, it should be chosen from database_list=["BFD", "ccECP", "BSE"]
        sleep_time (float): sleeping time for downloading (float)
        force (bool): if true, overwrite an existing database

    """

    basis_sets_output_dir = os.path.join(
        turbo_genius_tmp_dir, "basis_set", database
    )
    pseudo_potential_output_dir = os.path.join(
        turbo_genius_tmp_dir, "pseudo_potential", database
    )
    logger.debug(basis_sets_output_dir)
    logger.debug(pseudo_potential_output_dir)
    if database == "BFD":  # pseudo potential
        loader = BFD(
            basis_sets_output_dir=basis_sets_output_dir,
            pseudo_potential_output_dir=pseudo_potential_output_dir,
        )
        if os.path.isfile(
            os.path.join(basis_sets_output_dir, "completed")
        ) and os.path.isfile(
            os.path.join(pseudo_potential_output_dir, "completed")
        ):
            database_is_exist = True
        else:
            database_is_exist = False
    elif database == "ccECP":  # pseudo potential
        loader = ccECP(
            basis_sets_output_dir=basis_sets_output_dir,
            pseudo_potential_output_dir=pseudo_potential_output_dir,
        )
        if os.path.isfile(
            os.path.join(basis_sets_output_dir, "completed")
        ) and os.path.isfile(
            os.path.join(pseudo_potential_output_dir, "completed")
        ):
            database_is_exist = True
        else:
            database_is_exist = False
    elif database == "BSE":  # all-electron
        loader = BSE(
            basis_sets_output_dir=basis_sets_output_dir,
            pseudo_potential_output_dir=pseudo_potential_output_dir,
        )
        if os.path.isfile(os.path.join(basis_sets_output_dir, "completed")):
            database_is_exist = True
        else:
            database_is_exist = False
    else:
        logger.error(f"database = {database} is not implemented.")
        raise NotImplementedError

    if force or not database_is_exist:
        os.makedirs(basis_sets_output_dir, exist_ok=True)
        os.makedirs(pseudo_potential_output_dir, exist_ok=True)
        logger.info("Turbo-Genius database has not been downloaded yet.")
        logger.info("Downloading all the data from the web.")
        logger.info(
            f"Basis sets and PPs will be downloaded to {turbo_genius_tmp_dir}"
        )
        loader.all_to_file(sleep_time=sleep_time)
        with open(os.path.join(basis_sets_output_dir, "completed"), "w") as f:
            f.write("completed")
        with open(
            os.path.join(pseudo_potential_output_dir, "completed"), "w"
        ) as f:
            f.write("completed")

    else:
        logger.info("You have already downloaded the database")
        logger.info(
            "If you want to download the database again, switch on the force option."
        )


if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("DEBUG")
    logger_p = getLogger("pyturbo")
    logger_p.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)
    logger_p.addHandler(stream_handler)

    os.chdir(os.path.join(turbo_genius_root, "tests", "vmc"))

    database_setup(sleep_time=2, database="ccECP", force=False)
