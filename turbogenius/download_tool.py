#!python
# -*- coding: utf-8 -*-

"""
Standalone tool for downloading basis sets and pseudo potentials.

This tool provides a command-line interface for downloading basis sets and
pseudo potentials from various databases (BSE, ccECP, BFD).

Usage:
    python -m turbogenius.download_tool -db BSE
    python -m turbogenius.download_tool -db ccECP --force
    python -m turbogenius.download_tool -db BFD -s 2.0
"""

# python modules
import os
import sys
import argparse
from typing import Optional

# Logger
from logging import getLogger, StreamHandler, Formatter

# turbogenius modules
from turbogenius.pyturbo.utils.downloader import ccECP, BSE, BFD
from turbogenius.utils_workflows.env import turbo_genius_tmp_dir

logger = getLogger("Turbo-Genius").getChild(__name__)


def download_database(
    database: str = "",
    sleep_time: float = 1.5,
    force: bool = False,
    output_dir: Optional[str] = None,
) -> None:
    """
    Downloading basis set and pseudo potential database from the Internet.

    Args:
        database (str): name of database, it should be chosen from database_list=["BFD", "ccECP", "BSE"]
        sleep_time (float): sleeping time for downloading (float)
        force (bool): if true, overwrite an existing database
        output_dir (str, optional): output directory for downloaded files. 
                                     If None, uses default location (~/.turbo_genius_tmp/)

    """
    # Use specified output directory or default
    if output_dir is None:
        base_output_dir = turbo_genius_tmp_dir
    else:
        base_output_dir = os.path.abspath(output_dir)
    
    basis_sets_output_dir = os.path.join(
        base_output_dir, "basis_set", database
    )
    pseudo_potential_output_dir = os.path.join(
        base_output_dir, "pseudo_potential", database
    )
    logger.info(f"Database: {database}")
    logger.info(f"Basis sets output directory: {basis_sets_output_dir}")
    logger.info(f"Pseudo potential output directory: {pseudo_potential_output_dir}")

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
        raise NotImplementedError(
            f"Database '{database}' is not supported. Choose from: BFD, ccECP, BSE"
        )

    if force or not database_is_exist:
        os.makedirs(basis_sets_output_dir, exist_ok=True)
        os.makedirs(pseudo_potential_output_dir, exist_ok=True)
        logger.info("Turbo-Genius database has not been downloaded yet.")
        logger.info("Downloading all the data from the web.")
        logger.info(
            f"Basis sets and PPs will be downloaded to {base_output_dir}"
        )

        # Download data
        if database == "ccECP":
            logger.info("Cloning ccECP repository from GitHub...")
        elif database == "BFD":
            logger.info("Cloning BFD-ECP repository from GitHub...")
        elif database == "BSE":
            logger.info("Downloading basis sets from Basis Set Exchange...")
        loader.all_to_file(sleep_time=sleep_time)

        with open(os.path.join(basis_sets_output_dir, "completed"), "w") as f:
            f.write("completed")
        if database != "BSE":  # BSE doesn't have pseudo potentials
            with open(
                os.path.join(pseudo_potential_output_dir, "completed"), "w"
            ) as f:
                f.write("completed")
        logger.info("Download completed successfully!")

    else:
        logger.info("You have already downloaded the database")
        logger.info(
            "If you want to download the database again, switch on the force option."
        )


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Download basis sets and pseudo potentials from various databases.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -db BSE
  %(prog)s -db ccECP --force
  %(prog)s -db BFD -s 2.0
  %(prog)s -db BSE -o /path/to/custom/directory

Available databases:
  BSE:    All-electron basis sets from Basis Set Exchange
  ccECP:  ccECP basis sets and pseudo potentials
  BFD:    BFD basis sets and pseudo potentials

Downloaded files are stored in ~/.turbo_genius_tmp/ by default.
Use -o/--output-dir to specify a custom directory.
        """,
    )
    parser.add_argument(
        "-db",
        "--database",
        dest="database",
        required=True,
        choices=["BFD", "ccECP", "BSE"],
        help="Database to download: BFD, ccECP, or BSE",
    )
    parser.add_argument(
        "-f",
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force download even if database already exists",
    )
    parser.add_argument(
        "-s",
        "--sleep-time",
        dest="sleep_time",
        type=float,
        default=1.5,
        help="Sleep time between downloads (seconds, default: 1.5)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Enable verbose logging",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=str,
        default=None,
        help="Output directory for downloaded files (default: ~/.turbo_genius_tmp/)",
    )

    args = parser.parse_args()

    # Set up logger
    log_level = "DEBUG" if args.verbose else "INFO"
    logger_t = getLogger("Turbo-Genius")
    logger_t.setLevel(log_level)
    stream_handler = StreamHandler()
    stream_handler.setLevel(log_level)
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger_t.addHandler(stream_handler)

    # Download database
    try:
        download_database(
            database=args.database,
            sleep_time=args.sleep_time,
            force=args.force,
            output_dir=args.output_dir,
        )
    except Exception as e:
        logger.error(f"Error during download: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

