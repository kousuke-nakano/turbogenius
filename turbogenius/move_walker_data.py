#!/usr/bin/env python3
# coding: utf-8

"""

A utility to move walker data (positions and local energies)
to TREXIO file.

"""

try:
    from turbogenius._version import version as turbogenius_version
except (ModuleNotFoundError, ImportError):
    turbogenius_version = "unknown"

from logging import getLogger, StreamHandler, Formatter

def main():

    import argparse

    parser = argparse.ArgumentParser(
        description="This utility move walker data from fort.12(.new) to TREXIO file."
    )
    parser.add_argument("trexio_file", help="Name of TREXIO file")
    parser.add_argument(
        "-f",
        "--fort-12",
        help="Name of the fort.12 file",
        type=str,
        default="fort.12",
    )
    parser.add_argument(
        "-n",
        "--fort-12-new",
        help="Name of the fort.12.new file",
        type=str,
        default="fort.12.new",
    )
    parser.add_argument(
        "-log",
        "--loglevel",
        help="logger setlevel",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
    )
    args = parser.parse_args()

    logger = getLogger("Turbo-Genius").getChild(__name__)
    logger.setLevel(args.loglevel)
    stream_handler = StreamHandler()
    stream_handler.setLevel(args.loglevel)
    if args.loglevel in {"DEBUG"}:
        handler_format = Formatter(
            "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
        )
    else:
        handler_format = Formatter("%(message)s")
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    logger.info(f"turbogenius {turbogenius_version}")

    def gen(f):
        def fun(*args, **kwargs):
            try:
                while True:
                    yield f(*args, **kwargs)
            except scipy.io.FortranEOFError:
                return
        return fun

    import os
    import numpy as np
    import scipy

    try:
        import trexio
    except:
        logger.info(f"TREXIO library is missing")
        os.exit(1)

    with trexio.File(args.trexio_file, mode='r', back_end=trexio.TREXIO_HDF5) as trexio_file:
        nel = trexio.read_electron_num(trexio_file)

    walker_positions = []
    E_loc = []

    with scipy.io.FortranFile(args.fort_12) as f12:
        for ii, d in enumerate(gen(f12.read_reals)(dtype=np.float64)):
            E_loc.append(d[2])

    logger.info("Found {len(E_loc)} local energies")

    with scipy.io.FortranFile(args.fort_12_new) as f12:
        for ii, d in enumerate(gen(f12.read_reals)(dtype=np.float32)):
            walker_positions.append(d[2:3*nel+2])

    E_loc = np.array(E_loc, dtype=np.float64)
    walker_positions = np.array(walker_positions, dtype=np.float64)
    walker_positions = walker_positions.reshape(len(walker_positions), -1, 3)

    with trexio.File(args.trexio_file, mode='w', back_end=trexio.TREXIO_HDF5) as trexio_file:
        trexio.write_qmc_num(trexio_file, walker_positions.shape[0])
        trexio.write_qmc_e_loc(trexio_file, E_loc)
        trexio.write_qmc_point(trexio_file, walker_positions)

    logger.info("Done")

if __name__ == "__main__":
    logger = getLogger("Turbo-Genius")
    logger.setLevel("INFO")
    stream_handler = StreamHandler()
    stream_handler.setLevel("INFO")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    main()
