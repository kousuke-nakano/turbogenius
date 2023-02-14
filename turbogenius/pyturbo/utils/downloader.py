#!python
# -*- coding: utf-8 -*-

"""

pyturbo: downloader, for downloading basis set and pseudo potential information from www.

Todo:
    * docstrings are not completed.

"""

# python modules
import os
import re
from urllib.request import urlopen
from ase.data import chemical_symbols
import git
import tempfile
import pathlib
import itertools
import time
from tqdm import tqdm
from tqdm.contrib import itertools as titertools

# pymatgen
from pymatgen.core.periodic_table import Element

# basissetexchange
import basis_set_exchange as bse

# turbo-genius modules
from turbogenius.pyturbo.utils.env import pyturbo_root

# Logger
from logging import getLogger, StreamHandler, Formatter

logger = getLogger("pyturbo").getChild(__name__)


class ccECP:
    C_URL = "https://github.com/QMCPACK/pseudopotentiallibrary.git"
    list_of_basis_all = [f"cc-pV{s}Z" for s in "DTQ56"]
    list_of_basis_all += [f"ang-{x}" for x in list_of_basis_all]

    def __init__(
        self, basis_sets_output_dir=None, pseudo_potential_output_dir=None
    ):
        if basis_sets_output_dir is None:
            self.basis_sets_output_dir = None
        else:
            self.basis_sets_output_dir = basis_sets_output_dir
        if pseudo_potential_output_dir is None:
            self.ecp_output_dir = None
        else:
            self.ecp_output_dir = pseudo_potential_output_dir

    def to_file(self, element_list=None, basis_list=None):
        if self.basis_sets_output_dir is not None:
            os.makedirs(self.basis_sets_output_dir, exist_ok=True)
        if self.ecp_output_dir is not None:
            os.makedirs(self.ecp_output_dir, exist_ok=True)
        elements = {}

        def add_data(p, tempdir, elements=elements):
            if self.basis_sets_output_dir is None:
                return
            element = str(p.parent.parent.name)
            typ = str(p.parent.name)
            typ_core = typ.replace("ccECP", "")
            logger.info(typ)
            logger.info(typ_core)
            # element_path = p.parent.parent

            typ = str(p.parent.name)
            # typ_path = str(p.parent.name)

            if element_list is not None:
                if element not in element_list:
                    return

            # Load Basis sets
            # basis = []
            for r in (p.parent).glob("**/*"):
                if re.match("[A-z]{1,2}\.[A-z\-]*cc-.*\.gamess", r.name):
                    if re.match("[A-z]{1,2}\.[A-z\-]*cc-.*\.gamess", r.name):
                        basis_name = re.match(
                            "[A-z]{1,2}\.([A-z\-]*cc-.*)\.gamess", r.name
                        ).group(1)
                        if basis_list is not None:
                            if basis_name not in basis_list:
                                continue
                        logger.info(element)
                        name = f"{basis_name}{typ_core}.basis"
                        with open(r, "r") as fhandle:
                            with open(
                                os.path.join(
                                    self.basis_sets_output_dir,
                                    f"{element}_{name}",
                                ),
                                "w",
                            ) as fhandle_out:
                                fhandle_out.write(fhandle.read())

        tempdir = pathlib.Path(tempfile.mkdtemp())
        git.Repo.clone_from(self.C_URL, tempdir)
        # logger.debug(list((tempdir/"recipes").glob("**/*")))
        for p in (tempdir / "recipes").glob("**/*"):
            if re.match("[A-z]{1,2}\.ccECP\.gamess", p.name):
                logger.debug(p)
                add_data(p, tempdir)

                if self.ecp_output_dir is None:
                    return

                element = str(p.parent.parent.name)
                typ = str(p.parent.name)
                typ_core = typ.replace("ccECP", "")
                logger.info(typ)
                logger.info(typ_core)

                if element_list is not None:
                    if element not in element_list:
                        continue
                logger.info(p.name)
                logger.info(p.name.replace(".ccECP.gamess", ""))
                with open(p, "r") as fhandle:
                    with open(
                        os.path.join(
                            self.ecp_output_dir,
                            f"{p.name.replace('.ccECP.gamess', '_ccECP'+typ_core+'.pseudo')}",
                        ),
                        "w",
                    ) as fhandle_out:
                        fhandle_out.write(fhandle.read())

    def all_to_file(self, sleep_time=1):
        self.to_file(
            element_list=chemical_symbols, basis_list=self.list_of_basis_all
        )


class BSE:
    # URL: https://www.basissetexchange.org
    list_of_basis_all = [f"cc-pV{s}Z" for s in "DTQ56"]
    list_of_basis_all += [f"ang-{x}" for x in list_of_basis_all]

    def __init__(
        self, basis_sets_output_dir=None, pseudo_potential_output_dir=None
    ):
        if basis_sets_output_dir is None:
            self.basis_sets_output_dir = None
        else:
            self.basis_sets_output_dir = basis_sets_output_dir

    def to_file(self, element_list, basis_list, sleep_time=1.0):
        if self.basis_sets_output_dir is None:
            return
        os.makedirs(self.basis_sets_output_dir, exist_ok=True)
        for e, b in itertools.product(element_list, basis_list):
            # time.sleep(sleep_time + random.randint(sleep_time))
            try:
                basis_text = bse.get_basis(b, elements=e, fmt="gamess_us")
                E = Element(e)
                pat = re.compile(
                    "^.*?DATA.*" + E.long_name + "\n(.+)\$END",
                    re.DOTALL | re.IGNORECASE,
                )
                m = pat.match(basis_text)
                if m:
                    bas = m.group(1)
                    with open(
                        os.path.join(
                            self.basis_sets_output_dir, f"{e}_{b}.basis"
                        ),
                        "w",
                    ) as fhandle:
                        fhandle.write(bas)
            except KeyError:
                logger.debug(
                    f"element={e}, basis={b} do not exist in the database."
                )

    def all_to_file(self, sleep_time=1):
        self.to_file(
            element_list=chemical_symbols,
            basis_list=self.list_of_basis_all,
            sleep_time=sleep_time,
        )


class BFD:
    U = "http://burkatzki.com/pseudos/step4.2.php?format=gamess&element={e}&basis={b}"
    list_of_basis_all = [f"v{s}z" for s in "dtq56"]
    # list_of_basis_all += [f"{x}_ano" for x in list_of_basis_all]

    def __init__(
        self, basis_sets_output_dir=None, pseudo_potential_output_dir=None
    ):
        if basis_sets_output_dir is None:
            self.basis_sets_output_dir = None
        else:
            self.basis_sets_output_dir = basis_sets_output_dir
        if pseudo_potential_output_dir is None:
            self.ecp_output_dir = None
        else:
            self.ecp_output_dir = pseudo_potential_output_dir

    def to_file(self, element_list, basis_list, sleep_time=1.5):
        if self.basis_sets_output_dir is None and self.ecp_output_dir is None:
            return
        if self.basis_sets_output_dir is not None:
            os.makedirs(self.basis_sets_output_dir, exist_ok=True)
        if self.ecp_output_dir is not None:
            os.makedirs(self.ecp_output_dir, exist_ok=True)
        if "X" in element_list:
            element_list.remove("X")
        with tqdm(titertools.product(element_list, basis_list)) as pbar:
            for i, (e, b) in enumerate(pbar):
                logger.debug(f"b={b:s}")
                logger.debug(f"[BFD] e={e:s} b={b:s}")

                # check already downloaded or not.
                if self.basis_sets_output_dir is not None:
                    if self.ecp_output_dir is not None:
                        flag_downloaded = (
                            os.path.isfile(
                                os.path.join(
                                    self.ecp_output_dir, f"{e}_BFD.pseudo"
                                )
                            )
                            or os.path.isfile(
                                os.path.join(
                                    self.ecp_output_dir, f"{e}_BFD.NaN"
                                )
                            )
                        ) and (
                            os.path.isfile(
                                os.path.join(
                                    self.basis_sets_output_dir,
                                    f"{e}_{b}.basis",
                                )
                            )
                            or os.path.isfile(
                                os.path.join(
                                    self.basis_sets_output_dir, f"{e}_{b}.NaN"
                                )
                            )
                        )
                    else:
                        flag_downloaded = os.path.isfile(
                            os.path.join(
                                self.basis_sets_output_dir, f"{e}_{b}.basis"
                            )
                        ) or os.path.isfile(
                            os.path.join(
                                self.basis_sets_output_dir, f"{e}_{b}.NaN"
                            )
                        )
                else:
                    if self.ecp_output_dir is not None:
                        flag_downloaded = os.path.isfile(
                            os.path.join(
                                self.ecp_output_dir, f"{e}_BFD.pseudo"
                            )
                        ) or os.path.isfile(
                            os.path.join(self.ecp_output_dir, f"{e}_BFD.NaN")
                        )

                    else:
                        raise ValueError

                if not flag_downloaded:
                    logger.debug(f"b={b:s}")
                    pbar.set_description(f"[BFD] e={e:s} b={b:s}")
                    time.sleep(sleep_time)
                    logger.debug(self.U.format(b=b, e=e))
                    source = urlopen(self.U.format(b=b, e=e)).read()
                    source = str(source)
                    logger.debug(source)
                    pat = re.compile("^.*?(" + e + "\s0.*$)", re.M | re.DOTALL)
                    x = pat.sub("\g<1>", source)
                    x = re.sub("\<br\s*/\>", "\n", x)
                    x = re.sub("\&nbsp", "", x)
                    x = re.sub("\&nbsp", "", x)
                    x = re.sub(".*html.*$", "", x)
                    pat = re.compile(
                        "^.*?(" + e + "\s\d.*)(" + e + ".*QMC.*)$",
                        re.M | re.DOTALL,
                    )
                    m = pat.match(x)
                    if m:
                        bas = m.group(1)
                        bas = "\n".join(
                            [i for i in bas.split("\n")[1:]]
                        ).rstrip(
                            "\n"
                        )  # remove the first line.
                        ecp = m.group(2)  # ECP is OK.
                        if self.ecp_output_dir is not None:
                            with open(
                                os.path.join(
                                    self.ecp_output_dir, f"{e}_BFD.pseudo"
                                ),
                                "w",
                            ) as fo:
                                # logger.info(f"ecp -> {os.path.join(self.ecp_output_dir, f'{e}_BFD.pseudo')}")
                                fo.write(ecp)
                        if self.basis_sets_output_dir is not None:
                            if len(bas) > 15:
                                with open(
                                    os.path.join(
                                        self.basis_sets_output_dir,
                                        f"{e}_{b}.basis",
                                    ),
                                    "w",
                                ) as fo:
                                    # logger.info(f"basis -> {os.path.join(self.basis_sets_output_dir, f'{e}_{b}.basis')}")
                                    fo.write(bas)
                            else:
                                with open(
                                    os.path.join(
                                        self.basis_sets_output_dir,
                                        f"{e}_{b}.NaN",
                                    ),
                                    "w",
                                ) as fo:
                                    fo.write("NA in database")
                    else:
                        with open(
                            os.path.join(self.ecp_output_dir, f"{e}_BFD.NaN"),
                            "w",
                        ) as fo:
                            fo.write("NA in database")
                        with open(
                            os.path.join(
                                self.basis_sets_output_dir, f"{e}_{b}.NaN"
                            ),
                            "w",
                        ) as fo:
                            fo.write("NA in database")
                else:
                    logger.info(
                        f"e={e:s} b={b:s} already downloaded or NaN in the database"
                    )

    def all_to_file(self, sleep_time=1.5):
        logger.debug(self.list_of_basis_all)
        self.to_file(
            basis_list=self.list_of_basis_all,
            element_list=chemical_symbols,
            sleep_time=sleep_time,
        )


if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter(
        "%(name)s - %(levelname)s - %(lineno)d - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    downloader_test_dir = os.path.join(pyturbo_root, "tests", "downloader")
    os.chdir(downloader_test_dir)

    # moved to examples
