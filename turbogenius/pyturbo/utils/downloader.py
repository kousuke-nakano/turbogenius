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
import git
import tempfile
import pathlib
import itertools
from typing import Optional

# pymatgen
from pymatgen.core.periodic_table import Element

# ASE
from ase.data import chemical_symbols

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
        self,
        basis_sets_output_dir: Optional[str] = None,
        pseudo_potential_output_dir: Optional[str] = None,
    ):
        if basis_sets_output_dir is None:
            self.basis_sets_output_dir = None
        else:
            self.basis_sets_output_dir = basis_sets_output_dir
        if pseudo_potential_output_dir is None:
            self.ecp_output_dir = None
        else:
            self.ecp_output_dir = pseudo_potential_output_dir

    def to_file(
        self,
        element_list: Optional[list] = None,
        basis_list: Optional[list] = None,
    ):
        if element_list is None:
            element_list = []
        if basis_list is None:
            basis_list = []
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
            logger.debug(typ)
            logger.debug(typ_core)
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
                        logger.debug(element)
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
                logger.debug(typ)
                logger.debug(typ_core)

                if element_list is not None:
                    if element not in element_list:
                        continue
                logger.debug(p.name)
                logger.debug(p.name.replace(".ccECP.gamess", ""))
                with open(p, "r") as fhandle:
                    with open(
                        os.path.join(
                            self.ecp_output_dir,
                            f"{p.name.replace('.ccECP.gamess', '_ccECP'+typ_core+'.pseudo')}",
                        ),
                        "w",
                    ) as fhandle_out:
                        fhandle_out.write(fhandle.read())

    def all_to_file(self, sleep_time: float = 1):
        self.to_file(element_list=chemical_symbols, basis_list=self.list_of_basis_all)


class BSE:
    # URL: https://www.basissetexchange.org
    list_of_basis_all = [f"cc-pV{s}Z" for s in "DTQ56"]
    list_of_basis_all += [f"ang-{x}" for x in list_of_basis_all]

    def __init__(
        self,
        basis_sets_output_dir: Optional[str] = None,
        pseudo_potential_output_dir: Optional[str] = None,
    ):
        if basis_sets_output_dir is None:
            self.basis_sets_output_dir = None
        else:
            self.basis_sets_output_dir = basis_sets_output_dir

    def to_file(
        self,
        element_list: Optional[list] = None,
        basis_list: Optional[list] = None,
        sleep_time: int = 1.0,
    ):
        if element_list is None:
            element_list = []
        if basis_list is None:
            basis_list = []
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
                        os.path.join(self.basis_sets_output_dir, f"{e}_{b}.basis"),
                        "w",
                    ) as fhandle:
                        fhandle.write(bas)
            except KeyError:
                logger.debug(f"element={e}, basis={b} do not exist in the database.")

    def all_to_file(self, sleep_time: float = 1):
        self.to_file(
            element_list=chemical_symbols,
            basis_list=self.list_of_basis_all,
            sleep_time=sleep_time,
        )


class BFD:
    C_URL = "https://github.com/TREX-CoE/BFD-ECP.git"
    list_of_basis_all = [f"v{s}z" for s in "dtq56"]

    def __init__(
        self,
        basis_sets_output_dir: Optional[str] = None,
        pseudo_potential_output_dir: Optional[str] = None,
    ):
        if basis_sets_output_dir is None:
            self.basis_sets_output_dir = None
        else:
            self.basis_sets_output_dir = basis_sets_output_dir
        if pseudo_potential_output_dir is None:
            self.ecp_output_dir = None
        else:
            self.ecp_output_dir = pseudo_potential_output_dir

    def to_file(
        self,
        element_list: Optional[list] = None,
        basis_list: Optional[list] = None,
    ):
        if element_list is None:
            element_list = []
        if basis_list is None:
            basis_list = []
        if self.basis_sets_output_dir is not None:
            os.makedirs(self.basis_sets_output_dir, exist_ok=True)
        if self.ecp_output_dir is not None:
            os.makedirs(self.ecp_output_dir, exist_ok=True)

        # clone the git repository.
        tempdir = pathlib.Path(tempfile.mkdtemp())
        git.Repo.clone_from(self.C_URL, tempdir)

        # basis
        for p in (tempdir / "BASIS" / "GAMESS").glob("*/*"):
            if self.basis_sets_output_dir is None:
                break
            for element, basis in itertools.product(element_list, basis_list):
                if re.match(element, str(p.name)) and re.match(
                    basis, str(p.parent.name)
                ):
                    with open(p, "r") as fhandle:
                        with open(
                            os.path.join(
                                self.basis_sets_output_dir,
                                str(p.name) + f"_{str(p.parent.name)}.basis",
                            ),
                            "w",
                        ) as fhandle_out:
                            fhandle_out.write(fhandle.read())

        # ecp
        for p in (tempdir / "ECP" / "GAMESS").glob("*"):
            if self.ecp_output_dir is None:
                break
            for element in element_list:
                if re.match(element, str(p.name)):
                    with open(p, "r") as fhandle:
                        with open(
                            os.path.join(
                                self.ecp_output_dir,
                                str(p.name) + "_BFD.pseudo",
                            ),
                            "w",
                        ) as fhandle_out:
                            fhandle_out.write(fhandle.read())

    def all_to_file(self, sleep_time: float = 1):
        self.to_file(element_list=chemical_symbols, basis_list=self.list_of_basis_all)


if __name__ == "__main__":
    logger = getLogger("pyturbo")
    logger.setLevel("DEBUG")
    stream_handler = StreamHandler()
    stream_handler.setLevel("DEBUG")
    handler_format = Formatter("%(name)s - %(levelname)s - %(lineno)d - %(message)s")
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)

    downloader_test_dir = os.path.join(pyturbo_root, "tests", "downloader")
    os.chdir(downloader_test_dir)

    # moved to examples
