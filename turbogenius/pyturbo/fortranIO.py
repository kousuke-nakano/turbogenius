#!python -u
# -*- coding: utf-8 -*-

"""

pyturbo: fortranIO abstract class

"""

# python modules
from abc import ABC, abstractmethod

# set logger
from logging import getLogger

logger = getLogger("pyturbo").getChild(__name__)


# fortranIO abstract class
class FortranIO(ABC):

    """

    This class is an abstract class for TurboRVB python wrappers.

    Attributes:
    """

    def __init__(
        self,
        # namelist=Namelist()
    ):

        # self.namelist = namelist
        pass

    # commom methods
    def set_parameter(self, parameter, value, namelist=None):
        self.namelist.set_parameter(parameter=parameter, value=value, namelist=namelist)

    def get_parameter(self, parameter, namelist=None):
        return self.namelist.get_parameter(parameter=parameter, namelist=namelist)

    def comment_out(self, parameter):
        self.namelist.comment_out(parameter=parameter)

    def get_parameters(self):
        return self.namelist.get_parameters()

    # abstract methods
    @abstractmethod
    def sanity_check(self):
        pass

    @abstractmethod
    def generate_input(self):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def check_results(self):
        pass

    @abstractmethod
    def read_default_namelist(self):
        pass

    @abstractmethod
    def read_namelist_from_file(file):
        pass

    @classmethod
    @abstractmethod
    def parse_from_default_namelist(self):
        pass

    @classmethod
    @abstractmethod
    def parse_from_file(self):
        pass
