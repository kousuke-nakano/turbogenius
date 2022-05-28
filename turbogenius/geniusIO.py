#!python -u
# -*- coding: utf-8 -*-

#python modules
import os, sys
from abc import ABC, abstractmethod

#set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from namelist import Namelist

# fortranIO abstract class
class GeniusIO(ABC):

    def __init__(self,
                 ):
        pass

    #abstract methods
    @abstractmethod
    def run_all(self):
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
