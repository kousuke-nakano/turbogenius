#!python -u
# -*- coding: utf-8 -*-

"""

GeniusIO abstract class

ToDo:
    To mix manual (by hand) and automatic (turbogenius) calculations,
    a "parse" function should be implemented in the geniusIO (and its children).
    Indeed, a function to purse an input file to generate the corresponding genius
    instance and the pickled file is needed to edit (by turbogenius) input and output files
    generated manually. This is also important for turboworkflow pakcage.

"""

#python modules
import os, sys
from abc import ABC, abstractmethod

#set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('Turbo-Genius').getChild(__name__)

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# GeniusIO abstract class
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
