#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys

#set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo')
logger.setLevel("DEBUG")
stream_handler = StreamHandler()
stream_handler.setLevel("DEBUG")
handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from turbogenius.pyturbo.convertfort10mol import Convertfort10mol

# convertfort10mol
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "convertfort10mol"))
convertfort10mol = Convertfort10mol.parse_from_default_namelist(in_fort10="fort.10_in")
convertfort10mol.sanity_check()
convertfort10mol.generate_input(input_name="convertfort10mol.input")
convertfort10mol.run(input_name="convertfort10mol.input", output_name="out_mol")
flags = convertfort10mol.check_results(output_names=["out_mol"])
logger.info(f"check_results={flags}")