#!python
# -*- coding: utf-8 -*-

#python modules
import os, sys, shutil

#set logger
from logging import config, getLogger, StreamHandler, Formatter
logger = getLogger('pyturbo')
logger.setLevel("DEBUG")
stream_handler = StreamHandler()
stream_handler.setLevel("DEBUG")
handler_format = Formatter('%(name)s - %(levelname)s - %(lineno)d - %(message)s')
stream_handler.setFormatter(handler_format)
logger.addHandler(stream_handler)

#turbogenius modules
from turbogenius.pyturbo.convertfort10 import Convertfort10

from turbogenius.utils_workflows.env import turbo_genius_root

# convertfort10
prefix="convertfort10"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

logger.info(os.getcwd())
namelist=Convertfort10.read_default_namelist()
convertfort10=Convertfort10(
    in_fort10="fort.10_in",
    out_fort10="fort.10_out",
    namelist=namelist,
)
convertfort10.set_parameter("nx", 20, "&mesh_info")
convertfort10.set_parameter("ny", 20, "&mesh_info")
convertfort10.set_parameter("nz", 20, "&mesh_info")
convertfort10.sanity_check()
convertfort10.generate_input(input_name="convertfort10.input")
convertfort10.run(input_name="convertfort10.input", output_name="out_conv")

flags=convertfort10.check_results(output_names=["out_conv"])
logger.info(f"check_results={flags}")