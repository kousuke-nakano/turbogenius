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

#pyturbo modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from turbogenius.pyturbo.prep import Prep
from turbogenius.utils_workflows.env import turbo_genius_root

#prep
prefix="prep"
example_root_dir=os.path.join(turbo_genius_root, "examples", "pyturbo_examples")
if os.path.isdir(os.path.join(example_root_dir, prefix)):
    shutil.rmtree(os.path.join(example_root_dir, prefix))
shutil.copytree(os.path.join(example_root_dir, "all_input_files", prefix), os.path.join(example_root_dir, prefix))
os.chdir(os.path.join(example_root_dir, prefix))

namelist = Prep.read_default_namelist()
kpoints=[[
            [0.0, 0.0, 0.0, 1.0],
            [0.2, 0.2, 0.2, 0.5],
            [0.3, 0.3, 0.3, 0.4],
         ],
         [
            [-0.0, 0.0, 0.0, 1.0],
            [-0.2, 0.2, 0.2, 0.5],
            [-0.3, 0.3, 0.3, 0.4],
         ]
         ]

prep = Prep(
    in_fort10="fort.10",
    namelist=namelist,
    twist_average=2,
)
prep.sanity_check()
prep.manual_kpoints = kpoints
prep.generate_input(input_name="prep.input")
#prep.run(input_name="prep.input", output_name="out_prep")

#flags = prep.check_results(output_names=["out_prep"])
#logger.info(f"check_results={flags}")