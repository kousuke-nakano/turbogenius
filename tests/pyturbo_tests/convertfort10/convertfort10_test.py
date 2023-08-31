#!python
# -*- coding: utf-8 -*-
import os, sys

# turbogenius modules
from turbogenius.pyturbo.convertfort10 import Convertfort10

root_dir = os.path.dirname(__file__)
os.chdir(root_dir)

namelist = Convertfort10.read_default_namelist()
convertfort10 = Convertfort10(
    in_fort10="fort.10_in",
    out_fort10="fort.10_out",
    namelist=namelist,
)
convertfort10.set_parameter("nx", 20, "&mesh_info")
convertfort10.set_parameter("ny", 25, "&mesh_info")
convertfort10.set_parameter("nz", 30, "&mesh_info")
assert convertfort10.get_parameter("nx") == 20
assert convertfort10.get_parameter("ny") == 25
assert convertfort10.get_parameter("nz") == 30
convertfort10.generate_input(input_name="convertfort10.input")
convertfort10.run(input_name="convertfort10.input", output_name="out_conv")

flags = convertfort10.check_results(output_names=["out_conv"])
assert flags

os.chdir(root_dir)
