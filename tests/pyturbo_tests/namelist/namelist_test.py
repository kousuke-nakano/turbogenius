#!python
# -*- coding: utf-8 -*-
import os
from turbogenius.pyturbo.utils.env import pyturbo_data_dir
from turbogenius.pyturbo.namelist import Namelist


def test_namelist():
    # parse default value of prep
    prep_default_file = os.path.join(pyturbo_data_dir, "prep", "prep.input")
    namelist = Namelist.parse_namelist_from_file(prep_default_file)

    # check set and get parameters
    namelist.set_parameter(parameter="iopt", value=2)
    assert namelist.get_parameter(parameter="iopt") == 2

    # check set and get parameters
    namelist.comment_out("iopt")
    assert namelist.get_parameter(parameter="!iopt") == 2
