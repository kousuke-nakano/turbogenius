#!python
# -*- coding: utf-8 -*-

import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_convertfort10mol_examples():
    try:
        import convertfort10mol_examples
        assert True
    except:
        assert False

def test_lrdmc_examples():
    try:
        import lrdmc_examples
        assert True
    except:
        assert False

def test_makefort10_examples():
    try:
        import makefort10_examples
        assert True
    except:
        assert False

def test_prep_examples():
    try:
        import prep_examples
        assert True
    except:
        assert False

def test_sequential_workflow_examples():
    try:
        import sequential_workflow_examples
        assert True
    except:
        assert False

"""
def test_trexio_to_turborvb_examples():
    try:
        import trexio_to_turborvb_examples
        assert True
    except:
        assert False
"""

def test_vmc_examples():
    try:
        import vmc_examples
        assert True
    except:
        assert False

def test_vmcopt_examples():
    try:
        import vmcopt_examples
        assert True
    except:
        assert False