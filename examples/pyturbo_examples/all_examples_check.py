#!python
# -*- coding: utf-8 -*-
import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_basis_sets_examples():
    try:
        import basis_sets_examples
        assert True
    except:
        assert False


def test_convertfort10_examples():
    try:
        import convertfort10_examples
        assert True
    except:
        assert False


def test_convertfort10mol_examples():
    try:
        import convertfort10mol_examples
        assert True
    except:
        assert False


def test_fort10_examples():
    try:
        import fort10_examples
        assert True
    except:
        assert False


def test_lrdmc_examples():
    try:
        import lrdmc_examples
        assert True
    except:
        assert False


def test_lrdmcopt_examples():
    try:
        import lrdmcopt_examples
        assert True
    except:
        assert False


def test_makefort10_examples():
    try:
        import makefort10_examples
        assert True
    except:
        assert False

def test_namelist_examples():
    try:
        import namelist_examples
        assert True
    except:
        assert False


def test_prep_examples():
    try:
        import prep_examples
        assert True
    except:
        assert False


def test_pseudopotentials_examples():
    try:
        import pseudopotentials_examples
        assert True
    except:
        assert False


def test_structures_examples():
    try:
        import structures_examples
        assert True
    except:
        assert False


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