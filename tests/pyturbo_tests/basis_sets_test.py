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