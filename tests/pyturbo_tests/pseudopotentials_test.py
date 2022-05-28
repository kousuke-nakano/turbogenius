#!python
# -*- coding: utf-8 -*-
import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_pseudopotentials_examples():
    try:
        import pseudopotentials_examples
        assert True
    except:
        assert False