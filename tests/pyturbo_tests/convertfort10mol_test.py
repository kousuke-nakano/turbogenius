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