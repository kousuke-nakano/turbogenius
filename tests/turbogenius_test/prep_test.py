#!python
# -*- coding: utf-8 -*-

import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_prep_examples():
    try:
        import prep_examples
        assert True
    except:
        assert False