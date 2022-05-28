#!python
# -*- coding: utf-8 -*-

import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_trexio_to_turborvb_examples():
    try:
        import trexio_to_turborvb_examples
        assert True
    except:
        assert False