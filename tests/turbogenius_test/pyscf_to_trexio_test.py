#!python
# -*- coding: utf-8 -*-

import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_pyscf_to_trexio_examples():
    try:
        import pyscf_to_trexio_examples
        assert True
    except:
        assert False