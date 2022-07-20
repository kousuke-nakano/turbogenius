#!python
# -*- coding: utf-8 -*-

import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def test_sequential_workflow_examples():
    try:
        import sequential_workflow_examples
        assert True
    except:
        assert False