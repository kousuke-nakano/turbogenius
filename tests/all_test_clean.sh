#!/bin/bash
root_dir=`pwd`

cd $root_dir
cd pyturbo_tests; ./all_test_clean.sh

cd $root_dir
cd turbogenius_test; ./all_test_clean.sh