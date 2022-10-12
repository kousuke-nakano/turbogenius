#!/bin/bash
root_dir=`pwd`

cd $root_dir
cd pyturbo_examples; ./all_examples_clean.sh

cd $root_dir
cd turbogenius_examples; ./all_examples_clean.sh