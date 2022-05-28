#!/bin/bash

# files
find ./ -name "fort.10_ave*" | xargs rm
find ./ -name "fort.10_new*" | xargs rm
find ./ -name "fort.11" | xargs rm
find ./ -name "fort.12*" | xargs rm
find ./ -name "fort.12*" | xargs rm
find ./ -name "fort.2*" | xargs rm
find ./ -name "*.pseudo" | xargs rm
find ./ -name "*.basis" | xargs rm
find ./ -name "out_*" | xargs rm
find ./ -name "*.out" | xargs rm
find ./ -name "*.dat" | xargs rm
find ./ -name "parminimized.d" | xargs rm
find ./ -name "rand*" | xargs rm
find ./ -name "*.d" | xargs rm
find ./ -name "out_*" | xargs rm
find ./ -name "*ccECP*" | xargs rm
find ./ -name "run*.sh" | xargs rm
find ./ -name "*.png" | xargs rm
find ./ -name "structure.xsf" | xargs rm

# dirs
find ./ -name "evsa*" | xargs rm -r
find ./ -name "alat*" | xargs rm -r
find ./ -name "ave_temp" | xargs rm -r
find ./ -name "parameters_graphs" | xargs rm -r
