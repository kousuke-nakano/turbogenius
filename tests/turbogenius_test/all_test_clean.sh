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
find ./ -name "*.log" | xargs rm
find ./ -name "*.pkl" | xargs rm
find ./ -name "*.bin" | xargs rm
find ./ -name "submit.sh" | xargs rm
find ./ -name "turbogenius.o*" | xargs rm
find ./ -name "*.hdf5" | xargs rm
find ./ -name "*.output" | xargs rm
find ./ -name "*.chk" | xargs rm

# dirs
find ./ -name "evsa*" | xargs rm -r
find ./ -name "alat*" | xargs rm -r
find ./ -name "ave_temp" | xargs rm -r
find ./ -name "pkl" | xargs rm -r
find ./ -name "parameters_graphs" | xargs rm -r
find ./ -name "turborvb.scratch" | xargs rm -r