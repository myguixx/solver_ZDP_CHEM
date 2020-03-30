#! /bin/bash
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
gfortran dvode_f90_m.F90 zdplaskin_m.F90 test_2reac.F90 bolsig_x86_64.so -o zdplaskine