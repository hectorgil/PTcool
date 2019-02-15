#!/bin/bash

#compilation with gcc
gcc  main.c functions.c functions_rsd.c nlbias_integrals.c pt_integrals.c pt_kernels.c tns_integrals.c tns_kernels.c cubature.c structures.c  -O3 -lm -fopenmp -o file_RPT_gcc.out

./file_RPT_gcc.out params.ini

#compilation with icc

icc  main.c functions.c functions_rsd.c nlbias_integrals.c pt_integrals.c pt_kernels.c tns_integrals.c tns_kernels.c cubature.c structures.c  -O3 -lm -openmp -o file_RPT_icc.out

./file_RPT_icc.out params.ini

