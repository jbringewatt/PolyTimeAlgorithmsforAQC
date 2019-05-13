#!/bin/bash
gcc -c ./TBSolver_GroverWithPriorsFinal.c
gfortran -c dsygvic.f90 dlaasrt.f90 -llapack -lblas
gfortran -o TBSolver_GroverWithPriorsFinal TBSolver_GroverWithPriorsFinal.o dsygvic.o dlaasrt.f90 -lm -lblas -llapack
rm *.o