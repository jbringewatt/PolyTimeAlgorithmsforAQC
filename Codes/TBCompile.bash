#!/bin/bash
gcc -c ./TBSolver.c -g 
gfortran -c dsygvic.f90 dlaasrt.f90 -llapack -lblas -g  
gfortran -o TBSolver TBSolver.o dsygvic.o dlaasrt.f90 -lm -lblas -llapack -g
rm *.o