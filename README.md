# Polynomial Time Algorithms for Estimating Spectra of Adiabatic Hamiltonians
Codes to calculate eigenvalue gaps for adiabatic quantum computing problems with collections of individually Hamming symmetric wells. See the paper "Polynomial Time Algorithms for Estimating the Spectra of Adiabatic Hamiltonians" for details (https://arxiv.org/abs/1905.07461)  

There are two types of codes - exact solvers for one well, two wells and three wells, and tight binding solvers for 2 or more wells. Many of the codes have a version called codeName_GroverWithPriors.c. For these versions the first well in the input file will be the only one in the final potential of the adiabatic evolution and all other wells will be treated as "priors" which are in the initial Hamiltonian, but not the final Hamiltonian.  All codes take well input files of the same structure and unique parameter input files. See the comments at the top of the individual codes for details.  

Dependencies: Lapack, blas  

Compiler: gcc 

Uses code for implementation of Fix-Heiberger algorithm: https://github.com/cmjiang/xSYGVIC

For questions or comments email: jbringew@terpmail.umd.edu
