# TD1DIPSolid

A simple code to calculate time propagation of independent particle system, assumed to be Fermion, on a spatially periodic potential.
Most part is written in python. A Fortran library is attached to accelerate the most time-consuming part, time-propagation. The library is not mandatory for run, namely python3-numpy only environment is enough in principle.

# Requirements 
 
* numpy
* matplotlib
* A Fortran compiler, just checked with gfortran (gcc version 9.3.0)
  BLAS, LAPACK, FFTW are required for the Fortran library compile

# Author

* Yasushi SHINOHARA (篠原康)

 
# License
 
"TD1DIPSolid" is under [MIT license](https://en.wikipedia.org/wiki/MIT_License).
