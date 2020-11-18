#gfortran -O0 -shared -fPIC -o Fortlib.so Fortlib.f90
#gfortran -O3 -shared -fPIC -o Fortlib.so Fortlib.f90
#gfortran -mtune=native -march=native -O3 -shared -fPIC -o Fortlib.so Fortlib.f90
gfortran -mtune=native -march=native -fopenmp -O3 -shared -fPIC -o Fortlib.so Fortlib.f90
