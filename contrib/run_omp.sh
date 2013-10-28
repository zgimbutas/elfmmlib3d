#/bin/bash

#
#  Run an OpenMP program 
#
#  Examples:
#
#  run_omp.sh 1 ./int2
#  run_omp.sh 2 ./int2
#
#  Note that the OpenMP stack size is set to 1024MB in this script.
#

PATH=.:$PATH

ulimit -s

# needed for gfortran: libgomp default stack size is too small
export OMP_STACKSIZE=1024M

# needed for gfortran: use version 4.4.2 or above 
# PATH=/opt/gcc-4.4.2/bin:$PATH
# export LD_LIBRARY_PATH=/opt/gcc-4.4.2/lib64
# make -f testfmm3lt.make HOST=linux-gfortran-openmp -j4 


(export OMP_NUM_THREADS=$1; /usr/bin/time $2)

