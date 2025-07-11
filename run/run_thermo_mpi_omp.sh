#!/bin/bash

#nprocs=$1
> thermo.dat

mpirun -np 3 env OMP_NUM_THREADS=2 ./md2d_mpi_omp 100 0.001 10000 0.5 50 0 > thermo.dat

python3 thermoplot.py