#!/bin/bash

nprocs=$1
> thermo.dat

mpirun -np ${nprocs} ./md2d_mpi 100 0.001 10000 0.5 50 0 > thermo.dat

python3 thermoplot.py