#!/bin/bash

nproc=$1

# manual linking in run/ dir
mpicxx -O3 -c -o myio_mpi.o ../src/myio_mpi.cpp 
mpicxx -O3 -c -o md2d_mpi.o ../src/md2d_mpi.cpp

# manual compilation in run/ dir
mpicxx -O3 -o md2d_mpi md2d_mpi.o myio_mpi.o 