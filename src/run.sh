#!/bin/bash

> frame0.dump
g++ -O2 -c -o myio_serial.o myio_serial.cpp 
g++ -O2 -c -o md2d_serial.o md2d_serial.cpp
g++ -O2 -o md2d_serial md2d_serial.o myio_serial.o

./md2d_serial 100 0.001 100000 0.5 100 0





# debugging: 
g++ -O2 -g -o md2d.o md2d.cpp
./md2d.o 100 0.001 100000 0.5 100 0






# parallel:
mpicxx -o parallel_md.o parallel_md.cpp && mpirun -np 6 parallel_md.o 100 0.001 10000 0.5 100 0 > thermo