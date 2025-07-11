#!/bin/bash

# manual linking in run/ dir
g++ -O3 -c -o myio_serial.o ../src/myio_serial.cpp 
g++ -O3 -c -o md2d_serial.o ../src/md2d_serial.cpp

# manual compilation in run/ dir
g++ -O3 -o md2d_serial md2d_serial.o myio_serial.o