#!/bin/bash

# execution in run/ dir
./md2d_serial 100 0.001 10000 0.5 100 0 > thermo.dat

python3 thermoplot.py