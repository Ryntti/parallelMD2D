#!/bin/bash

avg_times="avg_wctimes.dat"
wctime_plotscript="wc-n_plots.py"

echo -n "" > $avg_times

for i in {1..6}; do
    total_time=0
    num_runs=10

    for j in {1..10}; do
        # run and capture wc time
        runtime=$(mpirun -np ${i} md2d_mpi 400 0.001 500 0.5 100 0 2>&1 | grep "Elapsed time" | awk '{print $3}')
        total_time=$(echo "$total_time + $runtime" | bc)
    done

    # avg time
    avg_time=$(echo "$total_time / $num_runs" | bc -l)
    echo "Average wall clock time in seconds using $i processes: $avg_time " >> $avg_times
done


cat << EOF > $wctime_plotscript
import matplotlib.pyplot as plt

nprocs, avgs = [], []
with open("$avg_times", "r") as file:
    for i, line in enumerate(file):
        avgt = float(line.split()[-1])
        np = int(line.split()[7])
        avgs.append(avgt)
        nprocs.append(np)

plt.figure()
plt.title('Wall clock time as a function of the number of processes \n using 16000 atoms')
plt.plot(nprocs, avgs, 'r-', label='wc(N)')
plt.xlabel('Number of processes')
plt.ylabel('Wall clock time (s)')
plt.legend()
plt.savefig('wc-n.png', dpi=100)
EOF

python3 ${wctime_plotscript}