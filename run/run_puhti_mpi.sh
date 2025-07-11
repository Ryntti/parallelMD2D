#!/bin/bash
#SBATCH --account=akuronen
#SBATCH -J run_thpc_project
#SBATCH -o out_md2d
#SBATCH -e error_md2d
#SBATCH -n 20
#SBATCH -t 0-00:15:00
#SBATCH --mem-per-cpu=100
#SBATCH -p test

module load gcc/11.3.0 openmpi/4.1.4 fftw/3.3.10-mpi-omp

avg_times="avg_puhtitimes.dat"
wctime_plotscript="wc-n_plots.py"

echo -n "" > $avg_times

for np in 2 4 6 8 10 12 14 16 18 20; do
    total_time=0
    num_runs=10

    for j in {1..10}; do
        runtime=$(srun -n ${np} md2d_mpi 100 0.001 10000 0.5 100 0 2>&1 | grep "Elapsed time" | awk '{print $3}')
        total_time=$(echo "$total_time + $runtime" | bc)
    done
    avg_time=$(echo "$total_time / $num_runs" | bc -l)
    echo "Average wall clock time in seconds using $np processes: $avg_time " >> $avg_times
done


cat << EOF > $wctime_plotscript
import matplotlib.pyplot as plt

nprocs, avgs = [], []
with open("$avg_times", "r") as file:
    for i, line in enumerate(file):
        avgt = float(line.split()[-1])
        avgs.append(avgt)
        np = int(line.split()[7])
        nprocs.append(np)

plt.figure()
plt.title('Wall clock time as a function of the number of processes \n using 10000 atoms')
plt.plot(nprocs, avgs, 'r-', label='wc(N)')
plt.xlabel('Number of processes')
plt.ylabel('Wall clock time (s)')
plt.legend()
plt.savefig('wc-n_puhti.png', dpi=100)
EOF

python3 ${wctime_plotscript}
