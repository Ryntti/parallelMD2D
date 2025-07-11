Serial program
--------------

execution (Manual):
    cd run/
    ./md2d_serial 100 0.001 10000 0.5 100 0

execution (Thermoplot example):  
    cd src/
    ./run_thermo_serial.sh


MPI program
--------------
execution (CSC Puhti):
    cd run/
    sbatch run_puhti_mpi.sh

execution (Manual):
    cd run/
    mpirun -np 6 ./md2d_mpi 100 0.001 10000 0.5 50 0

execution (Thermoplot example):
    cd run/
    ./run_thermo_mpi.sh 6

execution (gather data on time scaling and plot t(N)):  #This is the main time scaling analysis pipeline!
    cd run/
    ./parallel_timescaling.sh


Hybrid MPI/OpenMP program
--------------
execution (Manual):
    cd run/
    mpirun -np 3 env OMP_NUM_THREADS=2 ./md2d_mpi_omp 100 0.001 10000 0.5 50 0

execution (Thermoplot example):
    cd run/
    ./run_thermo_mpi_omp.sh


Powerlaw analysis
--------------
After running the puhti batch job (run_puhti_mpi.sh) which fills up the file "avg_puhtitimes.dat", powerlaw.py may be run.
execution:
    python3 powerlaw.py

This will create the powerlaw plot, powerlaw.png

