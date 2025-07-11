Serial program
--------------
Compilation option 1 (Makefile):
    make serial     # while in root dir

Compilation option 2 (Manual):    
    cd src/
    g++ -O3 -c -o ../run/myio_serial.o myio_serial.cpp         # compile serial io function module and save executable to run/
    g++ -O3 -c -o ../run/md2d_serial.o md2d_serial.cpp         # compile serial main program and save executable to run/
    g++ -O3 -o ../run/md2d_serial ../run/md2d_serial.o ../run/myio_serial.o     # link the above objects to create the executable



MPI program
--------------
Compilation option 1 (Makefile):
    make mpi     # while in root dir


Compilation option 2 (Manual):    
    cd src/
    mpicxx -O3 -c -o ../run/myio_mpi.o myio_mpi.cpp             # compile serial io function module and save executable to run/
    mpicxx -O3 -c -o ../run/md2d_mpi.o md2d_mpi.cpp             # compile serial main program and save executable to run/
    mpicxx -O3 -o ../run/md2d_mpi ../run/md2d_mpi.o ../run/myio_mpi.o   # link the above objects to create the executable



Hybrid MPI/OpenMP program
--------------
Compilation option 1 (Makefile):
    make mpi_omp     # while in root dir


Compilation option 2 (Manual):    
    cd src/
    mpicxx -O3 -c -o ../run/myio_mpi myio_mpi.cpp                # compile serial io function module and save executable to run/
    mpicxx -O3 -c -o ../run/md2d_mpi_omp.o md2d_mpi_omp.cpp      # compile serial main program and save executable to run/
    mpicxx -O3 -o ../run/md2d_mpi_omp ../run/md2d_mpi_omp.o ../run/myio_mpi_omp.o   # link the above objects to create the executable
