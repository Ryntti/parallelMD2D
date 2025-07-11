TARGET_mpi = run/md2d_mpi
TARGET_serial = run/md2d_serial
TARGET_mpi_omp = run/md2d_mpi_omp

OBJS_mpi = run/md2d_mpi.o run/myio_mpi.o	# object files (mpi)
OBJS_serial = run/md2d_serial.o run/myio_serial.o	# object files (serial)
OBJS_mpi_omp = run/md2d_mpi_omp.o run/myio_mpi_omp.o  # object files (mpi & OMP)

OBJS_ALL = $(OBJS_mpi) $(OBJS_serial) $(OBJS_mpi_omp)

CXX = g++
MPI_CXX = mpicxx
CXXFLAGS = -O3 
OMPFLAG = -fopenmp

OUT = thermo.dat dump.xyz frame0.dump wc-n.png wc-n_puhti.png avg_wctimes.dat

# Serial target build rule ----------------------------------
serial:	$(TARGET_serial)

$(TARGET_serial): $(OBJS_serial)
	$(CXX) $(CXXFLAGS) -o $(TARGET_serial) $(OBJS_serial)

# rules for creating the individual serial object files from the source files
run/md2d_serial.o: src/md2d_serial.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run/myio_serial.o: src/myio_serial.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAG) -c -o $@ $<
#--------------------------------------------------------------



# MPI with OpenMP ----------------------------------------------

mpi_omp: $(TARGET_mpi_omp)
	
$(TARGET_mpi_omp): $(OBJS_mpi_omp)
	$(MPI_CXX) $(CXXFLAGS) $(OMPFLAG) -o $(TARGET_mpi_omp) $(OBJS_mpi_omp)


# rules for creating the individual mpi object files from the source files
run/md2d_mpi_omp.o: src/md2d_mpi_omp.cpp
	$(MPI_CXX) $(CXXFLAGS) -c -o $@ $<

run/myio_mpi_omp.o:	src/myio_mpi.cpp
	$(MPI_CXX) $(CXXFLAGS) $(OMPFLAG) -c -o $@ $<
# -------------------------------------------------------------





# MPI target build rule ---------------------------------------
mpi: $(TARGET_mpi)

$(TARGET_mpi): $(OBJS_mpi)
	$(MPI_CXX) $(CXXFLAGS) -o $(TARGET_mpi) $(OBJS_mpi)

# rules for creating the individual mpi object files from the source files
run/md2d_mpi.o: src/md2d_mpi.cpp
	$(MPI_CXX) $(CXXFLAGS) -c -o $@ $<

run/myio_mpi.o: src/myio_mpi.cpp
	$(MPI_CXX) $(CXXFLAGS) -c -o $@ $<
# --------------------------------------------------------------



clean:
	rm -f $(OBJS_ALL) $(TARGET_mpi) $(TARGET_serial) $(TARGET_mpi_omp)