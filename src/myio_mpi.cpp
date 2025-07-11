#include <vector>
#include <numeric>
#include <iostream>
#include "myio_mpi.h"

void read_cmdline_inputs(int argc, char** argv, int *maxt, double *dt, int *nat, int *nuc, double *box, 
                         double *vsc, int *dumprate, int *thermorate, double *k2, double *k3, int *rpp, 
                         int *cpp, int *er, int *ec, int *dims) {
    /* Only executed in the master process. Read input params from the command line. */
    int arg_ind = 1;

    if (argc > arg_ind) *nuc = atoi(argv[arg_ind++]);
    if (argc > arg_ind) *dt = atof(argv[arg_ind++]);
    if (argc > arg_ind) *maxt = atoi(argv[arg_ind++]);
    if (argc > arg_ind) *vsc = atof(argv[arg_ind++]);

    // if not specified, use default values (thermo 100, dump 0)
    *thermorate = (argc > arg_ind) ? atoi(argv[arg_ind++]) : 100;
    *dumprate = (argc > arg_ind) ? atoi(argv[arg_ind++]) : 0;

    // if not specified, use default values
    *k2 = (argc > arg_ind) ? atoi(argv[arg_ind++]) : 1.0;
    *k3 = (argc > arg_ind) ? atoi(argv[arg_ind++]) : 0.1;

    int intnuc = *nuc;

    *nat = intnuc*intnuc;
    *box = static_cast<double>(intnuc);

    *rpp = intnuc / dims[0];
    *cpp = intnuc / dims[1];

    *er = intnuc % dims[0];
    *ec = intnuc % dims[1];
}


void print_energies(const std::vector<double> &local_ep, const std::vector<double> &local_ek, int rank, 
                    double dt, double t, MPI_Comm comm_2d) {
    /* Calculate the total kinetic and total potential energies by summing the local individual potential and kinetic energies of 
    the atoms together into local energy sums, and then sum all of these local sums into the global sums and print them */
    
    double local_epsum = std::accumulate(local_ep.begin(), local_ep.end(), 0.0);
    double local_eksum = std::accumulate(local_ek.begin(), local_ek.end(), 0.0);

    double global_epsum = 0.0; double global_eksum = 0.0;

    MPI_Reduce(&local_epsum, &global_epsum, 1, MPI_DOUBLE, MPI_SUM, 0, comm_2d);
    MPI_Reduce(&local_eksum, &global_eksum, 1, MPI_DOUBLE, MPI_SUM, 0, comm_2d);

    if (rank == 0) {
        std::cout << dt * t << " " << global_epsum + global_eksum << " " << global_epsum
                  << " " << global_eksum << std::endl;
    }
}   