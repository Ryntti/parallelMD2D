#include <vector>
#include <mpi.h>

struct vec2D {
    double x,y;
    vec2D(double x, double y) : x{x}, y{y} {}
    vec2D() : x{0.0}, y{0.0} {}
};

void read_cmdline_inputs(int argc, char** argv, int *maxt, double *dt, int *nat, int *nuc, double *box, 
                         double *vsc, int *dumprate, int *thermorate, double *k2, double *k3, int *rpp, 
                         int *cpp, int *er, int *ec, int *dims);

void print_energies(const std::vector<double> &local_ep, const std::vector<double> &local_ek, int rank, 
                    double dt, double t, MPI_Comm comm_2d);