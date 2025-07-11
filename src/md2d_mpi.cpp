#include <vector>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <mpi.h>
// #include <omp.h>
#include <numeric>
#include <random>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>

#include "myio_mpi.h"

const int n_neigh = 4; // number of neighbors always 4 due to 2d and pbcs
const double d = 1.0;

void bcast_inputdata(int master, int *maxt, double *dt, int *nat, int *nuc, double *box, double *vsc, int *dumprate, int *thermorate, 
                     double *k2, double *k3, int *rpp, int *cpp, int *er, int *ec, MPI_Comm comm_2d) {
    // Broadcast the read input values to other processes from the master process where they've been read in from the command line
    MPI_Bcast(maxt, 1, MPI_INT, master, comm_2d);
    MPI_Bcast(dt, 1, MPI_DOUBLE, master, comm_2d);
    MPI_Bcast(nat, 1, MPI_INT, master, comm_2d);
    MPI_Bcast(nuc, 1, MPI_INT, master, comm_2d);

    MPI_Bcast(box, 1, MPI_DOUBLE, master, comm_2d);
    MPI_Bcast(vsc, 1, MPI_DOUBLE, master, comm_2d);
    MPI_Bcast(dumprate, 1, MPI_INT, master, comm_2d);
    MPI_Bcast(thermorate, 1, MPI_INT, master, comm_2d);

    MPI_Bcast(k2, 1, MPI_DOUBLE, master, comm_2d);
    MPI_Bcast(k3, 1, MPI_DOUBLE, master, comm_2d);

    MPI_Bcast(rpp, 1, MPI_INT, master, comm_2d);
    MPI_Bcast(cpp, 1, MPI_INT, master, comm_2d);
    MPI_Bcast(er, 1, MPI_INT, master, comm_2d);
    MPI_Bcast(ec, 1, MPI_INT, master, comm_2d);
}

void create_system(std::vector<vec2D> &local_pos, std::vector<vec2D> &local_vel, double vsc, 
                   int seed, int num_cols, int num_rows, int start_col, int start_row) {
    /* nuc is the number of atoms in one dimension of the whole square lattice.
       pos and vel have been allocated already to include the atoms in the subdomain,
       as well as the "ghost zones", i.e. the edges from the neighbor processes. The ghost zones
       essentially function just as the receive buffers in this simulation */
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> fdistrib(-0.5, 0.5);

    for (int row = 0; row < num_rows; ++row) {
        for (int col = 0; col < num_cols; ++col) {

            int i_extended = (row + 1) * (num_cols + 2) + col + 1;  // linear index for the 1D vector

            local_pos[i_extended].x = start_col + col;     // Global x position
            local_pos[i_extended].y = start_row + row;     // Global y position

            int i = col + row * num_cols;

            local_vel[i].x = vsc * fdistrib(gen);   // assign random velocitites to all atoms
            local_vel[i].y = vsc * fdistrib(gen);
        }
    }
}


void get_local_comv(std::vector<vec2D> &local_vel, vec2D &local_vcom) {
    /* Obtain the total velocity sum of the subdomain */
    double local_velxsum = std::accumulate(local_vel.begin(), local_vel.end(), 0.0, [](double acc, const vec2D &v){return acc + v.x;});
    double local_velysum = std::accumulate(local_vel.begin(), local_vel.end(), 0.0, [](double acc, const vec2D &v){return acc + v.y;});
    local_vcom.x = local_velxsum;
    local_vcom.y = local_velysum;
}

void subtract_com_velocity(std::vector<vec2D>& local_vel, int nat, MPI_Comm comm_2d) {
    /* Calculate and subtract the c.o.m velocity of the system from the velocities of each atom */
    vec2D local_vcom, global_vcom;
    get_local_comv(local_vel, local_vcom);

    // compute global c.o.m. velocity by summing all local vcoms
    MPI_Reduce(&local_vcom.x, &global_vcom.x, 1, MPI_DOUBLE, MPI_SUM, 0, comm_2d);
    MPI_Reduce(&local_vcom.y, &global_vcom.y, 1, MPI_DOUBLE, MPI_SUM, 0, comm_2d);

    global_vcom.x /= static_cast<double>(nat);
    global_vcom.y /= static_cast<double>(nat);

    // Broadcast global c.o.m. velocity to all processes
    MPI_Bcast(&global_vcom.x, 1, MPI_DOUBLE, 0, comm_2d);
    MPI_Bcast(&global_vcom.y, 1, MPI_DOUBLE, 0, comm_2d);

    // Subtract the global center of mass velocity from each local velocity
    for (auto& v : local_vel) {
        v.x -= global_vcom.x;
        v.y -= global_vcom.y;
    }
}

void get_neigh_list(std::vector<std::array<int, n_neigh>>& local_nlist, int nuc, int num_cols, int num_rows, int start_col, int start_row) {
    // Assuming local nlist is already allocated appropriately
    /* Again, subdomain lattice ghost zones are included only to make indexing correspond to local_pos. In fact, the ghost zone indices
       of the nlist are never accessed nor assigned any values. */

    for (int row = 0; row < num_rows; ++row) {
        for (int col = 0; col < num_cols; ++col) {
            
            int i = (row + 1) * (num_cols + 2) + col + 1;   // only assign values to the "inner submatrix", which excludes the ghost zones

            local_nlist[i][0] = i - 1;  // Left neighbor
            local_nlist[i][1] = i + 1; // Right neighbor
            local_nlist[i][2] = i - num_cols - 2; // Lower neighbor
            local_nlist[i][3] = i + num_cols + 2; // Upper neighbor
        }
    }
}


void accel(std::vector<vec2D> &local_pos, int i_extended, double *u, vec2D &a, std::vector<std::array<int, n_neigh>> &local_nlist, double box, double k2, double k3, double d) {
    /* Calculate the potential energy u and acceleration a of atom whose index is i_extended in local_pos and nlist */
    a.x = 0.0; a.y = 0.0;
    *u = 0.0;   // reset the potential energy for a new calculation

    double u2 = 0.0; double u3 = 0.0;
    for (int j = 0; j < n_neigh; j++) {
        int in = local_nlist[i_extended][j];

        double dx = local_pos[in].x - local_pos[i_extended].x;
        double dy = local_pos[in].y - local_pos[i_extended].y;

        // apply periodic boundary conditions
        if (dx < -0.5 * box) dx += box; 
        if (dx >= 0.5 * box) dx -= box;
        if (dy < -0.5 * box) dy += box;
        if (dy >= 0.5 * box) dy -= box;

        double r = std::sqrt(dx*dx + dy*dy);
        // Calculate the potential energy. Both harmonic and anharmonic terms are used
        double c = r-d;
        u2 += 0.5 * k2 * c*c;
        u3 += 1.0 / 3.0 * k3 * c*c*c;  // 1.0 / 3.0

        double dxr_term = dx/r;    

        double fx2 = k2 * c * dxr_term;
        double fx3 = k3 * c*c * dxr_term;

        // update accelerations
        a.x += fx2 + fx3;

        double dyr_term = dy/r;

        double fy2 = k2 * c * dyr_term;
        double fy3 = k3 * c*c * dyr_term;

        a.y += fy2 + fy3;
    }
    *u = 0.5 * (u2 + u3); // eliminate double counting by dividing by two
}

void update_all_accelerations(std::vector<vec2D> &local_pos, std::vector<double> &local_ep, std::vector<vec2D> &local_a, 
                                std::vector<std::array<int, n_neigh>> &local_nlist, double box, double k2, double k3, double d, 
                                int num_cols, int num_rows) {
    /* Update all accelerations of a local subdomain */
    for (int row = 0; row < num_rows; ++row) {
        for (int col = 0; col < num_cols; ++col) {

            int i_extended = (row + 1) * (num_cols + 2) + col + 1;  // index of the atom for the arrays with extra padding 
            int i = col + row*num_cols;                             // index of the atom for normal arrays

            accel(local_pos, i_extended, &local_ep[i], local_a[i], local_nlist, box, k2, k3, d);
        }
    }
}

void handle_pbcs(std::vector<vec2D> &local_pos, int i, double box) {
    /* Apply periodic boundaries */
    if (local_pos[i].x < 0.0) local_pos[i].x += box;
    if (local_pos[i].y < 0.0) local_pos[i].y += box;
    if (local_pos[i].x >= box) local_pos[i].x -= box;
    if (local_pos[i].y >= box) local_pos[i].y -= box;
}

void leapfrog(std::vector<vec2D> &local_pos, std::vector<vec2D> &vel, std::vector<vec2D> &vel0, 
              std::vector<vec2D> &a, std::vector<double> &ep, std::vector<double> &ek, int nat, 
              double dt, double box, int num_cols, int num_rows) {
    /* Perform Leapfrog time itegration by one timestep. Also handle periodic boundaries 
       and update kinetic and potential energies. This is performed locally in each subprocess */
    for (int row = 0; row < num_rows; ++row) {
        for (int col = 0; col < num_cols; ++col) {

            int i = col + row*num_cols;

            vel[i].x += dt * a[i].x;
            vel[i].y += dt * a[i].y;

            int i_extended = (row + 1) * (num_cols + 2) + col + 1;

            local_pos[i_extended].x += dt * vel[i].x;
            local_pos[i_extended].y += dt * vel[i].y;

            handle_pbcs(local_pos, i_extended, box);

            double vx = 0.5 * (vel0[i].x + vel[i].x);
            double vy = 0.5 * (vel0[i].y + vel[i].y);

            ek[i] = 0.5 * (vx*vx + vy*vy);
        }
    }
}


void exchange_boundary_data(std::vector<vec2D> &local_pos, int num_cols, MPI_Datatype Vec2D_Type, MPI_Datatype Col_Type, MPI_Comm comm_2d) {
    /*  */
    int axis = 0;   // vertical axis message passing first
    // Get neighboring ranks
    int north, south; // the IDs of neighboring processes (in the 2d cart topology) are saved for message passing
    
    int stride = 1;
    MPI_Cart_shift(comm_2d, axis, stride, &south, &north);    // Neighbors in the vertical dimension

    MPI_Sendrecv(&local_pos[local_pos.size() - 1 - 2*num_cols - 2], num_cols, Vec2D_Type, north, 10,
                    &local_pos[1], num_cols, Vec2D_Type, south, 10,
                    comm_2d, MPI_STATUS_IGNORE);


    MPI_Sendrecv(&local_pos[1 + num_cols + 2], num_cols, Vec2D_Type, south, 100,
                    &local_pos[local_pos.size() - 1 - num_cols], num_cols, Vec2D_Type, north, 100,
                    comm_2d, MPI_STATUS_IGNORE);


    axis = 1;   // horizontal message passing

    int west, east;
    MPI_Cart_shift(comm_2d, axis, stride, &west, &east);  // Neighbors in the first dimension

    MPI_Sendrecv(&local_pos[1 + num_cols + 2], 1, Col_Type, west, 30,
                    &local_pos[2*num_cols + 3], 1, Col_Type, east, 30,
                    comm_2d, MPI_STATUS_IGNORE);


    MPI_Sendrecv(&local_pos[2*num_cols + 2], 1, Col_Type, east, 20,
                    &local_pos[num_cols + 2], 1, Col_Type, west, 20,
                    comm_2d, MPI_STATUS_IGNORE);
}



int main(int argc, char** argv) {

    int ntasks, master, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    master = 0;

    int dims[2] = {0, 0}; int periods[2] = {1, 1};    // apply periodic boundary conditions      
    int coords[2];      
    int reorder = 1;  

    // create the topology
    MPI_Dims_create(ntasks, 2, dims);
    MPI_Comm comm_2d;       // declare the communicator for 2d cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm_2d);

    MPI_Comm_size(comm_2d, &ntasks);
    MPI_Comm_rank(comm_2d, &rank);
    MPI_Cart_coords(comm_2d, rank, 2, coords);   // get the coordinates of the process
    
    //-------------------------------------------------------
    double vsum, vave, rn;
    vec2D local_vcom;
    double k2, k3;
    double vsc;          // mean initial velocity
    double box;          // system size (length)
    // double *box_cellvectors;
    double dt;           // time step
    int nat, nuc;        // nat = number of atoms, nuc the amount of unit cells
    int maxt;            // number of time steps simulated
    int seed = 1234 * 11 * rank;
    int thermorate, dumprate;

    int rows_per_proc, cols_per_proc, extra_rows, extra_cols;
    //-------------------------------------------------------


    if (rank == master) {
        read_cmdline_inputs(argc, argv, &maxt, &dt, &nat, &nuc, &box, &vsc, &dumprate, &thermorate, &k2, &k3, 
                    &rows_per_proc, &cols_per_proc, &extra_rows, &extra_cols, dims);
    }

    bcast_inputdata(master, &maxt, &dt, &nat, &nuc, &box, &vsc, &dumprate, &thermorate, &k2, &k3, 
                    &rows_per_proc, &cols_per_proc, &extra_rows, &extra_cols, comm_2d);

    
    // next, starting and ending rows and columns of the atoms of the current subdomain
    int start_row = coords[0] * rows_per_proc + (coords[0] < extra_rows ? coords[0] : extra_rows);
    int start_col = coords[1] * cols_per_proc + (coords[1] < extra_cols ? coords[1] : extra_cols);

    int end_row = start_row + rows_per_proc + (coords[0] < extra_rows ? 1 : 0) - 1;
    int end_col = start_col + cols_per_proc + (coords[1] < extra_cols ? 1 : 0) - 1;

    // number of rows and columns in the current subdomain
    int num_cols = end_col - start_col + 1;
    int num_rows = end_row - start_row + 1;
    int num_atoms = num_rows * num_cols;

    // int *atom_counts;
    // if (rank == master) {
    //     atom_counts = new int[ntasks];
    // }

    //MPI_Gather(&num_atoms, 1, MPI_INT, atom_counts, 1, MPI_INT, master, comm_2d);   // gather all subdomain atom counts to an array
    // MPI_Offset init_dump_offset = calculate_offset(rank);

    int extra_size = (num_cols + 2) * (num_rows + 2);
 
    std::vector<vec2D> local_pos(extra_size, vec2D{}), local_vel(num_atoms, vec2D{}), local_vel0(num_atoms, vec2D{}), local_a(num_atoms, vec2D{});
    std::vector<double> local_ep(num_atoms), local_ek(num_atoms);
    std::vector<std::array<int, n_neigh>> local_nlist(extra_size);

    create_system(local_pos, local_vel, vsc, seed, num_cols, num_rows, start_col, start_row);    // initializes positions and velocities of the atoms
    get_neigh_list(local_nlist, nuc, num_cols, num_rows, start_col, start_row);   // set up the neighbor lists

    subtract_com_velocity(local_vel, nat, comm_2d);

    // MPI_File dumpfile;
    // MPI_File_open(comm_2d, "dump.xyz", MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND, MPI_INFO_NULL, &dumpfile);


    // define custom MPI types for vec2Ds and columns of a 2D array stored as a 1D array of vec2D structs. The latter is for receive buffer only
    MPI_Datatype Vec2D_Type;
    MPI_Type_contiguous(2, MPI_DOUBLE, &Vec2D_Type);
    MPI_Type_commit(&Vec2D_Type);

    MPI_Datatype Col_Type;
    MPI_Type_vector(num_rows, 1, num_cols + 2, Vec2D_Type, &Col_Type);
    MPI_Type_commit(&Col_Type);

    auto t0 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < maxt; t++) {    // time integration
        local_vel0 = local_vel;

        exchange_boundary_data(local_pos, num_cols, Vec2D_Type, Col_Type, comm_2d);

        update_all_accelerations(local_pos, local_ep, local_a, local_nlist, box, k2, k3, d, num_cols, num_rows);

        leapfrog(local_pos, local_vel, local_vel0, local_a, local_ep, local_ek, nat, dt, box, num_cols, num_rows);


        if (t%thermorate == 0 && thermorate > 0) print_energies(local_ep, local_ek, rank, dt, t, comm_2d);

        // if (t%dumprate == 0 && dumprate > 0) write_dump(id, t, num_atoms, nat, init_dump_offset, Vec2D_Type, dumpfile, comm_2d);

        MPI_Barrier(comm_2d);
    }


    // delete heap-allocated arrays
    // delete[] atom_counts;
    // if (id == master) delete[] box_cellvectors;

    if (rank == master) {
        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> runtime = t1 - t0;
        std::cout << "Elapsed time: " << runtime.count() << " s\n";
    }

    // close dump file
    // MPI_File_close(&dumpfile);

    MPI_Type_free(&Vec2D_Type);
    MPI_Type_free(&Col_Type);
    MPI_Finalize();
    return 0;
}