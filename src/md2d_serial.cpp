#include <vector>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <numeric>
#include <random>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>
#include "myio_serial.h"

const int n_neigh = 4; // number of neighbors always 4 due to 2d and pbcs
const double d = 1.0;

void create_system(std::vector<vec2D> &pos, std::vector<vec2D> &vel, int nuc, double vsc, int seed) {
    /* nuc is the number of atoms in one dimension of the square lattice.
       pos and vel have been allocated already for the number of atoms in the system */
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> fdistrib(-0.5, 0.5);

    double tot_xvel_sum = 0.0; double tot_yvel_sum = 0.0; // for calculating the total v and subtracting the c.o.m velocity

    for (int i = 0; i < pos.size(); i++) {
        // lattice constant of 1 is used in the calculations, but its true value is factored in the final output
        pos[i].x = static_cast<double>(i % nuc);
        pos[i].y = static_cast<double>(i / nuc);

        vel[i].x = vsc * fdistrib(gen);
        vel[i].y = vsc * fdistrib(gen);

        tot_xvel_sum += vel[i].x;  // sum the components of velocities 
        tot_yvel_sum += vel[i].y;
    }

    tot_xvel_sum /= static_cast<double>(vel.size());  
    tot_yvel_sum /= static_cast<double>(vel.size());
    
    // subtract c.o.m. velocity components to eliminate bulk momentum
    for (int i = 0; i < vel.size(); i++) {
        vel[i].x -= tot_xvel_sum;
        vel[i].y -= tot_yvel_sum;
    }
}

void get_neigh_list(std::vector<std::array<int, n_neigh>>& nlist, int nat, int nuc) {
    // Assuming neigh is already sized appropriately: neigh[nat][n_neigh]
    for (int i = 0; i < nat; ++i) {
        int row = (i / nuc);
        int col = (i % nuc);

        int leftCol = (col < 1) ? nuc - 1 : col - 1;
        int rightCol = (col >= nuc - 1) ? 0 : col + 1;
        int lowerRow = (row < 1) ? nuc - 1 : row - 1;
        int upperRow = (row >= nuc - 1) ? 0 : row + 1;

        nlist[i][0] = leftCol + row * nuc;
        nlist[i][1] = rightCol + row * nuc;
        nlist[i][2] = col + lowerRow * nuc;
        nlist[i][3] = col + upperRow * nuc;
    }
}


void accel(std::vector<vec2D> &pos, int i, double *u, vec2D &a, std::vector<std::array<int, n_neigh>> &neigh, double box, double k2, double k3, double d) {
    // Calculate the potential energy u and acceleration a of atom whose index is i.
    double u2 = 0.0; double u3 = 0.0;
    a.x = 0.0; a.y = 0.0;

    *u = 0.0;   // reset the potential energy for a new calculation
    for (int j = 0; j < n_neigh; j++) {
        int in = neigh[i][j];

        double dx = pos[in].x - pos[i].x;
        double dy = pos[in].y - pos[i].y;

        // apply periodic boundary conditions
        if (dx < -0.5 * box) dx += box; 
        if (dx >= 0.5 * box) dx -= box;
        if (dy < -0.5 * box) dy += box;
        if (dy >= 0.5 * box) dy -= box;

        double r = std::sqrt(dx*dx + dy*dy);
        // Calculate the potential energy
        double c = r-d;
        u2 += 0.5 * k2 * c*c;
        u3 += 1.0 / 3.0 * k3 * c*c*c;  // 1.0 / 3.0

        double dxr_term = dx/r;    
        double dyr_term = dy/r;

        double fx2 = k2 * c * dxr_term;
        double fx3 = k3 * c*c * dxr_term;

        double fy2 = k2 * c * dyr_term;
        double fy3 = k3 * c*c * dyr_term;

        // update acceleration
        a.x += fx2 + fx3;
        a.y += fy2 + fy3;
    }
    *u = 0.5 * (u2 + u3); // eliminate double counting by dividing by two
}

void handle_pbcs(std::vector<vec2D> &pos, int i, double box) {
    /* Apply periodic boundaries */
    if (pos[i].x < 0.0) pos[i].x += box;
    if (pos[i].y < 0.0) pos[i].y += box;
    if (pos[i].x >= box) pos[i].x -= box;
    if (pos[i].y >= box) pos[i].y -= box;
}

void leapfrog(std::vector<vec2D> &pos, std::vector<vec2D> &vel, std::vector<vec2D> &vel0, 
              std::vector<vec2D> &a, std::vector<double> &ep, std::vector<double> &ek, int nat, 
              double dt, double box) {
    /* Perform Leapfrog time itegration by one timestep. Also handle periodic boundaries 
       and update kinetic and potential energies */
    for (int i = 0; i < nat; i++) {
        vel[i].x += dt * a[i].x;
        vel[i].y += dt * a[i].y;

        pos[i].x += dt * vel[i].x;
        pos[i].y += dt * vel[i].y;

        handle_pbcs(pos, i, box);

        double vx = 0.5 * (vel0[i].x + vel[i].x);
        double vy = 0.5 * (vel0[i].y + vel[i].y);
        ek[i] = 0.5 * (vx*vx + vy*vy);
    }

}

int main(int argc, char** argv) {
    int seed = 12345;
    int nuc = atoi(*++argv);
    double dt = atof(*++argv);
    int maxt = atoi(*++argv);
    double vsc = atof(*++argv);
    
    int thermorate = (argc > 5) ? atoi(*++argv) : thermorate;
    int dumprate = (argc > 6) ? atoi(*++argv) : dumprate;  

    double k2 = 1.0; double k3 = 0.1;   // default values
    k2 = (argc > 7) ? atof(*++argv) : k2;
    k3 = (argc > 8) ? atof(*++argv) : k3;  

    std::string dumpname = "mdrun.dump";

    int nat = nuc*nuc;
    double box = static_cast<double>(nuc);

    std::vector<vec2D> pos(nat, vec2D(0.0, 0.0)), vel(nat, vec2D(0.0, 0.0)), vel0(nat, vec2D(0.0, 0.0)), a(nat, vec2D(0.0, 0.0));
    std::vector<double> ep(nat), ek(nat);
    std::vector<std::array<int, n_neigh>> nlist(nat);

    create_system(pos, vel, nuc, vsc, seed);    // initializes positions and velocities of the atoms
    get_neigh_list(nlist, nat, nuc);   // set up the neighbor lists

    printcoords("frame0.dump", pos, box, ep, ek, 0.0);
    // Everything is ready for the time integration

    auto t0 = std::chrono::high_resolution_clock::now();

    for (int t = 0; t < maxt; t++) {    // time integration
        vel0 = vel;
        for (int i = 0; i < nat; i++) {
            accel(pos, i, &ep[i], a[i], nlist, box, k2, k3, d);
        }
        leapfrog(pos, vel, vel0, a, ep, ek, nat, dt, box);

        if (thermorate > 0 && t%thermorate == 0) {
            print_thermodata(ep, ek, dt, t);
        } 
        if (dumprate > 0) {
            if (t%dumprate == 0) {
                printcoords(dumpname, pos, box, ep, ek, t);
            }
        }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> runtime = t1 - t0;
    std::cout << "Elapsed time: " << runtime.count() << " s\n";
    return 0;
}













//double x0_left_neigh, double xn_right_neigh, double x0_top_neigh, double xn__neigh
