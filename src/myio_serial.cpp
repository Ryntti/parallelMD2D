#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include "myio_serial.h" 

void printcoords(std::string dumpfname, std::vector<vec2D> &pos, double box, std::vector<double> &ep, std::vector<double> &ek, double t) {
    int nat = pos.size();
    const double xsc = 1.0;
    std::ofstream dump;
    dump.open(dumpfname, std::ios::out | std::ios::app);
    dump << nat << "\n";
    dump << "Frame number " << t << " fs boxsize " << xsc*box << " " << xsc*box << " " << 10.0 << "\n"; 

    for (int i = 0; i < nat; i++) {
        dump << "Fe " << i + 1 << " " << std::scientific << std::setprecision(10) 
             << xsc*(pos[i].x - 0.5*box) << " " << xsc*(pos[i].y - 0.5*box) 
             << " 0.0000000000 " << ep[i] << " " << ek[i] << "\n";
    }
}

void print_thermodata(std::vector<double> &ep, std::vector<double> &ek, double dt, double t) {
    double epsum = std::accumulate(ep.begin(), ep.end(), 0.0);
    double eksum = std::accumulate(ek.begin(), ek.end(), 0.0);
    std::cout << dt * t << " " << epsum + eksum << " " << epsum << 
                 " " << eksum << std::endl;
}   