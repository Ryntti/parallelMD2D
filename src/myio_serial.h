#include <string>
#include <vector>

struct vec2D {
    double x,y;
    vec2D(double x, double y) : x{x}, y{y} {}
};

void printcoords(std::string dumpfname, std::vector<vec2D> &pos, double box, std::vector<double> &ep, std::vector<double> &ek, double t);

void print_thermodata(std::vector<double> &ep, std::vector<double> &ek, double dt, double t);