#include "header.h"

using namespace std;

int main()
{
    // Dynamical variables
    ABP_2d test;
    vec2d init_pos;
    double init_theta = 0.0;

    // Inizialize random seed to have repricability
    srand(0);

    // Set parameters
    unsigned _N_steps = 2000;
    double _dt = 0.0001;
    double _v = 1.0;
    double _D_r = 1.0;
    double _D_theta = 1.0;
    double _k = 1.0;
    double _L = 1.0;
    double _mu = 0.5;
    double _w = 0.0;

    // Initialize test active brownian particle
    test.__init__(init_pos, init_theta, _N_steps, _dt, _v, _D_r, _D_theta, _k, _L, _mu,_w);

    // Perform dynamics
    test.dynamics();

    // Print on file 
    test.print_dynamics("dynamics.txt");

    return  0;
}



