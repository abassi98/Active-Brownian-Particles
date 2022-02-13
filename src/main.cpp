#include "header.h"

using namespace std;

int main()
{
    // Inizialize random seed to have repricability
    srand(0);

    // Set parameters
    unsigned _N_steps = 4000;
    double _dt = 0.0001;
    double _v =2.0;
    double _D_r = 1.0;
    double _D_theta = 1.0;
    double _k = 1.0;
    double _L = 1.0;
    double _mu = 2;
    double _w = 0.0;

    // Set initial position to stay in a minimum and initial direction
    vec2d _r;
    double _theta = 0.0;
    _r.x = 3./16.*_L;
    _r.y = 3./16.*_L;

    // Initialize test active brownian particle
    ABP_2d test(_r, _theta, _N_steps, _dt, _v, _D_r, _D_theta, _k, _L, _mu,_w);


    // Perform dynamics
    test.dynamics();

    // Print on file 
    test.print_dynamics("dynamics.txt");

    return  0;
}



