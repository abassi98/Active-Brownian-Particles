#include "header.h"

using namespace std;

int main()
{
    // Inizialize random seed to have repricability
    srand(0);

    // Set parameters
    //unsigned _N_steps = 4000;
    double _dt = 0.0001;
    double _v =2.0;
    double _D_r = 1.0;
    double _D_theta = 1.0;
    double _k = 1.0;
    double _L = 1.0;
    double _mu = 2;
    double _w = 0.0;

    // Set initial position to stay in a minimum and initial direction
    point _r;
    double _theta = 0.0;
    _r.x = 3./16.*_L;
    _r.y = 3./16.*_L;

    // Initialize test active brownian particle
    ABP_2d test(_r, _theta, _dt, _v, _D_r, _D_theta, _k, _L, _mu,_w);


    //Sear target dynamics
    region target;
    target.x = _r.x + test.L/4;
    target.y = _r.y + test.L/4;
    target.radius = test.L/16.;
    unsigned max_num_steps = 10000;

    // Target search of ten particles
    unsigned num_particles = 10;
    for (unsigned i=0; i<num_particles; ++i){
        ABP_2d particle(_r, _theta, _dt, _v, _D_r, _D_theta, _k, _L, _mu,_w);
        particle.search_target(target, max_num_steps);
        // Print on file 
        particle.print_dynamics("dynamics.txt");
    }
   

    return  0;
}



