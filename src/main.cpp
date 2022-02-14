#include "header.h"

using namespace std;

int main()
{
    // Inizialize random generator with seed to have repricability
    default_random_engine engine;
    engine.seed(0);



    // Set parameters
    //unsigned _N_steps = 4000;
    double _dt = 0.0001;
    double _v =2.0;
    double _D_r = 1.0;
    double _D_theta = 1.0;
    double _k = 1.0;
    double _L = 1.0;
    double _mu = 0.5;
    double _w = 0.0;
    

    // Define reactant region
    region reactant;
    reactant.x = 3./16.*_L;
    reactant.y =  3./16.*_L;
    reactant.radius = _L/16.;

    //Define target region
    region target;
    target.x = reactant.x ;
    target.y = reactant.y + _L/4.;
    target.radius = _L/16.;


    // Maximum number of steps to stop and number of particles 
    unsigned max_num_steps = 100000;
    unsigned num_particles = 100;
    
    // Dynamics
    double mean_number_steps;
    mean_number_steps = mean_search_steps(reactant, target, num_particles, max_num_steps, _dt, _v, _D_r, _D_theta, _k, _L, _mu, _w);

    cout<<"Mean number of steps: "<<mean_number_steps<<endl;
    
    return  0;
}



