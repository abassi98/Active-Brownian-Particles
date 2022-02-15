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
    double _v = 0.;
    double _D_r = 1.0;
    double _D_theta = 1.0;
    double _k = 1.0;
    double _L = 1.0;
    double _mu = 2.; 
    double _w = 0.; 
    
    

    // Define reactant region
    region reactant;
    reactant.x = 0.0;
    reactant.y = 0.0;
    reactant.radius = _L/16.;

    //Define target region
    region target;
    target.x = reactant.x;
    target.y = reactant.y - _L/4.;
    target.radius = _L/16.;
    
    // Dynamics
    unsigned num_steps = 10000;
    string dyn = "dynamics.txt";
    string bool_dyn =  "bool_dynamics.txt";

    ABP_2d test(reactant, _dt,_v,_D_r,_D_theta,_k,_L,_mu,_w);
    test.print_bool_dynamics(target, num_steps,bool_dyn);
    test.print_dynamics(dyn);
    //Crossing
    cout<<"Quadrant x: "<<test.quadrant_x<<endl;
    cout<<"Quadrant y: "<<test.quadrant_y<<endl;
    
    return  0;
}



