#include "header.h"

using namespace std;

int main()
{
    
    // Inizialize random generator with seed to have repricability
    default_random_engine engine;
    engine.seed(0);



    // Set parameters
    double _dt = 0.0001;
    double _v = 0.0;
    double _D_r = 1.0;
    double _D_theta = 0.0;
    double _k = 1.0;
    double _L = 1.0;
    double _mu = 2.0; 
    double _w = 0.; 
    
    

    // Define reactant region
    region reactant;
    reactant.x = 0.0 ;
    reactant.y = 0.0;
    reactant.radius = _L/16.;
    
    //Define target region
    region target;
    int m = 1;
    int n = 0;
    target.x = m*_L/4;
    target.y = n*_L/4.;
    target.radius = _L/16.;
    

    // Check if they overlap
    if(target.distance_to_point(reactant)<target.radius+reactant.radius){
        cout<<"Error: target and reactant overalp"<<endl;
        abort();
    }
    
    // Get quadrants and traslate coordinates
    reactant.get_quadrant(_L); 
    reactant.translate_to_origin(_L);
    target.get_quadrant(_L);
    target.translate_to_origin(_L);




    // Dynamics
    unsigned num_steps = 1000000;
    string dyn = "dynamics.txt";
    string bool_dyn =  "bool_dynamics.txt";

    ABP_2d test(reactant, _dt,_v,_D_r,_D_theta,_k,_L,_mu,_w);
    test.thetas[0] = 0.0;
    test.print_bool_dynamics(target, num_steps,bool_dyn);
    test.print_dynamics(dyn);
    //cout<<"Initial theta: "<<test.thetas[0]*180/M_PI<<endl;
    
    
    return  0;
}



