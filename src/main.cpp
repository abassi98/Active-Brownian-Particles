#include "header.h"

using namespace std;

int main()
{
    
    // Inizialize random generator with seed to have repricability
    default_random_engine engine;
    engine.seed(0);

    cout<<"Insert Peclet number"<<endl;
    double pe;
    cin>>pe;

    cout<<"Insert second parameter"<<endl;
    double l_star;
    cin>>l_star;


    // Set fix constants
    double _dt = 0.0001;
    double _D_r = 1.0;
    double _k = 1.0;
    double _L = 1.0;
    double _mu = 1.0; 
    double _w = 0.; 

    // Set other in function of dimensionless numbers
    double v_max = 8*M_PI*_k*_mu/_L;
    double _v = pe*v_max;
    double _D_theta = pe*v_max*v_max/_D_r/l_star;  
    _D_theta = 0.0;

    

    // Define reactant region
    region reactant;
    reactant.x = 0.0 ;
    reactant.y = 0.0;
    reactant.radius = _L/16.;
    
    //Define target region
    region target = reactant;


    // Dynamics
    unsigned num_steps = 1000000;
    ABP_2d test(reactant, _dt,_v,_D_r,_D_theta,_k,_L,_mu,_w);
    test.thetas[0] = 0.0;
    test.dynamics(target, num_steps);

    //cout<<"Initial theta: "<<test.thetas[0]*180/M_PI<<endl;
    
    
    return  0;
}



