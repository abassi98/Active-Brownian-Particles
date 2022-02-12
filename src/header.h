#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <complex>
#include <cstdlib>

using namespace std;

struct vec2d {
    double x = 0;
    double y = 0;
};


class ABP_2d
{
    public:
    
    // Dynamical variables
    vector<vec2d> positions;
    vector<double> thetas;

    // Coefficients
    unsigned N_steps; 
    double dt;
    double v;
    double D_r;
    double D_theta;
    double k;
    double L;
    double mu;
    double w;

    ABP_2d();

    void __init__(vec2d, double, unsigned, double, double, double, double, double, double, double, double );

    vec2d compute_force();
    void position_step(double, double);
    void theta_step(double);
   
    void dynamics();

    void print_dynamics(string);

};



