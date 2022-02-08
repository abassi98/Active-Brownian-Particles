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
    vector<vec2d> forces;
    vector<double> thetas;

    // Coefficients
    double dt;
    double v;
    double D_r;
    double D_theta;
    double k;
    double L;
    double w;

    ABP_2d();
    ~ABP_2d();

    void __init__(vec2d, double);
    void theta_step();
    void compute_force();
    void position_step();
    void print_dynamics(string);

};



