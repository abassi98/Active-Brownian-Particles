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

class point{
    public:
    double x = 0;
    double y = 0;

    double distance_to_point(point&);
};

class region : public point  {
    public:
    double radius = 1;
};


point getVector(const point&, const point &);


class ABP_2d
{
    public:
    
    // Dynamical variables
    vector<point> positions;
    vector<double> thetas;

    // Debug file 
    ofstream debug;

    // Coefficients 
    double dt;
    double v;
    double D_r;
    double D_theta;
    double k;
    double L;
    double mu;
    double w;

    // Reactant region
    region reactant;

    // Random generator
    default_random_engine engine;

    ABP_2d(region&, double&, double&, double&, double&, double&, double&, double&, double& );
    ~ABP_2d();

   
    void apply_pbc_to_point(point&);
    void apply_pbc();
    double pbc_distance(const point&, const point&);

    double potential(const point&);
    point compute_force();
    void position_step(double&, double&);
    void theta_step(double&);

    
    void print_dynamics(string&);

   
    bool is_near_minimum(point&);

    bool is_inside_region(const region&);
    void print_bool_dynamics(const region&, unsigned&, string&);

};




