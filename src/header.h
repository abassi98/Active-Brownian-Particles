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

 class point {
    public:

    double x = 0;
    double y = 0;

    double distance(point&);
};

class region : public point {
    public:
    double radius = 1;
};

bool is_inside_target(point&, region&);

class ABP_2d
{
    public:
    
    // Dynamical variables
    vector<point> positions;
    vector<double> thetas;

    // Coefficients 
    double dt;
    double v;
    double D_r;
    double D_theta;
    double k;
    double L;
    double mu;
    double w;

    ABP_2d(point&, double&, double&, double&, double&, double&, double&, double&, double&, double& );
    ~ABP_2d();

    double potential(point&);
    point compute_force();
    void position_step(double&, double&);
    void theta_step(double&);
   
    void dynamics(unsigned&);

    
    void print_dynamics(string&);

   
    bool is_near_minimum(point&);

    unsigned search_target(region&, unsigned&, default_random_engine&);

};

double mean_search_steps(region&, region&, unsigned&, unsigned&, double&, double&, double&, double&, double&, double&, double&, double& );



