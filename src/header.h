#include <iostream>
#include <stdexcept>
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
        double x;
        double y;
        point();
        point(double _x, double _y);
        double distance_to_point(const point&);
};

point getVector(const point&,const point &);

class region : public point  {
    public:
        double radius;
        region();
        region(double, double, double);
};





class ABP_2d
{
    public:
    
        // Dynamical variables
        vector<double> position_x;
        vector<double> position_y;
        vector<double> theta;
        vector<bool> bool_reactant;
        vector<bool> bool_target;
        vector<bool> reactive_path;
        vector<bool> transition_path;

        // Coefficients 
        unsigned num_steps;
        double dt;
        double v;
        double D_r;
        double D_theta;
        double k;
        double L;
        double mu;
        double w;

        // Dynamics time statistics
        unsigned step;
        double time;

        // Reactant region
        region reactant;
        region target;

        // Random generator
        default_random_engine engine;

        ABP_2d(const region&, const region&, unsigned&, double&, double&, double&, double&, double&, double&, double&, double& );

    
        void apply_pbc(point&);
        void apply_pbc_to_theta(double&);
        double pbc_distance(const point&, const point&);

        double potential(const point&);
        point compute_force(const point&);
        void position_step(point&, const double &, const double&, const double&);
        void theta_step(double&, double&);
        bool is_near_minimum(point&);

        bool is_inside_region(const point&, const region&);
        void dynamics(bool, bool, bool, bool);
        void print_dynamics(string&);
        void print_bool_dynamics(string&);

};




