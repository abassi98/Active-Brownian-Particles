#include "header.h"

using namespace std;

ABP_2d::ABP_2d(){
    dt = 0;
    v = 0;
    D_r = 0;
    D_theta = 0;
    k = 0;
    L = 0;
    w = 0;
}

void ABP_2d::__init__(vec2d init_position, double init_theta){
    /**
     * @brief Initialize the position and the orientation of the particle
     * make sure to call it before any dynamics step
     * 
     * @param init_position is the initial position
     * @param init_theta is the initial orientation
     * 
     */

    // If not empty clear the vectors
    positions.clear();
    thetas.clear();

    // Add the first element
    positions.push_back(init_position);
    thetas.push_back(init_theta);
}


void ABP_2d::theta_step(){
    // Random generator
    default_random_engine engine;
    normal_distribution<double> normal; 

    double noise = normal(engine);
    double next_orientation = thetas.back() + w*dt + sqrt(2*D_theta*dt)*noise;
    thetas.push_back(next_orientation);
}

void ABP_2d::compute_force(){
    vec2d next_force;

    // Compute force
    next_force.x = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().x/L);
    next_force.y = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().y/L);

    forces.push_back(next_force);

}

void ABP_2d::position_step(){

    // Random generator
    default_random_engine engine;
    normal_distribution<double> normal; 
    vec2d next_position;

    double noise_x = normal(engine);
    double noise_y = normal(engine);
    

    // Compute next position
    next_position.x = positions.back().x + v*cos(thetas.back())*dt + sqrt(2*D_r*dt)*noise_x;
    next_position.y = positions.back().y + v*cos(thetas.back())*dt + sqrt(2*D_r*dt)*noise_y;

    positions.push_back(next_position);
}

void ABP_2d::print_dynamics(string filename){

    ofstream out(filename);
    for(unsigned i=0; i<positions.size(); ++i){
        out<<positions[i].x<<" "<<positions[i].y<<" "<<thetas[i]<<" "<<forces[i].x<<" "<<forces[i].y<<endl;
    }
    out.close();
}




