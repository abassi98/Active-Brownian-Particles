#include "header.h"

using namespace std;

ABP_2d::ABP_2d(){
    /**
     * @brief Default constructor, all variables set to zero
     * 
     */
    dt = 0;
    v = 0;
    D_r = 0;
    D_theta = 0;
    k = 0;
    L = 0;
    w = 0;
}

void ABP_2d::__init__(vec2d init_position, double init_theta,  unsigned _N_steps, double _dt, double _v, double _D_r, double _D_theta, double _k, double _L, double  _mu, double _w){
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

    // Set parameters
    N_steps = _N_steps;
    dt = _dt;
    v = _v;
    D_r = _D_r;
    D_theta = _D_theta;
    k = _k;
    L = _L;
    mu = _mu;
    w = _w;
}

void ABP_2d::compute_force(){
    vec2d next_force;

    // Compute force
    next_force.x = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().x/L);
    next_force.y = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().y/L);

    forces.push_back(next_force);

}

void ABP_2d::position_step(double noise_x, double noise_y){

    vec2d next_position;

    // Compute next position
    next_position.x = positions.back().x + v*cos(thetas.back())*dt + sqrt(2*D_r*dt)*noise_x + mu*forces.back().x*dt;
    next_position.y = positions.back().y + v*cos(thetas.back())*dt + sqrt(2*D_r*dt)*noise_y + mu*forces.back().y*dt;

    positions.push_back(next_position);
}

void ABP_2d::theta_step(double noise_theta){    
    double next_orientation = thetas.back() + w*dt + sqrt(2*D_theta*dt)*noise_theta;
    thetas.push_back(next_orientation);
}



void ABP_2d::dynamics(){
    /**
     * @brief Perform stochastic dynamics
     * 
     */

    // Random generator
    default_random_engine engine;
    normal_distribution<double> normal_x;
    normal_distribution<double> normal_y;
    normal_distribution<double> normal_theta;
     

    for(unsigned step=0; step<N_steps; ++step){
        // Generate white gaussian noise
        double noise_x = normal_x(engine);
        double noise_y = normal_y(engine);
        double noise_theta = normal_theta(engine);

        // Dynamics steps
        compute_force(); // Compute force acting on the particle due to potential, take last position which is the actual one
        position_step(noise_x, noise_y); // Update the position, appending the new posistion to the queu of the positions vector
        theta_step(noise_theta); // Update the angle theta, appending the new angle to the queu of the positions vector
    }
}

void ABP_2d::print_dynamics(string filename){
    ofstream out(filename);
    for(unsigned i=0; i<positions.size(); ++i){
        out<<positions[i].x<<" "<<positions[i].y<<" "<<thetas[i]<<" "<<forces[i].x<<" "<<forces[i].y<<endl;
    }
    out.close();
}



