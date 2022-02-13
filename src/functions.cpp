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

    // Print parameters on a file for plotting
    ofstream par("parameters.txt");
    par<<"N_steps"<<" "<<N_steps<<endl;
    par<<"dt"<<" "<<dt<<endl;
    par<<"v"<<" "<<v<<endl;
    par<<"D_r"<<" "<<D_r<<endl;
    par<<"D_theta"<<" "<<D_theta<<endl;
    par<<"k"<<" "<<k<<endl;
    par<<"L"<<" "<<L<<endl;
    par<<"mu"<<" "<<mu<<endl;
    par<<"w"<<" "<<w<<endl;
    par.close();
}

vec2d ABP_2d::compute_force(){
    /**
     * @brief Compute force acting on a particle due to potential, take its last position
     * To called after there is at least one element in positions
     * 
     */
    vec2d force;

    // Compute force
    force.x = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().x/L);
    force.y = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().y/L);

    return force;

}

void ABP_2d::position_step(double noise_x, double noise_y){
    /**
     * @brief Copmute next position and append to vector positions
     * To be called after __init__()
     * 
     */

    vec2d next_position;
    vec2d force;

    // Compute force acting on last position
    force = compute_force();

    // Compute next position
    next_position.x = positions.back().x + v*cos(thetas.back())*dt + sqrt(2*D_r*dt)*noise_x + mu*force.x*dt;
    next_position.y = positions.back().y + v*sin(thetas.back())*dt + sqrt(2*D_r*dt)*noise_y + mu*force.y*dt;

    // Append to vector positions
    positions.push_back(next_position);
}

void ABP_2d::theta_step(double noise_theta){   
    /**
     * @brief Compute next orientation
     * 
     */
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
     

    for(unsigned step=0; step<N_steps-1; ++step){
        // Generate white gaussian noise
        double noise_x = normal_x(engine);
        double noise_y = normal_y(engine);
        double noise_theta = normal_theta(engine);

        // Dynamics steps
        position_step(noise_x, noise_y); // Update the position, appending the new posistion to the queu of the positions vector
        theta_step(noise_theta); // Update the angle theta, appending the new angle to the queu of the positions vector
    }
}

void ABP_2d::print_dynamics(string filename){
    ofstream out(filename);
    for(unsigned i=0; i<positions.size(); ++i){
        out<<positions[i].x<<" "<<positions[i].y<<" "<<thetas[i]<<endl;
    }
    out.close();
}




