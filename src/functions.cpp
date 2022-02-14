#include "header.h"

using namespace std;

double point::distance(point &_r){
    return sqrt((x-_r.x)*(x-_r.x) + (y-_r.y)*(y-_r.y));
}

bool is_inside_target(point &A, region &target){
    /**
     * @brief Check if a point A is inside region target
     */
    bool is_near = false;
    if (A.distance(target) <  target.radius){
        is_near = true;
    }

    return is_near;
}

ABP_2d::ABP_2d(point &init_position, double &init_theta, double &_dt, double &_v, double &_D_r, double &_D_theta, double &_k, double &_L, double & _mu, double &_w){
    /**
     * @brief Constructor: initialize the position and the orientation of the particle
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

ABP_2d::~ABP_2d(){
    /**
     * @brief Destuctor, clear position and theta vectors
     * 
     */
    positions.clear();
    thetas.clear();
}

double ABP_2d::potential(point &_r){
    /**
     * @brief Compute the potential at position r
     * 
     */
    return k*(sin(8*M_PI*_r.x/L) + sin(8*M_PI*_r.y/L));
}

point ABP_2d::compute_force(){
    /**
     * @brief Compute force acting on a particle due to potential, take its last position
     * To called after there is at least one element in positions
     * 
     */
    point force;

    // Compute force
    force.x = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().x/L);
    force.y = - 8.0*M_PI/L*k*cos(8*M_PI*positions.back().y/L);

    return force;

}

void ABP_2d::position_step(double &noise_x, double &noise_y){
    /**
     * @brief Copmute next position and append to vector positions
     * To be called after __init__()
     * 
     */

    point next_position;
    point force;

    // Compute force acting on last position
    force = compute_force();

    // Compute next position
    next_position.x = positions.back().x + v*cos(thetas.back())*dt + sqrt(2*D_r*dt)*noise_x + mu*force.x*dt;
    next_position.y = positions.back().y + v*sin(thetas.back())*dt + sqrt(2*D_r*dt)*noise_y + mu*force.y*dt;

    // Append to vector positions
    positions.push_back(next_position);
}

void ABP_2d::theta_step(double &noise_theta){   
    /**
     * @brief Compute next orientation
     * 
     */
    double next_orientation = thetas.back() + w*dt + sqrt(2*D_theta*dt)*noise_theta;
    thetas.push_back(next_orientation);
}



void ABP_2d::dynamics(unsigned &N_steps){
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

void ABP_2d::print_dynamics(string &filename){
    ofstream out(filename);
    for(unsigned i=0; i<positions.size(); ++i){
        out<<positions[i].x<<" "<<positions[i].y<<" "<<thetas[i]<<endl;
    }
    out.close();
}

bool ABP_2d::is_near_minimum(point &_r){
    /**
     * @brief Check if the position is inside a minimum region
     * Defined to have potential energy <=-1.5
     * 
     */
    bool is_near = false;
    double u = potential(_r);
    if (u<-1.5){
        is_near = true;
    }

    return is_near;
}


unsigned ABP_2d::search_target(region &target, unsigned &max_num_steps, default_random_engine &engine){
    /**
     * @brief Run dynamics until the particle hit (is inside) the target region
     * 
     * @return Number of steps taken to hit the target if they are less than max_num_steps, max_num_steps otherwise
     */
    // Random generator
    normal_distribution<double> normal_x;
    normal_distribution<double> normal_y;
    normal_distribution<double> normal_theta;

    // Stop condition
    bool is_inside = false;

    // Step counter
    unsigned step = 0;

    while (is_inside==false && step < max_num_steps){
        // Generate white gaussian noise
        double noise_x = normal_x(engine);
        double noise_y = normal_y(engine);
        double noise_theta = normal_theta(engine);

        // Dynamics steps
        position_step(noise_x, noise_y); // Update the position, appending the new posistion to the queu of the positions vector
        theta_step(noise_theta); // Update the angle theta, appending the new angle to the queu of the positions vector

        // Update step counter and stop condition
        is_inside = is_inside_target(positions.back(), target);
        ++step;
    }

    return step;
}


double mean_search_steps(region &reactant, region &target, unsigned &num_particles, unsigned &max_num_steps, double &_dt, double &_v, double &_D_r, double &_D_theta, double &_k, double &_L, double & _mu, double &_w){
    /**
     * @brief Compute the average number of steps for particle to jump from the reactant region to target region
     * where the search is done for a max_num_steps
     * 
     * @param reactant is the starting regione, where the starting point is drawn uniformly
     * @param target is the target region for search
     * @param num_particles is the total number of trial particles
     * @param max_num_steps is the maximum number of search steps for each particle
     * 
     * @return the average number of steps taken by the particle to hit the target starting from rectant region
     * 
     */

    // Random generator 
    default_random_engine engine;
    uniform_real_distribution<double> uniform_r(0.0, sqrt(reactant.radius));
    uniform_real_distribution<double> uniform_theta(0.0, 2*M_PI);

    // Generate starting point inside reactant and starting angle 
    point start_point;
    double _r = uniform_r(engine);
    double _theta = uniform_theta(engine);
    double start_theta = uniform_theta(engine);
    start_point.x = reactant.x + _r*cos(_theta);
    start_point.y = reactant.y + _r*sin(_theta);
    
    // Mean number of steps
    double mean_num_steps = 0.0;
    unsigned print_step = num_particles/10;

    for (unsigned i=0; i<num_particles; ++i){
        ABP_2d particle(start_point, start_theta, _dt, _v, _D_r, _D_theta, _k, _L, _mu,_w);
        unsigned num_steps = particle.search_target(target, max_num_steps, engine);
        mean_num_steps += (double) num_steps;

        // Print dynamics
        if( i%print_step ==0 ){
            string filename = "dynamics"+to_string(i)+".txt";
            particle.print_dynamics(filename);      
        }
        
    }
    return mean_num_steps/((double) num_particles);
}







