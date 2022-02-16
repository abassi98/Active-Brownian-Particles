#include "header.h"

using namespace std;

double point::distance_to_point(point &_r){
    return sqrt((x-_r.x)*(x-_r.x) + (y-_r.y)*(y-_r.y));
}


point getVector(const point&A, const point& B){
  // Get the vector  pointing from B to A
  point C;
  C.x = A.x - B.x;
  C.y = A.y - B.y;
  return C;
}


ABP_2d::ABP_2d(region &_reactant, double &_dt, double &_v, double &_D_r, double &_D_theta, double &_k, double &_L, double & _mu, double &_w){
    /**
     * @brief Constructor: initialize the position and the orientation of the particle
     * make sure to call it before any dynamics step
     * 
     * @param reactant is the initial region, initital position randomly drawn inside it
     * @param init_theta is the initial orientation
     * 
     */

    // Random generator 
    uniform_real_distribution<double> uniform_r(0.0, 1.0);
    uniform_real_distribution<double> uniform_theta(0.0, 2*M_PI);

    // If not empty clear the vectors
    positions.clear();
    thetas.clear();

    // Set reactant region
    reactant = _reactant;

    // Open debug
    debug_file = "debug.txt";
    ofstream debug(debug_file);

    // Generate starting point inside reactant and starting angle 
    point start_point;
    double start_theta = uniform_theta(engine);
    start_point.x = _reactant.x - 0.02;
    start_point.y = _reactant.y;
    apply_pbc_to_point(start_point); // pbc
    
   
    
    // Add the first element
    positions.push_back(start_point);
    thetas.push_back(start_theta);

    

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

void ABP_2d::apply_pbc_to_point(point &A){
    if(A.x > L/8) {A.x -= L/4;}
    if(A.x < -L/8) {A.x += L/4;}
    if(A.y > L/8) {A.y -= L/4;}
    if(A.y < -L/8) {A.y += L/4;}
}

void ABP_2d::apply_pbc(){
    /**
     * @brief Apply periodic boundary conditions and update cross parameters
     * 
     */
    point &last = positions.back();
    if(last.x > L/8.) {
        last.x -= L/4.;}
    if(last.x < -L/8.) {
        last.x += L/4.;}
    if(last.y > L/8.) {
        last.y -= L/4.;}
    if(last.y < -L/8.){
        last.y += L/4.;}
}

double ABP_2d::pbc_distance(const point&A, const point &B){
    /**
     * @brief Compute the distance between points with periodic boundary conditions
     * 
     */
    point C = getVector(A,B);
    apply_pbc_to_point(C);
    
    return sqrt(C.x*C.x + C.y*C.y);
}

double ABP_2d::potential(const point &_r){
    /**
     * @brief Compute the potential at position r, define to have minimum at (0,0)
     * 
     */
    return k*(sin(8*M_PI*(_r.x+3./16.*L)/L) + sin(8*M_PI*(_r.y + 3./16.*L)/L));
}



point ABP_2d::compute_force(){
    /**
     * @brief Compute force acting on a particle due to potential, take its last position
     * To called after there is at least one element in positions
     * 
     */
    point force;

    // Compute force
    force.x = - 8.0*M_PI/L*k*cos(8*M_PI*(positions.back().x + 3./16.*L)/L);
    force.y = - 8.0*M_PI/L*k*cos(8*M_PI*(positions.back().y + 3./16.*L)/L);

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



bool ABP_2d::is_inside_region(const region &target){
    /**
     * @brief Check if the last point of the trajectories lies inside a region with pbc
     */
    point &A = positions.back();
    bool is_inside = false;

    double _r = pbc_distance(A, target);
    if(_r<target.radius){
        is_inside = true;
    }
    return is_inside;
}


void ABP_2d::print_bool_dynamics(const region &target, unsigned &max_num_steps, string &filename){
    /**
     * @brief Run dynamics for many steps and compute statistics
     * 
     * @param target is the target region to hit
     * @param max_num_steps is the maximum number of time steps allowed to run a particle dynamics
     * 
     * @return print on a file a unsigned sequence, 0 if the particle is inside the reactant, 
     * 1 if it is outside of everything, 2 if it is inside target and
     * return error if it is in both (impossible)
     */
    

    // Random generator
    normal_distribution<double> normal_x;
    normal_distribution<double> normal_y;
    normal_distribution<double> normal_theta;

    // Output
    ofstream out(filename);

    // Run a super dynamics of max_num_steps steps
    for (unsigned step=0; step< max_num_steps; ++step){

        // Update is inside reactant and target bools
        bool is_inside_reactant = is_inside_region(reactant);
        bool is_inside_target = is_inside_region(target);

        // Generate white gaussian noise
        double noise_x = normal_x(engine);
        double noise_y = normal_y(engine);
        double noise_theta = normal_theta(engine);

        // Print on file
        out<<is_inside_reactant<<" "<<is_inside_target<<endl;

        // Dynamics steps
        position_step(noise_x, noise_y); // Update the position, appending the new posistion to the queu of the positions vector
        theta_step(noise_theta); // Update the angle theta, appending the new angle to the queu of the positions vector
    

        // Apply periodic boundary conditions
        apply_pbc();
        
    }
    out.close();
}


