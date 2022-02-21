#include "header.h"

using namespace std;

point::point(){
    x = 0.0;
    y = 0.0;
}
point::point(double _x, double _y){
    x = _x;
    y = _y;
}

double point::distance_to_point(const point &_r){
    return sqrt((x-_r.x)*(x-_r.x) + (y-_r.y)*(y-_r.y));
}

region::region() : point() {
    radius = 0;
}
region::region(double _x, double _y, double _radius) : point(_x,_y) {
    if(_radius<0){
        throw invalid_argument("Received negative radius");
    }
    radius = _radius;
}

point getVector(const point& A, const point& B){
  /**
   * @brief Get the vector connecting B to point
   * 
   * @param B is the second point
   * 
   * @return a point connecting two points
   * 
   */
  point C;
  C.x = A.x - B.x;
  C.y = A.y - B.y;
  return C;
}


ABP_2d::ABP_2d(const region &_reactant,  const region &_target, unsigned &_num_steps, double &_dt, double &_v, double &_D_r, double &_D_theta, double &_k, double &_L, double & _mu, double &_w){
    /**
     * @brief Default constructor, accept parameters and perform validity checks
     * 
     * @param _reactant is the initial region, initital position chosen at zero, direction drawn randomly in [0,2*pi]
     * @param _target is the target region the particle aim to search
     * @param _num_steps is the total number of steps of the long dynamics
     * @param _dt is the timestep of discretized dynamics
     * @param _v is the module of the velocity
     * @param _D_r is the diffusion coefficient for updating the position of the particle
     * @param  _D_theta is the diffusion coefficient for the direction
     * @param _k is a constant measuring the strenght of the force
     * @param _L is the length scale oif the periodic potential (dynamics runs with periodic boundary conditions in the square [-L/8,L/8]**2)
     * @param _mu is a coefficient measuring how much the particle feels the potential
     * @param _w is the velocity of direction update (chirality of the particle)
     * 
     */



    // Allocate memory
    num_steps = _num_steps;
    position_x.clear();
    position_y.clear();
    theta.clear();
    bool_reactant.clear();
    bool_target.clear();
    reactive_path.clear();

    // Set reactant region
    reactant = _reactant;
    target = _target;
    // Apply pbc
    apply_pbc(reactant);
    apply_pbc(target);

    /*
    // Check if they overlap
    if(pbc_distance(reactant, target)<(reactant.radius+target.radius)){
        throw invalid_argument("Received overlapping regions");
    }*/

    // Check other parameters
    if(_dt<0){throw invalid_argument("Received negative timestep");}
    if(_v<0){throw invalid_argument("Received negative velocity");}
    if(_D_r<0){throw invalid_argument("Received negative diffusion coefficient");}
    if(_D_theta<0){throw invalid_argument("Received negative diffusion coefficient");}
    if(_L<0){throw invalid_argument("Received null volume");}


    // Generate starting point at the center of reactant
    point start_point(reactant.x,reactant.y);
    double start_theta = 0.0;
    apply_pbc(start_point); // pbc
    

    // Set initial conditions
    position_x.push_back(start_point.x);
    position_y.push_back(start_point.y);
    theta.push_back(start_theta);
    bool_reactant.push_back(is_inside_region(start_point, reactant));
    bool_target.push_back(is_inside_region(start_point, target));
    reactive_path.push_back(false);
    transition_path.push_back(false);

    // Set parameters
    dt = _dt;
    v = _v;
    D_r = _D_r;
    D_theta = _D_theta;
    k = _k;
    L = _L;
    mu = _mu;
    w = _w;

    // Dynamics statistics
    step = 0;
    time = 0.0;
    
}



void ABP_2d::apply_pbc(point &A){
    if(A.x > L/8) {A.x -= L/4;}
    if(A.x < -L/8) {A.x += L/4;}
    if(A.y > L/8) {A.y -= L/4;}
    if(A.y < -L/8) {A.y += L/4;}
}

void ABP_2d::apply_pbc_to_theta(double &theta){
    /**
     * @brief Translate theta in [0,2*pi] no matter its value
     * 
     * @param theta is the direction of the particle
     */

    // Translate in interval [0,2pi]
    if(theta>M_PI){
        theta -= 2*M_PI;
    }
    if(theta<-M_PI){
        theta += 2*M_PI;
    }
}



double ABP_2d::pbc_distance(const point&A, const point &B){
    /**
     * @brief Compute the distance between points with periodic boundary conditions
     * 
     */
    point C = getVector(A, B);
    apply_pbc(C);
    
    return sqrt(C.x*C.x + C.y*C.y);
}

double ABP_2d::potential(const point &_r){
    /**
     * @brief Compute the potential at position r, define to have minimum at (0,0)
     * 
     */
    return k*(sin(8*M_PI*(_r.x+3./16.*L)/L) + sin(8*M_PI*(_r.y + 3./16.*L)/L));
}



point ABP_2d::compute_force(const point& position){
    /**
     * @brief Compute force acting on a point due to potential
     * 
     */
    point force;

    // Compute force
    force.x = - 8.0*M_PI/L*k*cos(8*M_PI*(position.x + 3./16.*L)/L);
    force.y = - 8.0*M_PI/L*k*cos(8*M_PI*(position.y + 3./16.*L)/L);
    return force;
}

void ABP_2d::position_step(point &position, const double &theta, const double &noise_x, const double &noise_y){
    /**
     * @brief Update position according to stochastic dynamics
     * 
     * @param position is the actual position of the particle
     * @param theta is the direction of the particle velocity
     * @param noise_x is the white gaussian noise along x direction
     * @param noise_y is the white gaussian noise along y direction
     * 
     */


    // Compute force acting on last position
    point force = compute_force(position);

    // Compute next position
    position.x +=  v*cos(theta)*dt + sqrt(2.*D_r*dt)*noise_x + mu*force.x*dt;
    position.y += v*sin(theta)*dt + sqrt(2.*D_r*dt)*noise_y + mu*force.y*dt;

    // Apply periodic boundary conditinos
    apply_pbc(position);
}

void ABP_2d::theta_step(double &theta, const double &noise_theta){   
    /**
     * @brief Update orientation according to dynamics
     * 
     * @param theta is the direction of the velocity
     * @param noise_theta is the white gaussian noise for direction
     * 
     */

    // Compute next direction
    theta += w*dt + sqrt(2.*D_theta*dt)*noise_theta;

    // Apply periodic boundary conditions
    apply_pbc_to_theta(theta);
}



bool ABP_2d::is_near_minimum(point &_r){
    /**
     * @brief Check if the position is inside a minimum region
     * Defined to have potential energy <=-1.5
     * 
     * @param _r is the position to be checked
     * 
     */
    bool is_near = false;
    double u = potential(_r);
    if (u<-1.5){
        is_near = true;
    }

    return is_near;
}



bool ABP_2d::is_inside_region(const point &A,  const region &target){
    /**
     * @brief Check if the point lies inside a region with pbc
     * 
     * @param A is the position to be checked
     * @param target is the region
     * 
     * @return true if A is inside target, false otherwise
     */
    bool is_inside = false;

    double _r = pbc_distance(A, target);
    if(_r<target.radius){
        is_inside = true;
    }
    return is_inside;
}


void ABP_2d::dynamics(bool track_in_reactant=false, bool track_in_target=false, bool track_reactive_path=true, bool track_transition_path=true){
    /**
     * @brief Run dynamics for many steps and compute statistics
     *
     */
    
    // Computation clocks
    clock_t start, end;

    // Dynamical variables
    point position;
    double theta_dyn;
    bool is_inside_reactant;
    bool is_inside_target;
    bool reactive = false;
    bool transition = false;
    bool has_ever_visited_target = false;

    // Set initial values
    position.x = position_x.back();
    position.y = position_y.back();
    theta_dyn = theta.back();
    is_inside_reactant = bool_reactant.back();
    is_inside_target = bool_target.back();
    reactive = reactive_path.back();
    transition = transition_path.back();
    unsigned count_reactive = 0;

    // Random generator
    default_random_engine engine;
    normal_distribution<double> normal_x;
    normal_distribution<double> normal_y;
    normal_distribution<double> normal_theta;

    // Starting time
    start = clock();
    srand(0);
    // Run a super dynamics of max_num_steps steps
    for (unsigned i=0; i<num_steps; ++i){
        // Generate white gaussian noise
        double noise_x = normal_x(engine);
        double noise_y = normal_y(engine);
        double noise_theta = normal_theta(engine);

        // Dynamics stepsss
        position_step(position, theta_dyn, noise_x, noise_y); // Update the position, appending the new posistion to the queu of the positions vector
        theta_step(theta_dyn, noise_theta); // Update the angle theta, appending the new angle to the queu of the positions vector

        // Take next bools
        bool next_is_inside_reactant = is_inside_region(position, reactant);
        bool next_is_inside_target = is_inside_region(position, target);
        bool is_exing_reactant = is_inside_reactant && !next_is_inside_reactant;
        bool is_entering_target = !is_inside_target && next_is_inside_target;
        bool is_entering_reactant = !is_inside_reactant && next_is_inside_reactant;


        // Track reactive paths
        // Defined as those for which the particle reaches T region after exing R before visiting again R
        if(track_reactive_path){
            // Set reactive bool true when the particle is exing reactant region
            if(is_exing_reactant){
                reactive = true;}
            
            // Set reactive bool false when particle is the reactant again and set false the whole previous reactive path
            if(is_entering_reactant){
                reactive = false;
                for(unsigned i=0; i<count_reactive;++i){
                    reactive_path[reactive_path.size()-i-1] = false;
                }
                count_reactive = 0; // Set counter to zero
            }
            
             // Set reactive bool false when the particle is entering target and sety count to zero
            if(is_entering_target){
                reactive = false;
                count_reactive = 0;
            }
            
            // Update counter when reactive bool is true
            if(reactive){
                ++count_reactive;
            }
            // Append to reactive_path
            reactive_path.push_back(reactive);
        }
        
        // Track transition times
        // Defined as those for which the particle enters R after having visited T till it goes back to T
        if(track_transition_path){
            if(is_inside_target){has_ever_visited_target=true;}
            if(has_ever_visited_target){
                if(is_entering_reactant){transition=true;}
                if(is_entering_target){transition=false;}
            }
            transition_path.push_back(transition);
        }

        // Append positions and bools
        position_x.push_back(position.x);
        position_y.push_back(position.y);
        theta.push_back(theta_dyn);

        // Check track conditions and push_back
        if(track_in_reactant){bool_reactant.push_back(next_is_inside_reactant);}
        if(track_in_target){bool_target.push_back(next_is_inside_target);}

        // Set new states for bool conditions
        is_inside_reactant = next_is_inside_reactant;
        is_inside_target = next_is_inside_target;

        // Update dynamics counter
        ++step;
    }
    // End time
    end = clock();
    time += 1.0*(end-start)/CLOCKS_PER_SEC;
}


void ABP_2d::print_dynamics(string &filename){
    ofstream out(filename);
    for(unsigned i=0; i<num_steps+1; ++i){
        out<<position_x[i]<<" "<<position_y[i]<<" "<<theta[i]<<endl;
    }
    out.close();
}

void ABP_2d::print_bool_dynamics(string &filename){
    ofstream out(filename);
    for(unsigned i=0; i<num_steps+1; ++i){
        out<<bool_reactant[i]<<" "<<bool_target[i]<<endl;
    }
    out.close();
}