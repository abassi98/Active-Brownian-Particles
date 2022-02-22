### Import modules
import matplotlib.pyplot as plt
import numpy as np
import math
from abp import ABP_2d, point, region

### Set fixed parameters
L = 1.0
k=1.0
mu = 1.0
dt = 1e-5
D_r = 1.0
w = 0.0
v_max = 8*math.pi*k*mu/L

### Contour plot of potential
# Define potential function
def plot_potential():
    
    fig =  plt.figure(figsize=(4, 4))
    def potential(x, y):
        return k*(np.sin(8*math.pi*(x +3./16.*L)/L) + np.sin(8*math.pi*(y+3./16.*L)/L))

    # Mesh grid
    x = np.linspace(-L/8, L/8, 1000)
    y = np.linspace(-L/8, L/8, 1000)
    xx, yy = np.meshgrid(x, y)

    # Potential
    zz = potential(xx,yy)

    # Plot
    cf = plt.contourf(xx,yy, zz, cmap="RdBu_r")
    plt.colorbar()

    # Rendering
    plt.title("Potential")
    plt.xlabel("x")
    plt.ylabel("y")
    x_ticks = [-L/8, 0, L/8 ]
    x_labels = ["- L/8","0", "L/8" ]
    y_ticks = [0, L/8]
    y_labels = ["0", "L/8"]
    plt.xticks(x_ticks, x_labels)
    plt.yticks(y_ticks, y_labels)
    #plt.legend()
    plt.tight_layout()
    return fig

### Define some useful functions
# Plot dynamics examples
def plot_trajectory(pe_vec, lstar_vec, w=0.0, num_steps=10**3):
    # Define fake target and reactant regions
    reactant = region(x=0,y=0,radius=L/40.0)
    target = region(x=L/16,y=L/16,radius=L/40.0)

    fig =  plt.figure(figsize=(8, 8))
    sax = fig.add_axes([0.06, 0.06, 0.9, 0.9])

    # Super grid ticks
    s_xticks = [0.177, 0.4767, 0.7853]
    s_yticks = [0.78, 0.5, 0.2]
    sax.set_xticks(s_xticks, pe_vec)
    sax.set_yticks(s_yticks, lstar_vec)
    sax.grid(False)
    sax.set_xlabel("$Pe$")
    sax.set_ylabel("$\ell^*$", rotation = 0)


    x_ticks = [-L/8, 0, L/8 ]
    x_labels = ["- L/8","0", "L/8" ]
    y_ticks = [0, L/8]
    y_labels = ["0", "L/8"]

    for i in range(len(lstar_vec)):
        for j in range(len(pe_vec)):
            n = i*len(pe_vec)+j
            ax = fig.add_subplot(3,3,n+1)
            # Define parameters
            pe = pe_vec[j]
            l_star = lstar_vec[i]
            v = pe*v_max
            D_theta = float(pe*v_max*v_max/(1.*l_star))
            # Compute statistics
            particle = ABP_2d(reactant,target, num_steps = num_steps, dt=dt, v=v, D_r=D_r, D_theta= D_theta, k=k, L=L, mu=mu, w=w)
            particle.dynamics()
            # Retireve data
            x_dyn= np.array(particle.position_x)
            y_dyn = np.array(particle.position_y)
            
            #Plot
            ax.plot(x_dyn, y_dyn,c="black")
            ax.set_xticks(x_ticks, x_labels)
            ax.set_yticks(y_ticks, y_labels)

    #fig.tight_layout(pad=0.2)
    return fig

# Plot single passive particle dynamics
def plot_passive_particle(reactant, target,is_out=False, w=0.0, num_steps=10**7 ):
    fig =  plt.figure(figsize=(4, 4))
    x_ticks = [-L/8, 0, L/8 ]
    x_labels = ["- L/8","0", "L/8" ]
    y_ticks = [0, L/8]
    y_labels = ["0", "L/8"]


    # Define parameters]
    v = 0
    D_theta = 0
    # Compute statistics
    particle = ABP_2d(reactant,target, num_steps=num_steps, dt=dt, v=v, D_r=D_r, D_theta= D_theta, k=k, L=L, mu=mu, w=w)
    particle.dynamics()
    # Retireve data
    x_dyn= np.array(particle.position_x)
    y_dyn = np.array(particle.position_y)
    reactive_path = np.array(particle.reactive_path)
    if is_out==True:
        is_outside = np.logical_and(np.logical_not(particle.bool_reactant), np.logical_not(particle.bool_target))
        x_dyn = x_dyn[is_outside]
        y_dyn = y_dyn[is_outside]
    else:
        x_dyn = x_dyn[reactive_path]
        y_dyn = y_dyn[reactive_path]
        plt.text(target.x, target.y, "T",c = "white", size ="xx-large" , ha = "center", va="center")
    #Plot
    h = plt.hist2d(x_dyn, y_dyn, bins = 100,cmap="RdBu_r", density = True)
    plt.text(reactant.x, reactant.y, "R", c = "white", size ="xx-large", ha = "center" , va="center")
    plt.xticks(x_ticks, x_labels)
    plt.yticks(y_ticks, y_labels)
    plt.colorbar(h[3])
    return fig

# Plot transition probability density
def plot_transition_density(reactant, target, pe_vec, lstar_vec, is_out=False, w=0.0, num_steps=10**7):
    fig =  plt.figure(figsize=(8, 8))
    sax = fig.add_axes([0.06, 0.06, 0.9, 0.9])

    # Super grid ticks
    s_xticks = [0.177, 0.4767, 0.7853]
    s_yticks = [0.78, 0.5, 0.2]
    sax.set_xticks(s_xticks, pe_vec)
    sax.set_yticks(s_yticks, lstar_vec)
    sax.grid(False)
    sax.set_xlabel("$Pe$")
    sax.set_ylabel("$\ell^*$", rotation = 0)


    x_ticks = [-L/8, 0, L/8 ]
    x_labels = ["- L/8","0", "L/8" ]
    y_ticks = [0, L/8]
    y_labels = ["0", "L/8"]

    for i in range(len(lstar_vec)):
        for j in range(len(pe_vec)):
            n = i*len(pe_vec)+j
            ax = fig.add_subplot(3,3,n+1)
            # Define parameters
            pe = pe_vec[j]
            l_star = lstar_vec[i]
            v = pe*v_max
            D_theta = float(pe*v_max*v_max/(1.*l_star))
            # Compute statistics
            particle = ABP_2d(reactant,target, num_steps = num_steps, dt=dt, v=v, D_r=D_r, D_theta= D_theta, k=k, L=L, mu=mu, w=w)
            particle.dynamics()
            # Retireve data
            x_dyn= np.array(particle.position_x)
            y_dyn = np.array(particle.position_y)
            reactive_path = np.array(particle.reactive_path)
            if is_out==True:
                is_outside = np.logical_and(np.logical_not(particle.bool_reactant), np.logical_not(particle.bool_target))
                x_dyn = x_dyn[is_outside]
                y_dyn = y_dyn[is_outside]
            else:
                x_dyn = x_dyn[reactive_path]
                y_dyn = y_dyn[reactive_path]
                ax.text(target.x, target.y, "T",c = "white", size ="xx-large" , ha = "center", va="center")
                
            #Plot
            h = ax.hist2d(x_dyn, y_dyn, bins = 100,cmap="RdBu_r", density = True)
            ax.text(reactant.x, reactant.y, "R", c = "white", size ="xx-large", ha = "center" , va="center")
            ax.set_xticks(x_ticks, x_labels)
            ax.set_yticks(y_ticks, y_labels)
            fig.colorbar(h[3], ax=ax)

    #fig.tight_layout(pad=0.2)
    return fig


# Compute times from bool array
def calculate_reactive_times(bool_array):
    times = bool_array
    times[times==True] = 1
    times[times==False] = 0
    times = np.where(times==1)[0]
    times = np.array([times[i]-times[i-1]-1 for i in range(1,len(times)) ])
    times = times[times!=0]
    return times


### Compute average reactive path length
def plot_reactive_time(reactant, target, pe_vec, lstar_vec, is_out = False, w=0.0, num_steps=10**5):
    fig =  plt.figure(figsize=(4, 4))
    # Initialize the matrux to save times
    matrix = np.zeros((len(lstar_vec), len(pe_vec)))

    for i in range(len(lstar_vec)):
        for j in range(len(pe_vec)):
            pe = pe_vec[j]
            l_star = lstar_vec[i]
            v = pe*v_max
            D_theta = float(pe*v_max*v_max/l_star)
            # Compute statistics
            particle = ABP_2d(reactant,target, num_steps = num_steps, dt=dt, v=v, D_r=D_r, D_theta= D_theta, k=k, L=L, mu=mu, w=w)
            particle.dynamics()
            
            #Retrieve data
            if is_out==True:
                path = np.logical_and(np.logical_not(particle.bool_reactant), np.logical_not(particle.bool_target))
            else:
                path = np.array(particle.reactive_path)
           
            # Calculate times
            times = calculate_reactive_times(path)
            # Fill matrix
            matrix[i][j] = np.mean(times)    

    plt.xlabel("$Pe$")
    plt.ylabel("$\ell^*$", rotation = 0)
    plt.pcolormesh(pe_vec, lstar_vec, matrix, cmap="RdBu")
    plt.colorbar()
    return fig

### Compute transition rates
def plot_transition_rates(reactant, target, pe_vec, lstar_vec, is_out = False, w=0.0, num_steps=10**5):
    fig =  plt.figure(figsize=(4, 4))
    # Initialize the matrux to save times
    matrix = np.zeros((len(lstar_vec), len(pe_vec)))
    
    for i in range(len(lstar_vec)):
        for j in range(len(pe_vec)):
            pe = pe_vec[j]
            l_star = lstar_vec[i]
            v = pe*v_max
            D_theta = float(pe*v_max*v_max/l_star)
            # Compute statistics
            particle = ABP_2d(reactant,target, num_steps = num_steps, dt=dt, v=v, D_r=D_r, D_theta= D_theta, k=k, L=L, mu=mu, w=w)
            particle.dynamics()
            # Retireve data
            if is_out==True:
                path = np.logical_and(np.logical_not(particle.bool_reactant), np.logical_not(particle.bool_target))
            else:
                path = np.array(particle.transition_path)
            # Calculate times
            times = calculate_reactive_times(path)
            # Fill matrix
            matrix[i][j] = 1./np.mean(times)   

    plt.xlabel("$Pe$")
    plt.ylabel("$\ell^*$", rotation = 0)
    plt.pcolormesh(pe_vec, lstar_vec, matrix, cmap="RdBu_r")
    plt.colorbar()
    return fig

