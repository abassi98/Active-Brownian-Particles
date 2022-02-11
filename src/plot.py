import matplotlib.pyplot as plt
import numpy as np
import math

# Load files
par = open("parameters.txt")
dynamics = np.loadtxt("dynamics.txt")

# Create parameter dictionsry
par_dict = {}
for line in par:
    key, value = line.split()
    par_dict[key] = float(value)


### Contour plot of potential
# Define potential function
def potential(x, y):
    return par_dict["k"]*(np.sin(8*math.pi*x/par_dict["L"]) + np.sin(8*math.pi*y/par_dict["L"]))

# Take x and y positions
x_dyn = dynamics[:,0]
y_dyn = dynamics[:,1]

# Take start and finish points
x_start = x_dyn[0]
y_start = y_dyn[0]
x_end = x_dyn[-1]
y_end = y_dyn[-1]

# Mesh grid
x = np.linspace(min(x_dyn)*1.1, max(x_dyn)*1.1, 1000)
y = np.linspace(min(y_dyn)*1.1, max(y_dyn)*1.1, 1000)
xx, yy = np.meshgrid(x, y)

# Potential
zz = potential(xx,yy)


### Plotting
# Contourf and dynamics
cf = plt.contourf(xx,yy, zz, cmap ="RdBu")
plt.plot(x_dyn, y_dyn, c = "black")
plt.colorbar(cf)
plt.plot(x_start, y_start, c = "blue", marker = "o", ms = 10, ls  = "",label = "Start" )
plt.plot(x_end, y_end, c = "red", marker = "o", ms = 10, ls  = "",label = "End" )

# Rendering
plt.title("Active brownian particle dynamics")
plt.xlabel("x")
plt.ylabel("y")
plt.xticks([min(x_dyn)*1.1, (min(x_dyn) + max(x_dyn))*0.55,max(x_dyn) ])
plt.yticks([min(y_dyn)*1.1,  (min(y_dyn) + max(y_dyn))*0.55, max(y_dyn)*1.1 ])
plt.legend()
plt.tight_layout()
plt.savefig("dynamics.png")