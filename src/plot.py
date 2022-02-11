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

# Mesh grid
x = np.linspace(-par_dict["L"], par_dict["L"], 1000)
y = np.linspace(-par_dict["L"], par_dict["L"], 1000)
xx, yy = np.meshgrid(x, y)

# Potential
zz = potential(xx,yy)


# Plot
cf = plt.contourf(xx,yy, zz, cmap ="RdBu")
plt.plot(dynamics[:,0], dynamics[:,1], c = "black")
plt.colorbar(cf)

plt.title("Active brownian particle dynamics")
plt.xlabel("x")
plt.ylabel("y")
plt.xticks([-par_dict["L"],0, par_dict["L"] ])
plt.yticks([0, par_dict["L"] ])
plt.tight_layout()
plt.savefig("dynamics.png")