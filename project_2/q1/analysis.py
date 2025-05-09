import numpy as np
import matplotlib.pyplot as plt
import os

directory_q1 = ''

if not os.path.exists(directory_q1 + 'plots'):
    os.makedirs(directory_q1 + 'plots')

# open T.csv
T = np.genfromtxt(directory_q1 + 'Temperature.csv', delimiter=',')

# from input.txt obtain Lx, Ly, Nx, Ny
var_list = []
with open(directory_q1 + 'input.txt', 'r') as f:

    data = f.readlines()
    for i in range(len(data)):
        f_i = data[i].strip()
        # ignore empty lines and lines starting with @
        if f_i == '' or f_i[0] == '@':
            continue
        else:
            var_list.append(f_i)

# extract Lx, Ly, Nx, Ny
Lx, Ly, Nx, Ny = [int(var_list[i].split()[0]) for i in range(4)]

dx = Lx/(Nx + 1)
dy = Ly/(Ny + 1)

# create x and y arrays
x = np.arange(2, 2 + Lx + (0.1 * dx), dx)
y = np.arange(4, 4 + Ly+ (0.1 * dy), dy)


# plot the temperature distribution as a contour plot

plt.figure(figsize=(10, 6))
plt.contourf(x, y, T, cmap='hot', levels=40)
plt.axis('scaled')
plt.xticks(np.arange(2, 2 + Lx + 0.2, 0.5))
plt.yticks(np.arange(4, 4 + Ly + 0.2, 0.5))
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Temperature distribution', fontsize=14)
plt.tight_layout()
plt.savefig('plots/temperature_contour.png')
