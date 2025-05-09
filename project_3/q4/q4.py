import numpy as np
import matplotlib.pyplot as plt

# Constants
Ac = 0.1        # m^2
E = 200 * (10 ** 9)     # N/m^2
L = 10          # m
P = 100         # N/m
del_x = 0.5     # m

N_elements = int(L / del_x)
N_nodes = N_elements + 1

# Building the stiffness matrix A
K_stiffness = np.diag(2 * np.ones(N_nodes)) - np.diag(np.ones(N_nodes - 1), -1) - np.diag(np.ones(N_nodes - 1), 1)
K_stiffness = K_stiffness * (1 / del_x)
K_stiffness[0, 0] = 1
K_stiffness[N_nodes - 1, N_nodes - 1] = 1
K_stiffness[N_nodes - 1, N_nodes - 2] = 0
K_stiffness[0, 1] = 0           

# Building the matrix B
P_forcing = np.full(N_nodes, 2)
P_forcing[0] = 0
P_forcing[N_nodes-1] = 0
P_forcing = - P_forcing * (P * del_x / (2 * Ac * E))

# Solving the system
C_weights = np.linalg.solve(K_stiffness, P_forcing)

x = np.linspace(0, L, N_nodes)
# equal scaling
# plt.axis('equal')
plt.plot(x, C_weights, color='green', label='u(x)')
plt.xlabel('x (m)', fontsize=14)
plt.ylabel('Deflection (m)', fontsize=14)
plt.title('u(x) vs x', fontsize=16)
plt.grid()
plt.legend(fontsize=14)
plt.tight_layout()
# plt.show()
plt.savefig('q4/q4.png')
plt.close()