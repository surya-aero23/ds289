import numpy as np
import matplotlib.pyplot as plt
import os

# If running from terminal, keep the directory_q1 as ''
# If running from IDE, keep the directory_q1 as 'q1'

directory_q1 = 'q1'


N = ['N_64', 'N_1024']
schemes = ['CD2_scheme', 'upwind_scheme']
time_instances = [0, 0.0252, 0.05, 0.075]
sols = []

x_64, x_1024 = [], []
cd2_1024 = []
upwind_64, upwind_1024 = [], []

for n in N:
    for scheme in schemes:

        if scheme == 'CD2_scheme' and n == 'N_64':
                continue
        else:
            
            print('--- Output file: ' + directory_q1 + '/output_' + n + '_' + scheme + '.txt ---')
            data = np.genfromtxt(directory_q1 + '/output_' + n + '_' + scheme + '.txt', delimiter = ',')

            # get the x values at the very first row of the data
            x = data[0]
            x = x[~np.isnan(x)]


            if n == 'N_64':
                x_64 = x
            if n == 'N_1024':
                x_1024 = x
            
            # get all rows for time instances as the first entry
            for time in time_instances:
                entry = data[data[:,0] == time]
                # remove nan values
                entry = entry[~np.isnan(entry)]

                # remove the first entry
                entry = entry[1:]
                
                sols.append(entry)


            if scheme == 'CD2_scheme':
                if n == 'N_1024':
                    cd2_1024 = sols
            
            if scheme == 'upwind_scheme':
                if n == 'N_64':
                    upwind_64 = sols
                if n == 'N_1024':
                    upwind_1024 = sols


            sols = []



'Part (a):'
# compare upwind 64 and upwind 1024
# Use 4 subplot figure and plot the solution at t = 0, 0.0252, 0.05, 0.075

fig, axs = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)

for i in range(4):
    axs[i//2, i%2].plot(x_64, upwind_64[i], label='Upwind 64')
    axs[i//2, i%2].plot(x_1024, upwind_1024[i], label='Upwind 1024')
    axs[i//2, i%2].set_title('t = ' + str(time_instances[i]))
    axs[i//2, i%2].legend(fontsize=12)
    axs[i//2, i%2].grid()

    # set x and y labels
    if i == 2 or i == 3:
        axs[i//2, i%2].set_xlabel('x', fontsize=12) 
    
    if i == 0 or i == 2:
        axs[i//2, i%2].set_ylabel('U', fontsize=12)

plt.suptitle('Upwind Scheme: N = 64 vs N = 1024', fontsize=16)
plt.tight_layout()
plt.savefig(directory_q1 + '/upwind_64_vs_1024.png')
# plt.show()
plt.close()



'Part (b):'

# compare cd2 1024 and upwind 1024
# Use 4 subplot figure and plot the solution at t = 0, 0.0252, 0.05, 0.075

fig, axs = plt.subplots(2, 2, figsize=(10, 10), sharex=True, sharey=True)

for i in range(4):
    axs[i//2, i%2].plot(x_1024, cd2_1024[i], label='CD2 1024')
    axs[i//2, i%2].plot(x_1024, upwind_1024[i], label='Upwind 1024')
    axs[i//2, i%2].set_title('t = ' + str(time_instances[i]))
    axs[i//2, i%2].legend(fontsize=12)
    axs[i//2, i%2].grid()

    # set x and y labels
    if i == 2 or i == 3:
        axs[i//2, i%2].set_xlabel('x', fontsize=12) 
    
    if i == 0 or i == 2:
        axs[i//2, i%2].set_ylabel('U', fontsize=12)

plt.suptitle('CD2 Scheme vs Upwind Scheme: N = 1024', fontsize=16)
plt.tight_layout()
plt.savefig(directory_q1 + '/cd2_vs_upwind_1024.png')
# plt.show()
plt.close()