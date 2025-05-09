import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

if not os.path.exists('plots'):
	os.mkdir('plots')

# Read the CSV file
df = pd.read_csv("output.csv")

# Find indexes of where the value of last columns have non nan values
dt_indices = df.iloc[:, -1].notna()
dt_indices = dt_indices[dt_indices].index.tolist()
dt = df.iloc[dt_indices, -1]

# print the data frame. Use value at indexes[0] as the starting row and indexes[1]-1 as the ending row
dfs = [df.iloc[dt_indices[i]:dt_indices[i+1], :] for i in range(0, len(dt_indices)-1)]
dfs.append(df.iloc[dt_indices[-1]:, :])

# plot Time,Position,Velocity for all the dataframes
for i in range(0, len(dfs)):
    plt.plot(dfs[i].iloc[:, 0], dfs[i].iloc[:, 1], label='Position (m)')
    plt.plot(dfs[i].iloc[:, 0], dfs[i].iloc[:, 2], label='Velocity (m/s)')
    plt.xlabel('Time')
    plt.ylabel('Position, Velocity')
    plt.title(f'Time, Position and Velocity for dt = {dfs[i].iloc[0, -1]}')
    plt.legend()
    plt.grid()
    plt.savefig('plots/plot_dt=' + str(dfs[i].iloc[0, -1]) + '.png')
    plt.clf()


# Find min and max values of position and velocity and corresponding time
min_pos = [dfs[i].iloc[:, 1].min() for i in range(0, len(dfs))]
max_pos = [dfs[i].iloc[:, 1].max() for i in range(0, len(dfs))]

min_vel = [dfs[i].iloc[:, 2].min() for i in range(0, len(dfs))]
max_vel = [dfs[i].iloc[:, 2].max() for i in range(0, len(dfs))]

min_pos_time = [dfs[i].iloc[:, 0][dfs[i].iloc[:, 1].idxmin()] for i in range(0, len(dfs))]
max_pos_time = [dfs[i].iloc[:, 0][dfs[i].iloc[:, 1].idxmax()] for i in range(0, len(dfs))]

min_vel_time = [dfs[i].iloc[:, 0][dfs[i].iloc[:, 2].idxmin()] for i in range(0, len(dfs))] 
max_vel_time = [dfs[i].iloc[:, 0][dfs[i].iloc[:, 2].idxmax()] for i in range(0, len(dfs))]

# Create a output.txt and write the values
with open('output_2.csv', 'w') as f:
    f.write('dt, min_pos, max_pos, min_pos_time, max_pos_time, min_vel, max_vel, min_vel_time, max_vel_time\n')
    for i in range(0, len(dfs)):
        f.write(f'{dfs[i].iloc[0, -1]}, {min_pos[i]}, {max_pos[i]}, {min_pos_time[i]}, {max_pos_time[i]}, {min_vel[i]}, {max_vel[i]}, {min_vel_time[i]}, {max_vel_time[i]}\n')
    f.close()


# Plot min and max values of velocity for all the dataframes against dt
# Data stored in output_2.csv
num_grid_points = [len(dfs[i]) for i in range(0, len(dfs))]

plt.tight_layout()
plt.plot([1/i for i in num_grid_points], np.abs(min_vel), label='Min Velocity (m/s)')
plt.xlabel('1/N')
plt.ylabel('Velocity')
plt.yscale('log')
plt.xscale('log')
plt.title('Absolute Minimum Velocity for different dt')
plt.legend()
plt.grid()
plt.savefig('plots/minimum_vel_convergence.png')
plt.clf()

plt.tight_layout()
plt.plot([1/i for i in num_grid_points], np.abs(max_vel), label='Max Velocity (m/s)')
plt.xlabel('1/N')
plt.ylabel('Velocity')
plt.yscale('log')
plt.xscale('log')
plt.title('Absolute Maximum Velocity for different dt')
plt.legend()
plt.grid()
plt.savefig('plots/maximum_vel_convergence.png')
plt.clf()

# Plot min and max values of position for all the dataframes against dt
plt.tight_layout()
plt.plot([1/i for i in num_grid_points], np.abs(min_pos), label='Min Position (m)')
plt.xlabel('1/N')
plt.ylabel('Position')
plt.yscale('log')
plt.xscale('log')
plt.title('Absolute Minimum Position for different dt')
plt.legend()
plt.grid()
plt.savefig('plots/minimum_position_convergence.png')
plt.clf()

plt.tight_layout()
plt.plot([1/i for i in num_grid_points], np.abs(max_pos), label='Max Position (m)')
plt.xlabel('1/N')
plt.ylabel('Position')  
plt.yscale('log')
plt.xscale('log')
plt.title('Absolute Maximum Position for different dt')
plt.legend()
plt.grid()
plt.savefig('plots/maximum_position_convergence.png')
plt.clf()
