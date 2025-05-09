# read the data from the file and plot the data using pandas
import matplotlib.pyplot as plt
import numpy as np
import os


# for question 4
directory_q4 = ''

# Create plots directory if it does not exist
if not os.path.exists(directory_q4 + 'plots'):
    os.makedirs(directory_q4 + 'plots')

# Analytical solution:
def true_solution(x, t):
    alpha = 0.5
    return np.sin(x) * np.exp(-alpha * t) + np.sin(4*x) * np.exp(-16 * alpha * t)

# read the input.txt file
grid_sizes = []
rd_values = []
cdChoices = []
with open(directory_q4 + 'input.txt', 'r') as file:
    data = file.readlines()
    for i in range(len(data)):
        data[i] = data[i].strip()
        
        if (data[i] == ''):
            continue
        else:
            if (data[i] == '@ grid sizes'):
                i += 1
                while True:
                    if (data[i] == '@ rd values' or data[i] == '\n' or data[i][0] == '@'):
                        break
                    # remove any \n or newline from data[i]
                    data[i] = data[i].replace('\n', '')   
                    grid_sizes.append(data[i])
                    i += 1
        
            if (data[i] == '@ rd values'):
                i += 1
                while True:
                    if (data[i] == '@ CD scheme (2, 4 and 6)' or data[i] == '\n' or data[i][0] == '@'):
                        break
                    # remove any \n or newline from data[i]
                    data[i] = data[i].replace('\n', '')   
                    # 6 digit precision
                    rd_values.append("{:.6f}".format(float(data[i])))
                    i += 1
            
            if (data[i] == '@ CD scheme (2, 4 and 6)'):
                i += 1
                while True:
                    if (data[i] == '@ Diffusion Coefficient' or data[i] == '\n' or data[i][0] == '@'):
                        break
                    # remove any \n or newline from data[i]
                    data[i] = data[i].replace('\n', '')                    
                    cdChoices.append(data[i])
                    i += 1

print('grid_sizes: ', grid_sizes)
print('rd_values: ', rd_values)
print('cdChoices: ', cdChoices, '\n\n')




for rd in rd_values:
    # Open a new file for each rd to write the error values
    with open(directory_q4 + f'outputs/rd_{rd}_errors.txt', 'w') as file:
        
        for cd in cdChoices:
                file.write(f'@CD = {cd}\n')
                
                for Nx in grid_sizes:   
                    if True:     

                        # read the data from the csv file
                        file_name = directory_q4 + 'outputs/CD_' + cd  +'_Nx_' + Nx + '_rd_' + rd + '.csv'
                        print(f'--- Reading data from {file_name} ---')
                        
                        # read it as an np array
                        data_array = np.genfromtxt(file_name, delimiter=',')
                        data_array = data_array[~np.isnan(data_array).all(axis=1)]
                        data_array = data_array[:, ~np.isnan(data_array).all(axis=0)]
                        
                        # get the time values
                        x_values = data_array[0, :]
                        x_values = x_values[~np.isnan(x_values)]

                        # get the x values
                        time_values = data_array[1:, 0]

                        frequency = len(time_values) / 4
                        time_indices = [i for i in range(0, len(time_values), int(frequency))]
                        time_indices[-1] = len(time_values) - 1

                        # get the u values for the time indices
                        u_values = data_array[1:, 1:]

                        if Nx == grid_sizes[-1]:
                            # plot the data at the time indices
                            fig, ax = plt.subplots()
                            for i in range(len(time_indices)):
                                non_nan_u_values = u_values[time_indices[i]] [~np.isnan(u_values[time_indices[i]])]
                                plot_x_values = x_values[:len(non_nan_u_values)]
                                # mark only every few points
                                ax.plot(plot_x_values, non_nan_u_values, '-', label='t = ' + str(np.round(time_values[time_indices[i]], decimals=2)))
                            
                            ax.set_xlabel('x', fontsize='large')
                            ax.set_ylabel('u', fontsize='large')
                            ax.set_title('CD = ' + cd + ', Nx = ' + Nx + ', rd = ' + rd, fontsize='large')
                            ax.legend(fontsize='large')
                            plt.tight_layout()
                            plt.savefig(directory_q4 + 'plots/CD_' + cd + '_Nx_' + Nx + '_rd_' + rd + '.png')
                            plt.close()          

                        # True solution at x values for the last time index
                        true_u_values = true_solution(x_values, time_values[-1]) 

                        # Calculate the error
                        error = np.abs(true_u_values - u_values[-1])    
                        mean_error = np.mean(error)

                        # print(f'Nx = {Nx}, rd = {rd}, CD = {cd}\nMean error = {mean_error}')
                        
                        file.write(f'{Nx} {mean_error}\n')
                
                file.write('\n\n')
file.close()


# Open the error files and plot the error graphs of size (10, 16)
fig_err, ax_err = plt.subplots(1, 2, figsize=(10, 8), sharey=True, sharex=True)

for rd in rd_values:
    with open(directory_q4 + f'outputs/rd_{rd}_errors.txt', 'r') as file:
        # read the data from the file
        data = file.readlines()
        for i in range(len(data)):
            data[i] = data[i].strip()
        
        # For each cd value, store the error values in a different list
        nx_values = []
        cd_2_errors = []
        cd_4_errors = []
        cd_6_errors = []

        for i in range(len(data)):
            if data[i] == '@CD = 2':
                i += 1
                while True:
                    if data[i] == '@CD = 4' or data[i] == '\n' or len(data[i]) == 0:
                        break
                    data[i] = data[i].split()
                    nx_values.append(int(data[i][0]))
                    cd_2_errors.append(float(data[i][1]))
                    i += 1

            if data[i] == '@CD = 4':
                i += 1
                while True:
                    if data[i] == '@CD = 6' or data[i] == '\n' or len(data[i]) == 0:
                        break
                    data[i] = data[i].split()
                    cd_4_errors.append(float(data[i][1]))
                    i += 1

            if data[i] == '@CD = 6':
                i += 1
                while True:
                    if data[i] == '\n' or len(data[i]) == 0:
                        break
                    data[i] = data[i].split()
                    cd_6_errors.append(float(data[i][1]))
                    i += 1
            
    # close the file
    file.close()
            
        
    # Plot the error graphs
    ax_err[rd_values.index(rd)].plot(nx_values, cd_2_errors,'-', label='CD-2', color='blue', marker='o')
    ax_err[rd_values.index(rd)].plot(nx_values, cd_4_errors, '-', label='CD-4', color='green', marker='^', markersize='10')
    ax_err[rd_values.index(rd)].plot(nx_values, cd_6_errors, '-', label='CD-6', color='red', marker='o')

    # draw reference slopes as dashed lines, use two points to draw the line

    if rd == rd_values[0]:
        slope = -2
        xpoints = [nx_values[0], nx_values[-1]]
        ypoints = [cd_4_errors[0] - 0.2 * cd_2_errors[0] , (cd_4_errors[0] - 0.2 * cd_4_errors[0] ) * (xpoints[1]/xpoints[0])**slope]
        ax_err[rd_values.index(rd)].plot(xpoints, ypoints, '--', color='blue', label='slope = ' + str(slope))

        slope = -4
        ypoints = [cd_4_errors[0] - 0.9 * cd_4_errors[0] , (cd_4_errors[0] - 0.9 * cd_4_errors[0])  * (xpoints[1]/xpoints[0])**slope]
        ax_err[rd_values.index(rd)].plot(xpoints, ypoints, '--', color='green', label='slope = ' + str(slope))

        slope = -6  
        ypoints = [cd_2_errors[0] - 0.1 * cd_2_errors[0] , (cd_2_errors[0] - 0.1 * cd_2_errors[0])  * (xpoints[1]/xpoints[0])**slope]
        ax_err[rd_values.index(rd)].plot(xpoints, ypoints, '--', color='red', label='slope = ' + str(slope))

    if rd == rd_values[1]:
        slope = -2
        xpoints = [nx_values[0], nx_values[-1]]
        ypoints = [cd_2_errors[0] + 0.5 * cd_2_errors[0] , (cd_2_errors[0] + 0.5 * cd_2_errors[0] ) * (xpoints[1]/xpoints[0])**slope]
        ax_err[rd_values.index(rd)].plot(xpoints, ypoints, '--', color='blue', label='slope = ' + str(slope))

        slope = -4
        ypoints = [cd_4_errors[0] - 0.9 * cd_4_errors[0] , (cd_4_errors[0] - 0.9 * cd_4_errors[0])  * (xpoints[1]/xpoints[0])**slope]
        ax_err[rd_values.index(rd)].plot(xpoints, ypoints, '--', color='green', label='slope = ' + str(slope))

        slope = -6  
        ypoints = [cd_6_errors[0] - 0.9 * cd_6_errors[0], (cd_6_errors[0] -  0.9 * cd_4_errors[0])  * (xpoints[1]/xpoints[0])**slope]
        ax_err[rd_values.index(rd)].plot(xpoints, ypoints, '--', color='red', label='slope = ' + str(slope))
    

    ax_err[rd_values.index(rd)].set_xlabel('Nx', fontsize='large')
    ax_err[rd_values.index(rd)].set_ylabel('Mean Error', fontsize='large')
    ax_err[rd_values.index(rd)].set_title('rd = ' + rd, fontsize='large')
    ax_err[rd_values.index(rd)].legend(fontsize='large')

    ax_err[rd_values.index(rd)].set_xscale('log')
    ax_err[rd_values.index(rd)].set_yscale('log')



plt.tight_layout()
# plt.show()
plt.savefig(directory_q4 + 'plots/rd_error_vs_nx.png')
plt.close()
