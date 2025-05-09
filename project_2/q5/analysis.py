import numpy as np
import matplotlib.pyplot as plt
import os

if not os.path.exists('plots'):
	os.mkdir('plots')


# True solution
def true_solution(x, t):
    alpha = 0.5
    return np.sin(x) * np.exp(-alpha * t) + np.sin(4*x) * np.exp(-16 * alpha * t)


# N values
N = ['32', '64', '128', '256']
time = 0.4

all_errors_implicit = []
for n in N:
	# Read the CSV file into a numpy array
	data = np.genfromtxt('N' + n + '_output.csv', delimiter=',')

	# Separate the data into x and y
	x_approx = data[:, 0]
	y_approx = data[:, 1]

	# Calculate the true solution
	y_true = true_solution(x_approx, time)

	# Calculate the error mean
	error = np.abs(y_true - y_approx)
	error_mean = np.mean(error)

	# Append the error mean to the list
	all_errors_implicit.append(error_mean)
	

all_errors_explicit = [0.002115342476180218, 0.0005369162752563879, 0.00013442511904957266, 3.350996732802654e-05]


# Plot the errors
N = [int(n) for n in N]
plt.figure()
plt.plot(N, all_errors_implicit, label='Implicit - CD2', marker='o', color='r')
plt.plot(N, all_errors_explicit, label='Explicit - CD2', marker='o', color='b')

# reference slope of -2
slope = -2
xpoints = [N[0], N[-1]]
ypoints = [all_errors_explicit[0] + 0.5 * all_errors_explicit[0] , (all_errors_explicit[0] + 0.5 * all_errors_explicit[0] ) * (xpoints[1]/xpoints[0])**slope]
plt.plot(xpoints, ypoints, '--', label='slope = -2', color='g')
        

plt.xlabel('N', fontsize=14)
plt.ylabel('Error Mean', fontsize=14)

plt.title('Error Mean vs N', fontsize=16)
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.tight_layout()
plt.savefig('plots/error_Imp_vs_Exp.png')
# plt.show()
plt.close()