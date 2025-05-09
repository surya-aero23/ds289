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

# Calculate average and maximum true errors
average_error = [[] for _ in range(3)]  # 0: EE, 1: RK4, 2: AB2
max_error = [[] for _ in range(3)]

for i, df in enumerate(dfs):
    average_error[0].append(df["Explicit Euler"].sub(df["True Solution"]).abs().mean())
    average_error[1].append(df["RK4"].sub(df["True Solution"]).abs().mean())
    average_error[2].append(df["Adams-Bashforth 2nd order"].sub(df["True Solution"]).abs().mean())

    max_error[0].append(df["Explicit Euler"].sub(df["True Solution"]).abs().max())
    max_error[1].append(df["RK4"].sub(df["True Solution"]).abs().max())
    max_error[2].append(df["Adams-Bashforth 2nd order"].sub(df["True Solution"]).abs().max())


# Plot the true solution and the numerical solutions
for i, df in enumerate(dfs):
    plt.plot(df["Time"], df["True Solution"], label="True Solution", linestyle="-", color="black")
    plt.xlabel("t")
    plt.ylabel("y(t)")
    
    plt.grid()
    if i == len(dfs)-1:
        plt.title(f"Analytical Solution", fontsize="large")
        plt.legend(fontsize="medium")
        plt.tight_layout()
        plt.savefig(f"plots/analytical_solution.png")

    plt.title(f"Delta t: {df.iloc[0, -1]}", fontsize="large")
    plt.plot(df["Time"], df["Explicit Euler"], label="Explicit Euler", linestyle="--", color="blue")
    plt.plot(df["Time"], df["RK4"], label="RK4", linestyle="--", color="red")
    plt.plot(df["Time"], df["Adams-Bashforth 2nd order"], label="Adams-Bashforth 2nd order", linestyle="--", color="green")
    
    plt.tight_layout()
    plt.legend(fontsize="medium")
    plt.savefig(f"plots/comparison_dt={df.iloc[0, -1]}.png")
    plt.close()

'Average Error Plots'
# Use 4 subplots, one for all error comparison and 3 for each method's error comparison with reference lines
fig, axs = plt.subplots(2, 2, figsize=(10, 10), layout = 'constrained')
fig.suptitle("Average Error Comparison (Order of Accuracy)", fontsize="x-large")
# fig.constrained_layout = True

axs[0, 0].plot(dt, average_error[0], label="Explicit Euler", linestyle="-", color="blue")
axs[0, 0].plot(dt, average_error[1], label="RK4", linestyle="-", color="black")
axs[0, 0].plot(dt, average_error[2], label="Adams-Bashforth 2nd order", linestyle="-", color="green")
axs[0, 0].set_xscale("log")
axs[0, 0].set_yscale("log")
axs[0, 0].set_ylabel("Average True Error")
axs[0, 0].grid()
axs[0, 0].legend(fontsize="medium")
# axs[0, 0].set_title("All Methods", fontsize="medium")
axs[0, 0].set_xticks(dt)
axs[0, 0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# reference slope vs individual methods for all other plots, order 1 for explicit euler
axs[0, 1].plot(dt, average_error[0], label="Explicit Euler", linestyle="-", color="blue")
axs[0, 1].plot(dt, 1e-2*np.array(dt)**1, label="m = -1", linestyle="--", color="blue")
axs[0, 1].set_xscale("log")
axs[0, 1].set_yscale("log")
axs[0, 1].grid()
axs[0, 1].legend(fontsize="medium")
# axs[0, 1].set_title("Explicit Euler", fontsize="medium")
axs[0, 1].set_xticks(dt)
axs[0, 1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# order 4 for RK4
axs[1, 0].plot(dt, average_error[1], label="RK4", linestyle="-", color="black")
axs[1, 0].plot(dt, 1e-2*np.array(dt)**4, label="m = -4", linestyle="--", color="black")
axs[1, 0].set_xscale("log")
axs[1, 0].set_yscale("log")
axs[1, 0].set_xlabel("Delta t")
axs[1, 0].set_ylabel("Average True Error")
axs[1, 0].grid()
axs[1, 0].legend(fontsize="medium")
# axs[1, 0].set_title("RK4", fontsize="medium")
axs[1, 0].set_xticks(dt)
axs[1, 0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# order 2 for AB2
axs[1, 1].plot(dt, average_error[2], label="Adams-Bashforth 2nd order", linestyle="-", color="green")
axs[1, 1].plot(dt, 1e-2*np.array(dt)**2, label="m = -2", linestyle="--", color="green")
axs[1, 1].set_xscale("log")
axs[1, 1].set_yscale("log")
axs[1, 1].set_xlabel("Delta t")
axs[1, 1].grid()
axs[1, 1].legend(fontsize="medium")
# axs[1, 1].set_title("Adams Bashforth 2nd Order", fontsize="medium")
axs[1, 1].set_xticks(dt)
axs[1, 1].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.savefig(f"plots/average_true_error_comparison.png")
plt.close()

'Maximum Error Plots'
plt.plot(dt, max_error[0], label="Explicit Euler", linestyle="-", color="blue")
plt.plot(dt, max_error[1], label="RK4", linestyle="-", color="black")
plt.plot(dt, max_error[2], label="Adams-Bashforth 2nd order", linestyle="-", color="green")

plt.xscale("log")
plt.yscale("log")
plt.xticks(dt, dt)
plt.xlabel("Delta t")
plt.ylabel("Maximum True Error")
plt.title("Maximum True Error vs Delta t", fontsize="large")
plt.grid()
plt.legend(fontsize="medium")
plt.tight_layout()
plt.savefig(f"plots/maximum_true_error_comparison.png")
plt.close()


"""
plt.plot(dt, average_error[0], label=f"Explicit Euler", linestyle="-", color="blue")
plt.plot(dt, average_error[1], label=f"RK4", linestyle="-", color="black")
plt.plot(dt, average_error[2], label=f"Adams-Bashforth 2nd order", linestyle="-", color="green")

# draw reference slope lines to confirm order of accuracy
plt.plot(dt, 2e-2*np.array(dt)**1, label="m = -1", linestyle="--", color="blue")
plt.plot(dt, 2e-2*np.array(dt)**2, label="m = -2", linestyle="--", color="green")
plt.plot(dt, 2e-2*np.array(dt)**4, label = "m = -4", linestyle="--", color="black")

plt.xscale("log")
plt.yscale("log")
plt.xticks(dt, dt)
plt.xlabel("Delta t")
plt.ylabel("Average True Error")
plt.title("Average True Error vs Delta t", fontsize="large")
plt.grid()
plt.legend(fontsize="medium")
plt.savefig(f"plots/average_true_error_comparison.png")
plt.close()

"""