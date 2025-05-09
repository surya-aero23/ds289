import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

if not os.path.exists('plots'):
	os.mkdir('plots')

s = [100, 500, 1000, 1500]

# Read the data from the csv file (output_q5_S{s[i]}.csv)
for i in range(4):
    df = pd.read_csv(f'output_q5_S{s[i]}.csv')
    plt.plot(df['r'], df['T'], label=f'S = {s[i]}')
    plt.title('Temperature vs r')
    plt.xlabel('r')
    plt.ylabel('Temperature')
    plt.legend()
    plt.grid()
    plt.savefig(f'plots/Temperature_vs_r_S{s[i]}.png')
    plt.clf()
