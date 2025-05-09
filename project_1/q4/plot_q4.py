import pandas as pd
import matplotlib.pyplot as plt
import os

if not os.path.exists('plots'):
	os.mkdir('plots')

# Read the CSV file into a pandas DataFrame
df = pd.read_csv('output_q4.csv')

# Plot
plt.figure(figsize=(10, 6), layout='constrained')
plt.plot(df['Time'], df['True Solution'], linestyle='--', marker='o', label='True Solution', color='r', linewidth=2)
plt.plot(df['Time'], df['Implicit Euler'], linestyle='-', label='Implicit Euler', color='b')
plt.title('Implicit Euler Scheme, dt = 0.1', fontsize='x-large')
plt.xlabel('x', fontsize='large')
plt.ylabel('y(x)', fontsize='large')
plt.legend(fontsize='large')
plt.grid(True)
plt.savefig(f"plots/implicit_euler_sol.png")
plt.clf()