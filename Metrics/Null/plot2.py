#use with attempt 7

import numpy as np

import matplotlib.pyplot as plt



# Load the data (skip the header, treat missing values as NaN)

data = np.genfromtxt('output.csv', delimiter=',', skip_header=1)



# Split columns

timesteps = data[:, 0]

k1 = data[:, 1]

k2 = data[:, 2]

k3 = data[:, 3]

# Plot

plt.figure(figsize=(10, 6))

plt.plot(timesteps, k1, label='1st sign change (k1)', marker='o')

plt.plot(timesteps, k2, label='2nd sign change (k2)', marker='x')

plt.plot(timesteps, k3, label='3rd sign change (k3)', marker='o')

plt.xlabel("Timestep")

plt.ylabel("k-index")

plt.ylim(0,375)

plt.title("Sign Changes Positions Over Time")

plt.legend()

plt.grid(True)

plt.tight_layout()

plt.show()


