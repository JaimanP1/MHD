import numpy as np

import matplotlib.pyplot as plt



# File path for the output file

file = '/content/output.txt' 



# Load the data, skipping the header

data = np.loadtxt(file, skiprows=1)



# Extract columns: time steps and k values

time_steps = data[:, 0]

k_values = data[:, 1]



# Create the plot

plt.figure(figsize=(10, 6))

plt.plot(time_steps, k_values, marker='o', color='black')



# Add labels and title

plt.xlabel('Time')

plt.ylabel('Height')

plt.xticks(np.arange(0, 160, step=20))

plt.yticks(np.arange(0, 180, step=20))

plt.title('Null point evolution')

plt.grid(True, linestyle='--', alpha=0.7)



# Show the plot

plt.show()
