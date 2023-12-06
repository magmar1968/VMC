import numpy as np
import matplotlib.pyplot as plt

# Load data from the Fortran-formatted file
file_path = 'log.txt'
data = np.loadtxt(file_path)
# Calculate the average
average_value = np.mean(data)

# Plot histogram
plt.hist(data, bins=50, edgecolor='black')
plt.title('Histogram of Data Series')
plt.xlabel('Values')
plt.ylabel('Frequency')
plt.show()

# Print the average
print(f'Average: {average_value}')
