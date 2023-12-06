import numpy as np

# Load data from the Fortran-formatted file
file_path = 'log.txt'
data = np.loadtxt(file_path)

# Separate data into positive and negative entries in the first column
positive_data = data[data[:, 0] > 0]
negative_data = data[data[:, 0] < 0]

# Calculate average of the second column for positive and negative entries
average_positive = np.mean(positive_data[:, 1])
average_negative = np.mean(negative_data[:, 1])

# Print the results
print(f'Average of the second column for positive entries: {average_positive}')
print(f'Average of the second column for negative entries: {average_negative}')
