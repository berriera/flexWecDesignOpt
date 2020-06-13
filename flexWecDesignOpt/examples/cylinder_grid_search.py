import numpy as np

r = np.arange(start=0.2, stop=1.20, step=0.2)
h = np.arange(start=0.5, stop=1.20, step=0.2)

independent_vars = 2

row = 0
design_matrix = np.zeros((len(r) * len(h), independent_vars))
for i in range(0, len(r)):
    for j in range(0, len(h)):
        design_matrix[row, :] = [r[i], h[j]]
        row += 1

np.savetxt('cylinder_array.csv', design_matrix, delimiter=',')