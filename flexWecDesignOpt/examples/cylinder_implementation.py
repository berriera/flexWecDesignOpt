import numpy as np

r = np.arange(start=0.2, stop=1.20, step=0.2)
h = np.arange(start=0.5, stop=1.7, step=0.2)

independent_vars = 2
dependent_vars = 4

row = 0
design_matrix = np.zeros((len(r) * len(h), independent_vars + dependent_vars))
for i in range(0, len(r)):
    for j in range(0, len(h)):
        kx = ((1 / 12) * (3 * r[i] ** 2 + h[j] ** 2)) ** (1 / 2)
        ky = kx
        kz = ((1 / 2) * (r[i] ** 2)) ** (1 / 2)
        Cg = 0
        design_matrix[row, :] = [r[i], h[j], Cg, kx, ky, kz]
        row += 1

np.savetxt('cylinder_array.csv', design_matrix, delimiter=',')
