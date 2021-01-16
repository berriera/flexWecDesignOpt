def mass_matrix(mode_function, bounds, mode_count, orthogonal=False):
    import numpy as np
    import scipy.integrate

    modal_mass = np.zeros(mode_count)

    def mode_product(x, i, j):
        return mode_function(x, i) * mode_function(x, j)

    for mode_number in range(mode_count):
        x_start = bounds[0]
        x_end = bounds[1]

        # TODO: outer loop through i and j
        # TODO: transpose that transpose
        # TODO: finish modal product
        modal_mass[mode_number] = scipy.integrate.quad(func=mode_product, a=x_start, b=x_end,
                                                       args=(mode_number, mode_number))[0]

    return modal_mass


def barge_modes(x, i):
    import numpy as np
    import math

    eigenvalues = np.array([4.7300, 7.8532, 10.9956, 14.1372])
    eig = eigenvalues[i]
    length = 1
    k = eig / length

    constant = (math.cosh(eig) - math.cos(eig)) / (math.sinh(eig) - math.sin(eig))
    mode_shape = np.cosh(k * x) + np.cos(k * x) - constant * (np.sinh(k * x) + np.sin(k * x))
    end_displacement = np.cosh(eig) + np.cos(eig) - constant * (np.sinh(eig) + np.sin(eig))
    mode_shape = mode_shape / end_displacement

    return mode_shape


# import numpy as np

length = 1
# x_barge = np.linspace(start=0.0, stop=length, num=1001)
# y = barge_modes(x_barge, i=0)
# print(y)

barge_modal_mass = mass_matrix(mode_function=barge_modes, bounds=(0, length), mode_count=4, orthogonal=False)
print('\nBarge modal mass matrix:')
print(barge_modal_mass)
