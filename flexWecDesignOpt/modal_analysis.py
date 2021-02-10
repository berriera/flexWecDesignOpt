def boundary_condition_frequency_solver(function_name, eigenvalue_count=1, h=1.0, freq_limit=50, reltol=1e-5, eps=1e-3,
                                        extra_args=()):
    """
    This function finds a user-specified number of eigenvalue solutions to a device's characteristic equation.

    Args:
        function_name: the name of the system to be iteratively solved
        eigenvalue_count (int): number of unique roots to f(x) = 0 be found
        h (float): increment used to guess new roots; default is 1.0
        freq_limit (float): maximum allowable frequency to be found; used as a break in case of unsuccessful root
                            finding or to only search a certain space; default is 100
        reltol (float): relative tolerance used to find if a new found root is unique; default is 1e-5
        eps(float): used to as a minimum limit on found roots to avoid solutions too close to 0; default is 1e-3
        extra_args (tuple): used to pass extra constant values to fsolve

    Returns:
        nonzero_eigenvalues (numpy array): list of non-zero roots to the equation; shape is 1 row by eigenvalue_count
                                            columns or smaller if less roots in the specified range were found
    """

    import numpy as np
    import scipy.optimize

    x0 = reltol  # used to avoid division by zero errors in solving certain problems with fsolve
    found_eigenvalue_count = 0
    eigenvalues = np.zeros(eigenvalue_count, )
    while found_eigenvalue_count < eigenvalue_count and x0 < freq_limit:
        root = scipy.optimize.fsolve(function_name, x0, extra_args)[0]
        if root > eps and not np.any(abs(eigenvalues - root) / root < reltol):
            eigenvalues[found_eigenvalue_count] = root
            found_eigenvalue_count = found_eigenvalue_count + 1
        x0 = x0 + h
    nonzero_eigenvalues = eigenvalues[eigenvalues != 0]
    return nonzero_eigenvalues


# TODO: rewrite to travel along function_name and make brackets so there's no need for h

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
