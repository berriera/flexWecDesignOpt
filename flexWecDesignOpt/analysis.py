def run_wamit(bem_command, modes_command=''):
    """
    This function runs the boundary element solver command in the current case directory.

    Args:
        bem_command (str): location of the boundary element method executable
        modes_command (str):location of the defmod executable
    Returns:
        None
    """
    import subprocess
    print('\tRunning BEM...')
    subprocess.run([bem_command])

    # Runs user created generalized body modes application and then reruns WAMIT with the created gdf.mod file
    if not modes_command:
        subprocess.run([modes_command])
        subprocess.run([bem_command])
    return

def boundary_condition_frequency_solver(function_name, eigenvalue_count=1, h=1.0, freq_limit=50, reltol=1e-5, eps=1e-3,
                                        extra_args=()):
    """
    This function finds a user-specified number of eigenvalue solutions to a device's characteristic equation.

    Args:
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

    from scipy.optimize import fsolve
    import numpy as np
    x0 = reltol  # used to avoid division by zero errors in solving certain problems with fsolve
    found_eigenvalue_count = 0
    eigenvalues = np.zeros(eigenvalue_count, )
    while found_eigenvalue_count < eigenvalue_count and x0 < freq_limit:
        root = fsolve(function_name, x0, extra_args)[0]
        if root > eps and not np.any(abs(eigenvalues - root) / root < reltol):
            eigenvalues[found_eigenvalue_count] = root
            found_eigenvalue_count = found_eigenvalue_count + 1
        x0 = x0 + h
    nonzero_eigenvalues = eigenvalues[eigenvalues != 0]
    return nonzero_eigenvalues


def boundary_condition_frequency_solver__multidimensional(function_name, eigenvalue_count=1, h=1.0, var_count=1,
                                                          var_increment=0, bounds=None, freq_limit=100, reltol=1e-5,
                                                          eps=1e-2):
    """
    This function finds a user-specified number of eigenvalue solutions to a device's characteristic equation.

    Args:
        function_name (object): name of the equation f to be solved; format of the function should return a value
                                  f(x) for the equation f(x) = 0
        eigenvalue_count (int): number of unique roots to f(x) = 0 be found
        h (float): increment used to guess new roots; default is 1.0
        var_count(int): dimension of the nonlinear system
        var_increment (int): array index of the variable to be incremented by h; used for solving boundary conditions
                            involving more than one variable
        freq_limit (float): maximum allowable frequency to be found; used as a break in case of unsuccessful root
                            finding or to only search a certain space; default is 100
        reltol (float): relative tolerance used to find if a new found root is unique; default is 1e-5
        eps(float): used to as a minimum limit on found roots to avoid solutions too close to 0; default is 1e-3

    Returns:
        nonzero_eigenvalues (numpy array): list of non-zero roots to the equation; shape is n rows by eigenvalue_count
                                            columns or smaller if less roots in the specified range of frequencies were
                                            found to the n dimensional system
    """
    # TODO: add multidimensional test functions
    from scipy.optimize import fsolve
    import numpy as np
    import warnings

    x0 = np.zeros(var_count) + reltol
    found_roots_count = 0
    eigenvalues = np.zeros((eigenvalue_count, var_count))
    i = 0
    while found_roots_count < eigenvalue_count and x0[var_increment] < freq_limit:
        if bounds is not None:
            x0 = bounds[:, 0] + np.random.random(size=var_count, ) * bounds[:, 1]
        x0[var_increment] = i * h
        root = fsolve(function_name, x0)
        if np.all(root - eps > np.zeros(var_count)) and not np.any(abs(eigenvalues - root) / root < reltol):
            eigenvalues[found_roots_count, :] = root
            found_roots_count = found_roots_count + 1
        i = i + 1
    if var_count == 1:
        nonzero_eigenvalues = eigenvalues[eigenvalues != 0.0]
    else:
        nonzero_eigenvalues = eigenvalues[eigenvalues[:, var_increment].argsort()]
    return nonzero_eigenvalues
