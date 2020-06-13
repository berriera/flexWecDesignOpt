def run_wamit(bem_command):
    """
    This function runs the boundary element solver command in the current case directory.

    Args:
        bem_command (str):

    Returns:
        None
    """
    import subprocess
    print('\tRunning BEM...')
    subprocess.run([bem_command])
    return


def eigenvalue_solver(function_name, eigenvalue_count=1, h=1.0, freq_limit=100, reltol=1e-5, eps=1e-3):
    """
    This function finds a user-specified number of eigenvalue solutions to a device's characteristic equation.

    Args:
        function_name (object): name of the equation f to be solved; format of the function should return a value
                                  f(x) for the equation f(x) = 0
        eigenvalue_count (int): number of unique roots to f(x) = 0 be found
        h (float): increment used to guess new roots; default is 1.0
        freq_limit (float): maximum allowable frequency to be found; used as a break in case of unsuccessful root
                            finding or to only search a certain space; default is 100
        reltol (float): relative tolerance used to find if a new found root is unique; default is 1e-5
        eps(float): used to as a minimum limit on found roots to avoid solutions too close to 0; default is 1e-3

    Returns:
        nonzero_eigenvalues (numpy array): list of non-zero roots to the equation; shape is 1 row by eigenvalue_count
                                            columns or smaller if less roots in the specified range were found
    """
    from scipy.optimize import fsolve
    import numpy as np

    x0 = 0
    found_eigenvalue_count = 0
    eigenvalues = np.zeros(eigenvalue_count, )
    while found_eigenvalue_count < eigenvalue_count and x0 < freq_limit:
        root = fsolve(function_name, x0)[0]
        if root > eps and not np.any(abs(eigenvalues - root) / root < reltol):
            eigenvalues[found_eigenvalue_count] = root
            found_eigenvalue_count = found_eigenvalue_count + 1
        x0 = x0 + h
    nonzero_eigenvalues = eigenvalues[eigenvalues != 0]
    return nonzero_eigenvalues

