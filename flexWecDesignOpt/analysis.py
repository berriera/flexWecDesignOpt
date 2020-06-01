def run_wamit(output_folder, bem_command):
    """
    This function runs the boundary element solver command in the current case directory.

    Args:
        output_folder(str):
        bem_command (str):

    Returns:
        None
    """
    import os
    import subprocess
    os.chdir(output_folder)
    subprocess.run([bem_command])
    return


def eigenvalue_solver(function_name, eigenvalue_count=1, h=1, freq_limit=100, reltol=1e-5, eps=1e-3):
    from scipy.optimize import fsolve
    import numpy as np
    import sys

    solver_function = getattr(sys.modules[__name__], function_name)
    x0 = 0
    found_eigenvalue_count = 0
    eigenvalues = np.zeros(eigenvalue_count, )
    while found_eigenvalue_count < eigenvalue_count:
        root = fsolve(solver_function, x0)[0]
        if root > eps and not np.any(abs(eigenvalues - root) / root < reltol):
            eigenvalues[found_eigenvalue_count] = root
            found_eigenvalue_count = found_eigenvalue_count + 1
        x0 = x0 + h
    return eigenvalues



eigs = eigenvalue_solver('free_free_barge', eigenvalue_count=8)
print(eigs)
