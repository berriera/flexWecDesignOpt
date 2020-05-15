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