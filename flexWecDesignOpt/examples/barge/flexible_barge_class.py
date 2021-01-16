class FlexibleBarge(object):
    import os

    # output_directory = os.path.abspath(os.path.join('examples', 'output'))
    gmsh_exe_location = r'C:\Users\13365\Documents\gmsh-4.5.6-Windows64\gmsh-4.5.6-Windows64\gmsh'
    defmod_exe_location = os.path.abspath(r'defmod.exe')
    run_wamit_command = r'C:\WAMITv7\wamit'
    exe_locations = [run_wamit_command, defmod_exe_location, run_wamit_command]

    input_file_directory = os.path.abspath(r'input_files')
    output_file_directory = r'C:\Users\13365\Desktop\optimization_output'
    verbosity = True

    # TODO: build absolute file paths from relative path in separate file file_experimentation.py and then fix here

    def __init__(self, design_vars):
        import numpy as np
        import pygmsh

        self.name = 'barge'
        self.analysis = 'wec-sim'

        # Define design variables from input array
        length = design_vars[0]
        width = design_vars[1]
        height = design_vars[2]
        elasticity_modulus = design_vars[3]  # modulus of elasticity, Pa (30.720e6)
        density = design_vars[4]  # density, kg / (m ** 3) (500)

        # Constant material properties
        nu = 0.3
        mode_count = 8

        # Barge geometry properties
        inertia_area = (1 / 12) * width * height ** 3
        cross_section_area = width * height

        # Barge inertial properties
        mass = density * length * width * height
        Cg = 0
        kx = ((1 / 12) * (length ** 2 + height ** 2)) ** (1 / 2)
        ky = ((1 / 12) * (width ** 2 + height ** 2)) ** (1 / 2)
        kz = ((1 / 12) * (length ** 2 + width ** 2)) ** (1 / 2)
        Ix = mass * kx ** 2
        Iy = mass * ky ** 2
        Iz = mass * kz ** 2

        # Build entire mass matrix from inertial and modal mass matrices
        mass_matrix = np.array([mass, mass, mass, Ix, Iy, Iz])
        modal_mass_matrix = 0.25 * mass * np.ones(shape=(mode_count,))
        external_mass_matrix = np.diag(np.concatenate((mass_matrix, modal_mass_matrix), axis=0))

        # Build entire stiffness matrix
        stiffness_matrix = np.zeros(shape=(6,))
        boundary_conditions_roots = np.array([4.7300, 7.8532, 10.9956, 14.1372, 17.2788, 20.4204, 23.5619, 26.7035])
        kappa = boundary_conditions_roots / length
        modal_frequencies = (((elasticity_modulus * inertia_area) / (density * cross_section_area)) * (kappa ** 4)) ** (
                    1 / 2)
        modal_stiffness_matrix = (modal_frequencies ** 2) * modal_mass_matrix
        external_stiffness_matrix = np.diag(np.concatenate((stiffness_matrix, modal_stiffness_matrix), axis=0))

        # Define substitution dictionary for filling out input files
        self.substitutions = {'mass_matrix': external_mass_matrix,
                              'stiffness_matrix': external_stiffness_matrix,
                              'Cg': Cg
                              }

        # Define .stl file meshing
        barge_geometry = pygmsh.opencascade.Geometry()
        barge_geometry.add_box(x0=[-1/2 * length, -1/2 * width, -1/2 * height], extents=[length, width, height],
                               char_length=2.50)
        self.mesh = barge_geometry

    def objective(self, output_folder_path):  # TODO: finish objective function
        import os
        os.chdir(output_folder_path)
        return


class BargeWamitMesh(object):

    def __init__(self, design_vars):
        self.device_name = 'barge'
        self.L = design_vars[0]
        self.w = design_vars[1]
        self.h = design_vars[2]
        self.Cg = 0
        self.kx = ((1 / 12) * (self.L ** 2 + self.h ** 2)) ** (1 / 2)
        self.ky = ((1 / 12) * (self.w ** 2 + self.h ** 2)) ** (1 / 2)
        self.kz = ((1 / 12) * (self.L ** 2 + self.w ** 2)) ** (1 / 2)

    def substitutions(self):
        return {'L': self.L, 'w': self.w, 'h': self.h, 'Cg': self.Cg,
                'L_half': self.L / 2, 'w_half': self.w / 2, 'h_half': self.h / 2,
                'kx': self.kx, 'ky': self.ky, 'kz': self.kz}

    def geometry(self):
        return None
