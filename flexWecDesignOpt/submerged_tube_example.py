import pygmsh
import os
import numpy as np
import math

# Internal imports # TODO: change imports
from file_mgmt import create_case_directory
from substitution import create_case_files
from mesh import create_mesh_file
from mesh import submerged_mesh
from analysis import run_wamit
from output import write_dict_to_text_file
from output import convert_array_to_dict

# Executable locations and file setting variables
run_wamit_command = 'C:\WAMITv7\wamit'
gmsh_exe_location = 'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'
common_file_directory = os.path.abspath(os.path.join('examples', 'flexible_tube'))
output_directory = os.path.abspath(os.path.join('examples', 'output'))


class FlexibleTube(object):

    def __init__(self, design_vars):
        from modal_analysis import mass_matrix

        # Simulation parameters
        self.water_rho = 1000
        self.seafloor_depth = 5.0
        self.resonant_mode_count = 10
        self.degrees_of_freedom = 6 + self.resonant_mode_count

        # Design variable array
        self.device_name = 'tube'
        self.depth = design_vars[0]
        self.radius = design_vars[1]
        self.thickness = design_vars[2]
        self.length = design_vars[3]
        self.distensibility = design_vars[4]
        self.eta = design_vars[5]
        self.towhead_mass = design_vars[6]
        self.fiber_pretension = 1.8770e4
        self.mooring_stiffness = 510
        self.mooring_pretension = 443.4

        # Inertial properties of displaced water
        self.tube_displaced_volume = math.pi * (self.radius ** 2) * self.length
        self.displaced_water_mass = self.water_rho * self.tube_displaced_volume
        self.water_Ix = (1 / 2) * self.displaced_water_mass * (self.radius ** 2)
        self.water_Iy = (1 / 12) * self.displaced_water_mass * (3 * (self.radius ** 2) + (self.length ** 2))
        self.water_Iz = self.water_Iy

        # Inertial properties of tube # TODO: go back and check this math
        self.tube_elastic_material_rho = 532.6462876469635  # Calculated from Babarit et al. 2017
        self.tube_volume = (2 * np.pi * self.radius * self.thickness) * self.length
        self.tube_mass = self.tube_elastic_material_rho * self.tube_volume
        self.mass = self.tube_mass + 2 * self.towhead_mass
        self.tube_Ix = self.tube_mass * (self.radius ** 2) + 2 * ((1 / 2) * self.towhead_mass * (self.radius ** 2))
        self.tube_Iy = self.tube_mass * (1 / 2 * self.radius ** 2 + 1 / 12 * self.length ** 2) \
                       + 2 * self.towhead_mass * ((1 / 4) * (self.radius ** 2 + self.length ** 2))
        self.tube_Iz = self.tube_Iy

        # Complete mass matrix for rigid body
        self.mass_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        for i in range(3):
            self.mass_matrix[i][i] = self.displaced_water_mass
        self.mass_matrix[3][3] = self.water_Ix + self.tube_Ix
        self.mass_matrix[4][4] = self.water_Iy + self.tube_Iy
        self.mass_matrix[5][5] = self.water_Iz + self.tube_Iz

        # Mass matrix for flexible body modes
        for i in range(self.resonant_mode_count):
            for j in range(self.resonant_mode_count):
                self.mass_matrix[6 + i][6 + j] = 0  # TOOO: modal mass matrix

        # Damping matrix for rigid and flexible body modes
        self.damping_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        # TODO: modal damping matrix

        # Stiffness matrix for rigid and flexible body modes
        self.stiffness_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        # TODO: modal stiffness matrix

    def substitutions(self, substitution_type='WAMIT'):
        # Variable substitutions for input files
        if substitution_type == 'WAMIT':
            return {'z_s': self.depth, 'mode_count': self.degrees_of_freedom, 'seafloor_depth': self.seafloor_depth,
                    'mass_matrix': self.mass_matrix, 'damping_matrix': self.damping_matrix,
                    'stiffness_matrix': self.stiffness_matrix
                    }
        # Variable substitutions to be written to a .txt file for reading into custom defmod.f file
        elif substitution_type == 'defmod':
            return {'NEWMDS': self.resonant_mode_count, 'L': self.length, 'R': self.radius, 'Di': self.distensibility,
                    'Ts': self.fiber_pretension, 'rho': self.water_rho, 'zs': self.depth
                    }

    def geometry(self):
        # Mesh of a sunken tube
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_cylinder(x0=[-self.length / 2, 0, self.depth], axis=[self.length, 0, 0],
                              radius=self.radius, angle=2 * np.pi, char_length=1)
        return geometry

    def modes(self):
        # Contains radial displacement modal information for computing modal mass, damping, and stiffness values
        import math
        import numpy as np
        import scipy.io
        from analysis import boundary_condition_frequency_solver

        # x positions for use in tube radial displacement functions
        x = np.arange(-self.length / 2, self.length / 2, 0.001)

        # Variable renames from above for ease in reading math here
        L = self.length
        Di = self.distensibility
        r_s = self.radius
        Ss = math.pi * (self.radius ** 2)
        rho = self.water_rho
        Ts = self.fiber_pretension
        maximum_modal_radial_displacement = 0.005  # for modal normalization purposes
        mode_type_count = self.resonant_mode_count // 2

        # From Babarit et al. 2017. Zeroes of these two functions are used to find modal frequencies
        def mode_type_1__boundary_conditions(w, L, Di, Ts, rho):
            lowk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
            uppk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
            return math.fabs((lowk * L / 2) * math.tanh(uppk * L / 2) - (uppk * L / 2) * math.tan(lowk * L / 2))

        def mode_type_2__boundary_conditions(w, L, Di, Ts, rho, r_s, M, K_a):
            S_s = math.pi * (r_s ** 2)
            lowk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
            uppk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
            return (uppk * L / 2) * math.tanh(uppk * L / 2) + (lowk * L / 2) * math.tan(lowk * L / 2) \
                   - (((w ** 2) * rho * S_s * L) / (-M * (w ** 2) + 2 * K_a)) * (uppk / lowk + lowk / uppk) \
                   * math.tanh(uppk * L / 2) * math.tan(lowk * L / 2)

        # TODO: potentially change from math to numpy
        # TODO: change variable names to lowercase

        # Call modal solver to find roots of nonlinear system. Roots are used entirely to find modal frequency values
        w_type_1 = boundary_condition_frequency_solver(mode_type_1__boundary_conditions, self.resonant_mode_count // 2,
                                                       h=0.5, freq_limit=30 * math.pi,
                                                       extra_args=(self.length, self.distensibility,
                                                                   self.fiber_pretension, self.water_rho))
        w_type_2 = boundary_condition_frequency_solver(mode_type_2__boundary_conditions, self.resonant_mode_count // 2,
                                                       h=0.5, freq_limit=30 * math.pi,
                                                       extra_args=(self.length, self.distensibility,
                                                                   self.fiber_pretension, self.water_rho,
                                                                   self.radius, self.mass, self.mooring_stiffness))

        # Combine the two frequency arrays into one
        w = np.concatenate((w_type_1, w_type_2), axis=None)

        # Find lk and uk values for each modal frequency
        lk_vector = np.sqrt(((2 * np.pi) / (Di * Ts)) * (np.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
        uk_vector = np.sqrt(((2 * np.pi) / (Di * Ts)) * (np.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))

        # Setup vectors to calculate normalization factors for each modal shape.
        # Analysis used to calculate normalization factors is based on a Taylor Series expansion of a function used to
        # define a maximum radial displacement in the mode shape; since we define a small maximum displacement and
        # expect a small normalization factor, only up through linear terms are necessary in the expansion
        nrm_type_1 = np.zeros(self.resonant_mode_count // 2)
        nrm_type_2 = np.zeros(self.resonant_mode_count // 2)

        # Use modal shape information to calculate radial displacement function for each mode
        # Each row in dr is for each mode
        dr = np.zeros((np.size(w), np.size(x)))
        for i in range(mode_type_count):
            lk = lk_vector[i]
            uk = uk_vector[i]
            shape_type_1 = (lk * np.tanh(uk * L / 2) * np.cos(lk * x) / np.cos(lk * L / 2)
                            - uk * np.tan(lk * L / 2) * np.cosh(uk * x) / np.cosh(uk * L / 2))
            nrm_type_1[i] = (2 / np.max(np.abs(shape_type_1))) * (maximum_modal_radial_displacement / r_s)
            S = Ss - Ss * nrm_type_1[i] * shape_type_1
            dr[i, :] = np.sqrt(S / math.pi) - r_s

        for i in range(mode_type_count, 2 * mode_type_count):
            lk = lk_vector[i]
            uk = uk_vector[i]
            shape_type_2 = lk * uk * (-math.tanh(uk * L / 2) * np.sin(lk * x) / math.cos(lk * L / 2)
                                      + math.tan(lk * L / 2) * np.sinh(uk * x) / math.cosh(uk * L / 2))
            nrm_type_2[i - mode_type_count] = (2 / np.max(np.abs(shape_type_2))) \
                                              * (maximum_modal_radial_displacement / r_s)
            S = Ss - Ss * nrm_type_2[i - mode_type_count] * shape_type_2
            dr[i, :] = np.sqrt(S / math.pi) - r_s

        normalization_factors = np.concatenate((nrm_type_1, nrm_type_2), axis=None)

        # Save dr matrix to be read into Matlab for use in calculating objective function values
        scipy.io.savemat('radial_displacement.mat', {'radial_displacement_array': dr})
        return w, normalization_factors


# Geometry name and meshing variables
device_name = 'tube'
mesh_refinement_factor = 0.50

# Folder creation
print('Case: ', str(1))
create_case_directory(output_directory, 1)

# Creates a submerged flexible tube object
design_variables = np.asarray([-1, 0.274, 0.01, 10, 1.12e-4, 0, 110])
tube = FlexibleTube(design_variables)
tube_substitutions = tube.substitutions()
tube_geometry = tube.geometry()
frequencies, mode_normalization_factors = tube.modes()

# Variable information to be passed to custom Fortran file defmod.f
mode_variable_dict = tube.substitutions(substitution_type='defmod')
w_dict = convert_array_to_dict('w', frequencies, starting_index=1)
nrm_dict = convert_array_to_dict('nrm', mode_normalization_factors, starting_index=1)
combined_defmod_dict = {**mode_variable_dict, **w_dict, **nrm_dict}
write_dict_to_text_file(combined_defmod_dict, file_name='defmod.txt')

# Create input files in directory and then run WAMIT in that directory
create_case_files(common_file_directory, tube_substitutions)
create_mesh_file(tube_geometry, device_name, gmsh_exe_location, mesh_refinement_factor)
vertices = submerged_mesh(device_name)
# run_wamit(run_wamit_command, flexible_bool=True)

# Output logger
write_dict_to_text_file(tube_substitutions)
print('Done.')
