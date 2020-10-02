import pygmsh
import os
import numpy as np
import math
from file_mgmt import create_case_directory
from substitution import create_case_files
from mesh import create_mesh_file
from mesh import submerged_mesh
from analysis import run_wamit

run_wamit_command = 'C:\WAMITv7\wamit'
gmsh_exe_location = 'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'
common_file_directory = os.path.abspath(os.path.join('examples', 'flexible_tube'))
output_directory = os.path.abspath(os.path.join('examples', 'output'))
device_name = 'tube'
mesh_refinement_factor = 0.50


class FlexibleTube(object):

    def __init__(self, design_vars):
        # Simulation parameters
        self.seafloor_depth = 5.0
        self.resonant_mode_count = 1
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
        self.water_rho = 1000
        self.mooring_stiffness = 510
        self.mooring_pretension = 443.4

        # Resonant frequency information # TODO: take this out
        w = np.asarray(
            [2.1114360764249511, 4.5933899504176336, 7.7089354146509326, 11.6086922395909120, 16.3723249238653743,
             0.5489434835683205, 1.4244105181770936, 1.9121722822149714, 3.9035022567988560, 6.6024612335389854])
        nrm = np.asarray(
            [0.0522554138385799, 0.0241890418770306, 0.0152970427382452, 0.0110450820083542, 0.0086118105373954,
             -0.0814725919983839, -0.0218307770211154, -0.0230966022927663, -0.0104376725023571, -0.0058682851582704])

        # Mass properties
        self.rho = 532.6462876469635  # Calculated from Babarit et al. 2017
        self.tube_volume = 2 * np.pi * self.radius * self.thickness * self.length
        self.tube_displaced_volume = math.pi * (self.radius ** 2) * self.length
        self.tube_mass = self.rho * self.tube_volume
        self.mass = self.tube_mass + 2 * self.towhead_mass
        self.displaced_water_mass = self.water_rho * self.tube_displaced_volume
        self.Ix = self.tube_mass * (self.radius ** 2) + 2 * ((1 / 2) * self.towhead_mass * (self.radius ** 2))
        self.Iy = self.tube_mass * (1 / 2 * self.radius ** 2 + 1 / 12 * self.length ** 2) \
                  + 2 * self.towhead_mass * ((1 / 4) * (self.radius ** 2 + self.length ** 2))
        self.Iz = self.Iy

        # Complete mass matrix for rigid body
        self.mass_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        for i in range(3):
            self.mass_matrix[i][i] = self.displaced_water_mass
        self.mass_matrix[3][3] = self.Ix
        self.mass_matrix[4][4] = self.Iy
        self.mass_matrix[5][5] = self.Iz

        for i in range(self.resonant_mode_count):
            for j in range(self.resonant_mode_count):
                self.mass_matrix[6 + i][6 + j] = 2 # TOOO: modal mass matrix

        self.damping_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom)) # TODO: modal damping matrix
        self.stiffness_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom)) # TODO: modal stiffness matrix

    def substitutions(self):
        return {'z_s': self.depth, 'mode_count': self.degrees_of_freedom, 'seafloor_depth': self.seafloor_depth,
                'mass_matrix': self.mass_matrix, 'damping_matrix': self.damping_matrix,
                'stiffness_matrix': self.stiffness_matrix
                }

    def geometry(self):
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_cylinder(x0=[-self.length / 2, 0, self.depth], axis=[self.length, 0, 0],
                              radius=self.radius, angle=2 * np.pi, char_length=1)
        return geometry

    def modes(self):
        import math
        import numpy as np
        from analysis import boundary_condition_frequency_solver

        x = np.arange(-self.length / 2, self.length / 2, 0.001)

        L = self.length
        Di = self.distensibility
        r = self.radius
        Ss = math.pi * (self.radius ** 2)
        rho = 1000
        Ts = self.fiber_pretension

        def mode_type_1__boundary_conditions(w, L, Di, Ts, rho):
            lk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
            uk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
            return math.fabs((lk * L / 2) * math.tanh(uk * L / 2) - (uk * L / 2) * math.tan(lk * L / 2))

        def mode_type_2__boundary_conditions(w, L, Di, Ts, rho, r_s, M, K_a):
            S_s = math.pi * (r_s ** 2)
            lk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
            uk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
            return (uk * L / 2) * math.tanh(uk * L / 2) + (lk * L / 2) * math.tan(lk * L / 2) \
                   - (((w ** 2) * rho * S_s * L) / (-M * (w ** 2) + 2 * K_a)) * (uk / lk + lk / uk) \
                   * math.tanh(uk * L / 2) * math.tan(lk * L / 2)

        # TODO: potentially change from math to numpy

        lk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
        uk = math.sqrt(((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))

        def mode_1__shape(lk, uk, L, Ss, r, nrm):
            # First mode type
            S = Ss - Ss * nrm * (lk * math.tanh(uk * L / 2) * np.cos(lk * x) / math.cos(lk * L / 2)
                                 - uk * math.tan(lk * L / 2) * np.cosh(uk * x) / math.cosh(uk * L / 2))
            dr = np.sqrt(S / math.pi) - r
            return dr

        def mode_2__shape(lk, uk, L, Ss, r, nrm):
            S = Ss - Ss * nrm * lk * uk * (-math.tanh(uk * L / 2) * np.sin(lk * x) / math.cos(lk * L / 2)
                                           + math.tan(lk * L / 2) * np.sinh(uk * x) / math.cosh(uk * L / 2))
            dr = np.sqrt(S / math.pi) - r
            return dr


        w_type_1 = boundary_condition_frequency_solver(mode_type_1__boundary_conditions, 5, h=0.5, freq_limit=20 * math.pi,
                                                       extra_args=(self.length, self.distensibility,
                                                                   self.fiber_pretension, self.water_rho))
        w_type_2 = boundary_condition_frequency_solver(mode_type_2__boundary_conditions, 5, h=0.5, freq_limit=20 * math.pi,
                                                       extra_args=(self.length, self.distensibility,
                                                                   self.fiber_pretension, self.water_rho,
                                                                   self.radius, self.mass, 510))


# Design variables
print('Case: ', str(1))
design_variables = np.asarray([-1, 0.274, 0.01, 10, 0, 0, 110])
tube = FlexibleTube(design_variables)
tube_substitutions = tube.substitutions()
tube_geometry = tube.geometry()
create_case_directory(output_directory, 1)
create_case_files(common_file_directory, tube_substitutions)
create_mesh_file(tube_geometry, device_name, gmsh_exe_location, mesh_refinement_factor)
vertices = submerged_mesh(device_name)
# run_wamit(run_wamit_command, flexible_bool=True)
print('Done.')
