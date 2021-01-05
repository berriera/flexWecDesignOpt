class FlexibleTube(object):

    import os

    # File information
    defmod_location = r'C:\Users\13365\Desktop\defmod.exe'  # TODO: fix forward and back slashes
    gmsh_exe_location = r'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'
    run_wamit_command = r'C:\WAMITv7\wamit'
    input_file_directory = os.getcwd()
    output_file_directory = r'C:/Users/13365/Desktop/optimization_output'
    mesh_refinement_factor = 0.40

    def __init__(self, design_vars): # TODO: clean up self. points
        # from modal_analysis import mass_matrix
        import math
        import numpy as np
        import scipy.io
        from analysis import boundary_condition_frequency_solver

        # Simulation parameters
        self.water_rho = 1000.0
        self.seafloor_depth = 5.0
        self.resonant_mode_count = 10
        self.degrees_of_freedom = 6 + self.resonant_mode_count
        self.maximum_modal_radial_displacement = 0.005
        self.froude_scaling = 10.0

        # Design variable array
        self.name = 'tube'
        self.depth = design_vars[0]
        self.radius_s = design_vars[1]
        self.thickness = design_vars[2]
        self.length = design_vars[3]
        self.distensibility = design_vars[4]
        self.eta = design_vars[5]
        self.towhead_mass = design_vars[6]
        self.fiber_pretension = 1.8770e4
        self.mooring_stiffness = 510.0
        self.mooring_angle = 30 * (math.pi / 180)
        self.mooring_pretension = 443.4

        # Inertial properties of displaced water
        tube_displaced_volume = math.pi * (self.radius_s ** 2) * self.length
        displaced_water_mass = self.water_rho * tube_displaced_volume
        water_inertia_x = (1 / 2) * displaced_water_mass * (self.radius_s ** 2)
        water_inertia_y = (1 / 12) * displaced_water_mass * (3 * (self.radius_s ** 2) + (self.length ** 2))
        water_inertia_z = water_inertia_y

        # Inertial properties of tube # TODO: go back and check this math
        self.free_vibration_modes = np.ones((self.resonant_mode_count, ), dtype=int)
        tube_elastic_material_rho = 532.6462876469635  # Calculated from Babarit et al. 2017
        tube_volume = (2 * np.pi * self.radius_s * self.thickness) * self.length
        tube_mass = tube_elastic_material_rho * tube_volume
        mass = tube_mass + 2 * self.towhead_mass
        tube_inertia_x = tube_mass * (self.radius_s ** 2) + 2 * ((1 / 2) * self.towhead_mass * (self.radius_s ** 2))
        tube_inertia_y = tube_mass * (1 / 2 * self.radius_s ** 2 + 1 / 12 * self.length ** 2) \
                       + 2 * self.towhead_mass * ((1 / 4) * (self.radius_s ** 2 + self.length ** 2))
        tube_inertia_z = tube_inertia_y

        # Complete mass matrix for rigid body
        mass_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        for i in range(3):
            mass_matrix[i][i] = displaced_water_mass
        mass_matrix[3][3] = water_inertia_x + tube_inertia_x
        mass_matrix[4][4] = water_inertia_y + tube_inertia_y
        mass_matrix[5][5] = water_inertia_z + tube_inertia_z

        # Mass matrix for flexible body modes
        for i in range(self.resonant_mode_count):
            for j in range(self.resonant_mode_count):
                mass_matrix[6 + i][6 + j] = 0  # TOOO: modal mass matrix

        # Damping matrix for rigid and flexible body modes
        damping_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        # TODO: modal damping matrix

        # Stiffness matrix for rigid and flexible body modes
        stiffness_matrix = np.zeros(shape=(self.degrees_of_freedom, self.degrees_of_freedom))
        stiffness_matrix[0][0] = 2 * self.mooring_stiffness * math.cos(self.mooring_angle)
        stiffness_matrix[1][1] = 0.0
        stiffness_matrix[2][2] = 2 * self.mooring_stiffness * math.sin(self.mooring_angle)
        # TODO: modal stiffness matrix
        # TODO: MODE variable for fixed and free modes in .frc WAMTI file



        # Contains radial displacement modal information for computing modal mass, damping, and stiffness values


        # x positions for use in tube radial displacement functions
        x = np.arange(-self.length / 2, self.length / 2, 0.001)

        # Variable renames from above for ease in reading math here
        length = self.length
        dist = self.distensibility
        r_s = self.radius_s
        Ss = math.pi * (self.radius_s ** 2)
        rho = self.water_rho
        ts = self.fiber_pretension
        self.mode_type_count = self.resonant_mode_count // 2

        # From Babarit et al. 2017. Zeroes of these two functions are used to find modal frequencies
        def mode_type_1__boundary_conditions(w, L, Di, Ts, rho):  # TODO: change variable names for scoping
            lowk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
            uppk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
            return math.fabs((lowk * L / 2) * math.tanh(uppk * L / 2) - (uppk * L / 2) * math.tan(lowk * L / 2))

        def mode_type_2__boundary_conditions(w, L, Di, Ts, rho, r_s, M, K_a):  # TODO: remove extra variable here
            lowk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) - 1))
            uppk = math.sqrt(
                ((2 * np.pi) / (Di * Ts)) * (math.sqrt(1 + (Ts * rho * (Di ** 2) * (w ** 2) / math.pi)) + 1))
            return (uppk * L / 2) * math.tanh(uppk * L / 2) + (lowk * L / 2) * math.tan(lowk * L / 2) \
                   - (((w ** 2) * rho * Ss * L) / (-M * (w ** 2) + 2 * K_a)) * (uppk / lowk + lowk / uppk) \
                   * math.tanh(uppk * L / 2) * math.tan(lowk * L / 2)

        # TODO: potentially change from math to numpy
        # TODO: change variable names to lowercase

        # Call modal solver to find roots of nonlinear system. Roots are used entirely to find modal frequency values
        w_type_1 = boundary_condition_frequency_solver(mode_type_1__boundary_conditions,
                                                       self.resonant_mode_count // 2,
                                                       h=0.5,
                                                       freq_limit=30 * math.pi,
                                                       extra_args=(self.length,
                                                                   self.distensibility,
                                                                   self.fiber_pretension,
                                                                   self.water_rho))
        w_type_2 = boundary_condition_frequency_solver(mode_type_2__boundary_conditions,
                                                       self.resonant_mode_count // 2,
                                                       h=0.5,
                                                       freq_limit=30 * math.pi,
                                                       extra_args=(self.length,
                                                                   self.distensibility,
                                                                   self.fiber_pretension,
                                                                   self.water_rho,
                                                                   self.radius_s, mass,
                                                                   self.mooring_stiffness))
        w_vector = np.concatenate((w_type_1, w_type_2), axis=None)

        # Find lk and uk values for each modal frequency
        lk_vector = np.sqrt(((2 * np.pi) / (dist * ts)) * (np.sqrt(1 + (ts * rho * (dist ** 2) * (w_vector ** 2) / math.pi)) - 1))
        uk_vector = np.sqrt(((2 * np.pi) / (dist * ts)) * (np.sqrt(1 + (ts * rho * (dist ** 2) * (w_vector ** 2) / math.pi)) + 1))

        # Setup vectors to calculate normalization factors for each modal shape.
        # Analysis used to calculate normalization factors is based on a Taylor Series expansion of a function used to
        # define a maximum radial displacement in the mode shape; since we define a small maximum displacement and
        # expect a small normalization factor, only up through linear terms are necessary in the expansion
        nrm_type_1 = np.zeros(self.resonant_mode_count // 2)
        nrm_type_2 = np.zeros(self.resonant_mode_count // 2)

        # Use modal shape information to calculate radial displacement function for each mode
        # Each row in dr is for each mode
        dr = np.zeros((np.size(w_vector), np.size(x)))
        for i in range(self.mode_type_count):
            lk = lk_vector[i]
            uk = uk_vector[i]
            shape_type_1 = (lk * np.tanh(uk * length / 2) * np.cos(lk * x) / np.cos(lk * length / 2)
                            - uk * np.tan(lk * length / 2) * np.cosh(uk * x) / np.cosh(uk * length / 2))
            nrm_type_1[i] = (2 / np.max(np.abs(shape_type_1))) * (self.maximum_modal_radial_displacement / r_s)
            S = Ss - Ss * nrm_type_1[i] * shape_type_1
            dr[i, :] = np.sqrt(S / math.pi) - r_s

        for i in range(self.mode_type_count, 2 * self.mode_type_count):
            lk = lk_vector[i]
            uk = uk_vector[i]
            shape_type_2 = lk * uk * (-math.tanh(uk * length / 2) * np.sin(lk * x) / math.cos(lk * length / 2)
                                      + math.tan(lk * length / 2) * np.sinh(uk * x) / math.cosh(uk * length / 2))
            nrm_type_2[i - self.mode_type_count] = (2 / np.max(np.abs(shape_type_2))) \
                                              * (self.maximum_modal_radial_displacement / r_s)
            S = Ss - Ss * nrm_type_2[i - self.mode_type_count] * shape_type_2
            dr[i, :] = np.sqrt(S / math.pi) - r_s
            # TODO: figure out why f modes are similar shapes for different w values

        self.nrm_vector = np.concatenate((nrm_type_1, nrm_type_2), axis=None)

        # Save dr matrix to be read into Matlab for use in calculating objective function values
        scipy.io.savemat('radial_displacement.mat', {'radial_displacement_array': dr})

        # Store initialized design values
        self.mass_matrix = mass_matrix
        self.damping_matrix = damping_matrix
        self.stiffness_matrix = stiffness_matrix
        self.w_vector = w_vector

    def substitutions(self):
        # Variable substitutions for input files
        return {'z_s': self.depth,
                'mode_count': self.degrees_of_freedom,
                'seafloor_depth': self.seafloor_depth,
                'mass_matrix': self.mass_matrix,
                'damping_matrix': self.damping_matrix,
                'stiffness_matrix': self.stiffness_matrix,
                'resonant_mode_count': self.resonant_mode_count,
                'length': self.length,
                'radius': self.radius_s,
                'distensibility': self.distensibility,
                'fiber_pretension': self.fiber_pretension,
                'water_rho': self.water_rho,
                'depth': self.depth,
                'resonant_mode_vector': self.w_vector,
                'normalization_factors': self.nrm_vector,
                'free_vibration_modes': self.free_vibration_modes
                }

    def geometry(self):
        # Mesh of a sunken tube
        import pygmsh
        import math
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_cylinder(x0=[-self.length / 2, 0, self.depth], axis=[self.length, 0, 0],
                              radius=self.radius_s, angle=2 * math.pi, char_length=1)
        return geometry

    def objective(self, output_path):  # TODO: figure out how to call power_calculation from here
        import matlab.engine
        eng = matlab.engine.start_matlab()
        annual_power = eng.power_calculation(output_path, self.froude_scaling, self.maximum_modal_radial_displacement,
                                             float(self.radius_s), 0.2108, 0.323, self.water_rho, 1.0,
                                             float(self.length))
        # TODO: figure out r0, eta damping values
        eng.quit()
        return annual_power
