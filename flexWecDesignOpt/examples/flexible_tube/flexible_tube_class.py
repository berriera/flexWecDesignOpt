class FlexibleTube(object):
    import os

    # File information
    defmod_location = r'C:\Users\13365\Desktop\defmod.exe'
    gmsh_exe_location = r'C:\Users\13365\Documents\gmsh-4.5.6-Windows64\gmsh-4.5.6-Windows64\gmsh'
    run_wamit_command = r'C:\WAMITv7\wamit'
    input_file_directory = os.getcwd()
    output_file_directory = r'C:\Users\13365\Desktop\optimization_output'

    def __init__(self, design_vars):  # TODO: clean up self. points
        # from modal_analysis import mass_matrix
        import math
        import numpy as np
        import scipy.io
        import pygmsh
        import matlab.engine

        from analysis import boundary_condition_frequency_solver

        # Simulation parameters
        water_rho = 1000.0
        resonant_mode_count = 10
        degrees_of_freedom = 6 + resonant_mode_count
        maximum_modal_radial_displacement = 0.005
        froude_scaling = 10.0

        # Design variable array
        self.name = 'tube'
        depth = design_vars[0]
        radius_s = design_vars[1]
        thickness = design_vars[2]
        length = design_vars[3]
        distensibility = design_vars[4]
        eta = design_vars[5]
        towhead_mass = design_vars[6]
        fiber_pretension = 1.8770e4
        mooring_stiffness = 510.0
        mooring_angle = 30 * (math.pi / 180)
        mooring_pretension = 443.4

        # Inertial properties of displaced water
        tube_displaced_volume = math.pi * (radius_s ** 2) * length
        displaced_water_mass = water_rho * tube_displaced_volume  # TODO: make correction for depth or add bounds for fully submerged
        water_inertia_x = (1 / 2) * displaced_water_mass * (radius_s ** 2)
        water_inertia_y = (1 / 12) * displaced_water_mass * (3 * (radius_s ** 2) + (length ** 2))
        water_inertia_z = water_inertia_y

        # Inertial properties of tube # TODO: go back and check this math
        free_vibration_modes = np.ones((resonant_mode_count,), dtype=int)
        tube_elastic_material_rho = 532.6462876469635  # Calculated from Babarit et al. 2017
        tube_volume = (2 * np.pi * radius_s * thickness) * length
        tube_mass = tube_elastic_material_rho * tube_volume
        mass = tube_mass + 2 * towhead_mass
        tube_inertia_x = tube_mass * (radius_s ** 2) + 2 * ((1 / 2) * towhead_mass * (radius_s ** 2))
        tube_inertia_y = tube_mass * (1 / 2 * radius_s ** 2 + 1 / 12 * length ** 2) \
                         + 2 * towhead_mass * ((1 / 4) * (radius_s ** 2 + length ** 2))
        tube_inertia_z = tube_inertia_y

        # Complete mass matrix for rigid body
        inertia_x = water_inertia_x + tube_inertia_x
        inertia_y = water_inertia_y + tube_inertia_y
        inertia_z = water_inertia_z + tube_inertia_z

        mass_matrix = np.array([displaced_water_mass, displaced_water_mass, displaced_water_mass,
                                inertia_x, inertia_y, inertia_z
                                ])

        modal_mass_matrix = np.zeros(shape=resonant_mode_count, )
        external_mass_matrix = np.diag(np.concatenate((mass_matrix, modal_mass_matrix), axis=0))
        # TODO: finish modal mass matrix

        # Damping matrix for rigid and flexible body modes
        damping_matrix = np.zeros(shape=(degrees_of_freedom, degrees_of_freedom))
        # TODO: modal damping matrix

        # Stiffness matrix for rigid and flexible body modes
        stiffness_matrix = np.zeros(shape=(degrees_of_freedom, degrees_of_freedom))
        stiffness_matrix[0][0] = 2 * mooring_stiffness * math.cos(mooring_angle)
        stiffness_matrix[1][1] = 0.0
        stiffness_matrix[2][2] = 2 * mooring_stiffness * math.sin(mooring_angle)
        # TODO: modal stiffness matrix
        # TODO: MODE variable for fixed and free modes in .frc WAMTI file

        # Contains radial displacement modal information for computing modal mass, damping, and stiffness values

        # x positions for use in tube radial displacement functions
        x = np.arange(-length / 2, length / 2, 0.001)

        # Variable renames from above for ease in reading math here
        Ss = math.pi * (radius_s ** 2)  # TODO: check this vs 2*pi*r*t
        mode_type_count = resonant_mode_count // 2

        # From Babarit et al. 2017. Zeroes of these two functions are used to find modal frequencies
        def mode_type_1__boundary_conditions(w, L, Di, ts, density):  # TODO: change variable names for scoping
            low_k = math.sqrt(
                ((2 * np.pi) / (Di * ts)) * (math.sqrt(1 + (ts * density * (Di ** 2) * (w ** 2) / math.pi)) - 1))
            upp_k = math.sqrt(
                ((2 * np.pi) / (Di * ts)) * (math.sqrt(1 + (ts * density * (Di ** 2) * (w ** 2) / math.pi)) + 1))
            return math.fabs((low_k * L / 2) * math.tanh(upp_k * L / 2) - (upp_k * L / 2) * math.tan(low_k * L / 2))

        def mode_type_2__boundary_conditions(w, L, di, ts, density, s_s, m, k_a):
            low_k = math.sqrt(
                ((2 * np.pi) / (di * ts)) * (math.sqrt(1 + (ts * density * (di ** 2) * (w ** 2) / math.pi)) - 1))
            upp_k = math.sqrt(
                ((2 * np.pi) / (di * ts)) * (math.sqrt(1 + (ts * density * (di ** 2) * (w ** 2) / math.pi)) + 1))
            return (upp_k * L / 2) * math.tanh(upp_k * L / 2) + (low_k * L / 2) * math.tan(low_k * L / 2) \
                   - (((w ** 2) * density * s_s * L) / (-m * (w ** 2) + 2 * k_a)) * (upp_k / low_k + low_k / upp_k) \
                   * math.tanh(upp_k * L / 2) * math.tan(low_k * L / 2)

        # TODO: change variable names to lowercase

        # Call modal solver to find roots of nonlinear system. Roots are used entirely to find modal frequency values
        w_type_1 = boundary_condition_frequency_solver(mode_type_1__boundary_conditions,
                                                       resonant_mode_count // 2,
                                                       h=0.5,
                                                       freq_limit=30 * math.pi,
                                                       extra_args=(length,
                                                                   distensibility,
                                                                   fiber_pretension,
                                                                   water_rho))
        w_type_2 = boundary_condition_frequency_solver(mode_type_2__boundary_conditions,
                                                       resonant_mode_count // 2,
                                                       h=0.5,
                                                       freq_limit=30 * math.pi,
                                                       extra_args=(length,
                                                                   distensibility,
                                                                   fiber_pretension,
                                                                   water_rho,
                                                                   Ss,
                                                                   mass,
                                                                   mooring_stiffness))
        w_vector = np.concatenate((w_type_1, w_type_2), axis=None)

        # Find lk and uk values for each modal frequency
        lk_vector = np.sqrt(
            ((2 * np.pi) / (distensibility * fiber_pretension))
            * (np.sqrt(1 + (fiber_pretension * water_rho * (distensibility ** 2) * (w_vector ** 2) / math.pi)) - 1))
        uk_vector = np.sqrt(
            ((2 * np.pi) / (distensibility * fiber_pretension))
            * (np.sqrt(1 + (fiber_pretension * water_rho * (distensibility ** 2) * (w_vector ** 2) / math.pi)) + 1))

        # Setup vectors to calculate normalization factors for each modal shape.
        # Analysis used to calculate normalization factors is based on a Taylor Series expansion of a function used to
        # define a maximum radial displacement in the mode shape; since we define a small maximum displacement and
        # expect a small normalization factor, only up through linear terms are necessary in the expansion
        nrm_type_1 = np.zeros(resonant_mode_count // 2)
        nrm_type_2 = np.zeros(resonant_mode_count // 2)

        # Use modal shape information to calculate radial displacement function for each mode
        # Each row in dr is for each mode
        dr = np.zeros((np.size(w_vector), np.size(x)))
        for i in range(mode_type_count):
            lk = lk_vector[i]
            uk = uk_vector[i]
            shape_type_1 = (lk * np.tanh(uk * length / 2) * np.cos(lk * x) / np.cos(lk * length / 2)
                            - uk * np.tan(lk * length / 2) * np.cosh(uk * x) / np.cosh(uk * length / 2))
            nrm_type_1[i] = (2 / np.max(np.abs(shape_type_1))) * (maximum_modal_radial_displacement / r_s)
            S = Ss - Ss * nrm_type_1[i] * shape_type_1
            dr[i, :] = np.sqrt(S / math.pi) - radius_s

        for i in range(mode_type_count, 2 * mode_type_count):
            lk = lk_vector[i]
            uk = uk_vector[i]
            shape_type_2 = lk * uk * (-math.tanh(uk * length / 2) * np.sin(lk * x) / math.cos(lk * length / 2)
                                      + math.tan(lk * length / 2) * np.sinh(uk * x) / math.cosh(uk * length / 2))
            nrm_type_2[i - mode_type_count] = (2 / np.max(np.abs(shape_type_2))) \
                                              * (maximum_modal_radial_displacement / radius_s)
            S = Ss - Ss * nrm_type_2[i - mode_type_count] * shape_type_2
            dr[i, :] = np.sqrt(S / math.pi) - radius_s
            # TODO: figure out why f modes are similar shapes for different w values (check for negative uk and lk values as a start?)

        nrm_vector = np.concatenate((nrm_type_1, nrm_type_2), axis=None)

        # Save dr matrix to be read into Matlab for use in calculating objective function values
        scipy.io.savemat('radial_displacement.mat', {'radial_displacement_array': dr})

        # Variable substitutions for input_files files

        self.substitutions = {'z_s': depth,
                              'mode_count': degrees_of_freedom,
                              'mass_matrix': external_mass_matrix,
                              'damping_matrix': damping_matrix,
                              'stiffness_matrix': stiffness_matrix,
                              'resonant_mode_count': resonant_mode_count,
                              'length': length,
                              'radius': radius_s,
                              'distensibility': distensibility,
                              'fiber_pretension': fiber_pretension,
                              'water_rho': water_rho,
                              'depth': depth,
                              'resonant_mode_vector': w_vector,
                              'normalization_factors': nrm_vector,
                              'free_vibration_modes': free_vibration_modes
                              }

        mesh = pygmsh.opencascade.Geometry()
        mesh.add_cylinder(x0=[-length / 2, 0, depth], axis=[length, 0, 0],
                          radius=radius_s, angle=2 * math.pi, char_length=0.10)
        self.mesh = mesh

    def mesh_tube(self):
        pass

    def tube_substitutions(self):
        pass

    def objective(self, output_path, design_parameters):
        import matlab.engine
        eng = matlab.engine.start_matlab()

        froude_scaling = design_parameters[0]
        maximum_modal_radial_displacement = design_parameters[1]
        radius_s = design_parameters[2]
        water_rho = design_parameters[3]
        length = design_parameters[4]

        annual_power = eng.power_calculation(output_path, froude_scaling, maximum_modal_radial_displacement,
                                             float(radius_s), 0.2108, 0.323, water_rho, 1.0,
                                             float(length))
        # TODO: figure out r0, eta damping values
        eng.quit()
        return annual_power
