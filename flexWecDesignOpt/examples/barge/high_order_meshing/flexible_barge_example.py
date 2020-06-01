# input yaml vars at top
import numpy as np

# compute m, k matrices with numpy integration
# list first 8 eigenvalues OR change input files to be for 4 mode RAOs
# substitution classes
# device design np array
import pygmsh

device_name = 'FlexibleBarge'
run_wamit_command = 'C:\WAMITv7\wamit'
common_file_directory = 'C:/Users/13365/Documents/GitHub/flexWecDesignOpt/flexWecDesignOpt/examples/barge/high_order_meshing'
output_directory = 'C:/Users/13365/Desktop/barge_output_test'
cases_file = 'C:/Users/13365/Documents/GitHub/flexWecDesignOpt/flexWecDesignOpt/examples/barge_array.csv'
gmsh_exe_location = 'C:/Users/13365/Documents/gmsh-4.5.6-Windows64/gmsh-4.5.6-Windows64/gmsh'  # take these away
mesh_refinement_factor = 0.25  # take this away


# Free free roots of beam equation
def characteristic_equation(x):
    import math
    return math.cos(x) * math.cosh(x) - 1


class FlexibleBarge(object):

    def __init__(self, design_vars):
        self.device_name = 'barge'
        self.L = design_vars[0]
        self.w = design_vars[1]
        self.h = design_vars[2]

        self.Cg = 0
        self.kx = ((1 / 12) * (self.L ** 2 + self.h ** 2)) ** (1 / 2)
        self.ky = ((1 / 12) * (self.w ** 2 + self.h ** 2)) ** (1 / 2)
        self.kz = ((1 / 12) * (self.L ** 2 + self.w ** 2)) ** (1 / 2)

        self.rho = 500
        self.mass = self.rho * self.L * self.w * self.h
        self.E = 30.720 * 10 ** 6  # Pa
        self.inertia_area = (1 / 12) * self.w * self.h ** 2
        self.nu = 0.3
        self.M = np.asarray([0.25, 0.25, 0.25, 0.25])

        self.eigenvalues = np.asarray([4.7300, 7.8532, 10.9956, 14.1372])
        self.omega = self.eigenvalues / self.L
        self.C = (self.omega ** 2) * self.M



    def substitutions(self):
        return {'L': self.L, 'w': self.w, 'h': self.h, 'Cg': self.Cg,
                'kx': self.kx, 'ky': self.ky, 'kz': self.kz}

    def geometry(self):
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_box(x0=[-1 / 2 * self.L, -1 / 2 * self.w, -1 / 2 * self.h],
                         extents=[self.L, self.w, self.h])
        return geometry
