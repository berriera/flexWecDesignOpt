import numpy as np
import pygmsh
from flexWecDesignOpt.analysis import eigenvalue_solver


class Device(object):

    def __init__(self, design_vars):
        self.device_name = 'device_name'
        self.var1 = design_vars[0]
        self.var2 = design_vars[1]
        self.var3 = design_vars[2]
        # ...

    def substitutions(self):
        return {'var1name': self.var1, 'var2name': self.var2, 'var3name': self.var3}

    def geometry(self):
        geometry = pygmsh.opencascade.Geometry()
        # TODO: make genericized meshing idea here
        return geometry


class Cylinder(object):

    def __init__(self, design_vars):
        self.device_name = 'cylinder'
        self.r = design_vars[0]
        self.h = design_vars[1]
        self.Cg = 0
        self.kx = ((1 / 12) * (3 * self.r ** 2 + self.h ** 2)) ** (1 / 2)
        self.ky = self.kx
        self.kz = ((1 / 2) * (self.r ** 2)) ** (1 / 2)

    def substitutions(self):
        return np.array([self.r, self.h, self.Cg, self.kx, self.ky, self.kz])

    def geometry(self):
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_cylinder(x0=[0.0, 0.0, 1 / 2 * self.h],
                              axis=[0.0, 0.0, -1 * self.h],
                              radius=self.r,
                              angle=2 * np.pi,
                              char_length=1)
        return geometry


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


def barge_free_free_characteristic_equation(x):
    import math
    return math.cos(x) * math.cosh(x) - 1


# Finds 8 roots of beam equation
eigenvalues = eigenvalue_solver(barge_free_free_characteristic_equation, 8)


class FlexibleBarge(object):

    def __init__(self, design_vars):
        self.device_name = 'barge'
        self.L = design_vars[0]
        self.w = design_vars[1]
        self.h = design_vars[2]

        # Barge material properties
        self.E = 30.720e6  # Modulus of elasticity, Pa
        self.rho = 500  # density, kg / (m ** 3)
        self.nu = 0.3

        # Barge geometry properties
        self.inertia_area = (1 / 12) * self.w * self.h ** 3
        self.cross_section_area = self.h * self.w

        # Barge inertial properties
        self.mass = self.rho * self.L * self.w * self.h
        self.Cg = 0
        self.kx = ((1 / 12) * (self.L ** 2 + self.h ** 2)) ** (1 / 2)
        self.ky = ((1 / 12) * (self.w ** 2 + self.h ** 2)) ** (1 / 2)
        self.kz = ((1 / 12) * (self.L ** 2 + self.w ** 2)) ** (1 / 2)
        self.Ix = self.mass * self.kx ** 2
        self.Iy = self.mass * self.ky ** 2
        self.Iz = self.mass * self.kz ** 2

        # Barge structural deformation mass M and stiffness C matrices
        self.M = self.mass * np.asarray([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25])
        self.eigenvalues = eigenvalues
        self.kappa = self.eigenvalues / self.L
        self.omega = (((self.E * self.inertia_area) / (self.rho * self.cross_section_area * self.L ** 4)) *
                      self.eigenvalues ** 4) ** (1 / 2)
        self.C = (self.omega ** 2) * self.M
        # TODO: figure out where last value M(5,5) in barge.frc comes from

    def substitutions(self):
        # Dictionary of values
        return {'L': self.L, 'w': self.w, 'h': self.h, 'Cg': self.Cg,
                'mass': self.mass,
                'M11': self.M[0], 'M22': self.M[1], 'M33': self.M[2], 'M44': self.M[3],
                'M55': self.M[4], 'M66': self.M[5], 'M77': self.M[6], 'M88': self.M[7],
                'C11': self.C[0], 'C22': self.C[1], 'C33': self.C[2], 'C44': self.C[3],
                'C55': self.C[4], 'C66': self.C[5], 'C77': self.C[6], 'C88': self.C[7]
                }

    def geometry(self):
        # Mesh properties
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_box(x0=[-1 / 2 * self.L, -1 / 2 * self.w, -1 / 2 * self.h],
                         extents=[self.L, self.w, self.h])
        return geometry


class RM3(object):

    def __init__(self, design_vars):
        self.r1 = design_vars[0]
        self.r2 = design_vars[1]
        self.d1 = design_vars[2]
        self.d2 = design_vars[3]

    def substitutions(self):
        return {'r1': self.r1, 'r2': self.r2, 'd1': self.d1, 'd2': self.d2}

    def geometry(self):  # TODO: finish geometry definition here
        return None
