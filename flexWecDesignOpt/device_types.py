import numpy as np
import pygmsh


class Device(object):

    def __init__(self, design_vars):
        self.device_name = 'device_name'
        self.var1 = design_vars[0]
        self.var2 = design_vars[1]
        self.var3 = design_vars[2]
        # ...

    def substitutions(self):
        return

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


class BargeHighOrderMesh(object):

    def __init__(self, design_vars):
        self.device_name = 'barge'
        self.L = design_vars[0]
        self.w = design_vars[1]
        self.h = design_vars[2]
        self.Cg = 0
        self.kx = ((1 / 12) * (self.L ** 2 + self.h ** 2)) ** (1/2)
        self.ky = ((1 / 12) * (self.w ** 2 + self.h ** 2)) ** (1 / 2)
        self.kz = ((1 / 12) * (self.L ** 2 + self.w ** 2)) ** (1 / 2)

    def substitutions(self):
        return {'L': self.L, 'w': self.w, 'h': self.h, 'Cg': self.Cg,
                'L_half': self.L / 2, 'w_half': self.w / 2, 'h_half': self.h / 2,
                'kx': self.kx, 'ky': self.ky, 'kz': self.kz}

    def geometry(self):
        return None


class BargeLowOrderMesh(object):

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
                'kx': self.kx, 'ky': self.ky, 'kz': self.kz}

    def geometry(self):
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
        return

    def geometry(self): # TODO: finish geometry definition here
        return
