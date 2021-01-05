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
        import numpy as np
        return np.array([self.r, self.h, self.Cg, self.kx, self.ky, self.kz])

    def geometry(self):
        import pygmsh
        import math

        geometry = pygmsh.opencascade.Geometry()
        geometry.add_cylinder(x0=[0.0, 0.0, 1 / 2 * self.h],
                              axis=[0.0, 0.0, -1 * self.h],
                              radius=self.r,
                              angle=2 * math.pi,
                              char_length=1)
        return geometry