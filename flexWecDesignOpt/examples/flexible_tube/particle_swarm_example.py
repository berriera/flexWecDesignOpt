import numpy as np
from examples.flexible_tube import Fl
from flexible_tube import FlexibleTube
import analysis

particle_count = 10
np.random.seed(42)
lower_bounds = np.asarray([-5, 0.1, 0.005, 5, 1.10e-4, 0, 110])
upper_bounds = np.asarray([0, 0.5, 0.025, 30, 1.30e-4, 0, 110])

analysis.evaluate_device(FlexibleTube, np.array([-2, 0.275, 0.003, 15, 1.14e-4, 0.320, 100]))
