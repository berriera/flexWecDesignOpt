import numpy as np
from flexible_tube import FlexibleTube
import analysis

analysis.evaluate_device(FlexibleTube, np.array([-2, 0.275, 0.003, 15, 1.14e-4, 0.320, 100]))
