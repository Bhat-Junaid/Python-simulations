############################# Planetary System ################################
import numpy as np
from numba import njit
@njit
def Init():
    GM = [
        1.989e30 * 6.6743e-11,  # star
        1.89813e27 * 6.6743e-11,    # Jupiter
        9.393e20 * 6.6743e-11,    # Ceres
    ]
    return GM


