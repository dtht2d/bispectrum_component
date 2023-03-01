import numpy as np
import cmath

class Bispectrum:
    """
    Calculate bispectrum components
    """
    def __init__(self, j1, j2, j,Rcut,neighbour_list):
        neighbour_list = theta_0
        self.j = j
        self.m = m
        self.mp = mp
        self.theta_0 = theta_0
        self.theta = theta
        self.phi = phi
    def _param_generate_(self):
        """
        Generate parameters for bispectrum calculation
        """
        m = np.linspace(-j, j, int(2 * j + 1)).tolist()
        mp = np.linspace(-j, j, int(2 * j + 1)).tolist()
        m1 = np.linspace(-j1, j1, int(2 * j1 + 1)).tolist()
        m1p = np.linspace(-j1, j1, int(2 * j1 + 1)).tolist()
        m2 = np.linspace(-j2, j2, int(2 * j2 + 1)).tolist()
        m2p = np.linspace(-j2, j2, int(2 * j2 + 1)).tolist()
    def __input_check(self):
        """
        Check the validity of input parameters
        """
        if j < 0 or not np.isclose(j, int(j)) or (j % 1 == 0.5 and (m % 1 != 0 or mp % 1 != 0)):
            raise ValueError("Invalid input parameters: j must be a non-negative integer or half-integer, "
                             "m and mp must be between -j and j.")
        if theta_0 < 0 or theta_0 > np.pi or theta < 0 or theta > np.pi or phi < 0 or phi > 2 * np.pi:
            raise ValueError(
                "Invalid input parameters: theta_0, theta, and phi must be within [0, pi] and [0, 2pi], respectively.")