"""
Compare the Wigner D matrix calculation with the one from Sympy, PyXtal_FF and my own
"""
import numpy as np
from methods.Bi_SO4 import Wigner_D as wd
j, m, mp, theta_0, theta,phi= 1, 1, 0, np.pi, np.pi/2, 0

wD = wd.wigner_D((j, m, mp, theta_0, theta, phi))
print(wD)
