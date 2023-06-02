"""
This file is used to study the example of SO(4) bispectrum calculation.
A compound where the central atom has only two symmetrical neighbors of the same atom type is known as a
diatomic molecule. In diatomic molecules, the central atom is bonded to two identical atoms.
"""
from bispectrum.methods.calc.SO4 import Bispectrum
import numpy as np
data = {
    'w_ik': np.array([1.0, 1.0]),
    'delta': np.array([1.0, 1.0]),
    'r_ik': np.array([0.5, 0.5]),
    'r_cut': np.array([0.6, 0.6]),
    'theta_0': np.array([np.pi/6, np.pi/6]),
    'theta': np.array([np.pi/3, np.pi/3]),
    'phi': np.array([np.pi/5, np.pi/5])
}

B= Bispectrum(j=3, j1=5/2, j2=1/2, params=data)
B1 = B.evaluate(j=3, j1=5/2, j2=1/2, params=data)
print(B1)