"""
Compare the Wigner D matrix calculation from class function, numba vs SymPy
"""

import numpy as np
from bispectrum.methods.calc.WignerD import Wigner_D, U_rot, u
from sympy import *
from sympy.physics.quantum.spin import Rotation
import timeit
j, m, mp, theta_0, theta,phi= 1, 1, 0, np.pi, np.pi/2, 0

t0=timeit.default_timer()
#Wigner_D function from calc.WignerD
WD = Wigner_D(j, m, mp, theta_0, theta,phi)
wd = WD.wigner_D()
print ("Wigner_D calculation from our function", wd)
t1=timeit.default_timer()
print("Execution time for Wigner_D function from calc.WignerD:", t0 - t1, "seconds")

t2=timeit.default_timer()
#Calculate the Wigner D matrix with Sympy
rot1 = Rotation.D(j, m, mp, theta_0, theta,phi)
wignerD_sympy = rot1.doit()
print("Wigner_D calculation from our SymPy", N(wignerD_sympy))
t3=timeit.default_timer()

print("Execution time for Wigner_D function from Sympy:", t2 - t3, "seconds")
print("Execution time for Wigner_D calculation using class method is", round((t2 - t3)/(t0-t1)),
      "times faster than Sympy function")

#Calculate the Rotational matrix U
U = U_rot(j, m, mp, theta_0, theta, phi)
print(U)
