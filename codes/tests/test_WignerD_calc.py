"""
Compare the Wigner D matrix calculation from my function vs SymPy
"""
import numpy as np
from methods.Bi_SO4 import Wigner_D
from sympy import *
from sympy.physics.wigner import wigner_d
from sympy.physics.quantum.spin import Rotation
import timeit
j, m, mp, theta_0, theta,phi= 1, 1, 0, np.pi, np.pi/2, 0

t0=timeit.default_timer()
#Wigner_D function from Bi_SO4
wigner = Wigner_D(j, m, mp, theta_0, theta, phi)
wd = wigner.wigner_D()
print(wd)
t1=timeit.default_timer()
print("Execution time for Wigner_D function from Bi_SO4:", t0 - t1, "seconds")

t2=timeit.default_timer()
#Calculate the Wigner D matrix with Sympy
rot1 = Rotation.D(j, m, mp, theta_0, theta,phi)
wignerD_sympy = rot1.doit()
print(N(wignerD_sympy))
t3=timeit.default_timer()
print("Execution time for Wigner_D function from Sympy:", t2 - t3, "seconds")
print("Execution time for Wigner_D calculation from our function is", round((t2 - t3)/(t0-t1)),
      "times faster than Sympy function")
