"""
Compare the Wigner D matrix calculation from class function
"""

import numpy as np
from bispectrum.methods.calc.WignerD import *
#from sympy import *
#from sympy.physics.quantum.spin import Rotation
import timeit
j, m, mp, theta_0, theta,phi= 3/2, 3/2, -3/2, np.pi/4, np.pi/3, np.pi/4

t0=timeit.default_timer()
#Wigner_D function from calc.WignerD
WD = Wigner_D(j, m, mp, theta_0, theta,phi)
wd = WD.wigner_D()
print ("Wigner_D calculation from our function", wd)
t1=timeit.default_timer()
print("Execution time for Wigner_D function from calc.WignerD:", t1 - t0, "seconds")

#test D(phi,theta, -phi) D(theta,-theta, -phi)
WD_1 = Wigner_D(j, m, mp, phi, theta,-phi)
wd_1 = WD_1.wigner_D()
print ("Wigner_D term 1", wd_1)

WD_1p = Wigner_D(j, m, mp, phi, -theta,-phi)
wd_1p = WD_1p.wigner_D()
print ("Wigner_D term 3", wd_1p)

#test d_small
WD_1= Wigner_D(j, m, mp, theta_0, theta,phi)
d_small = WD_1.compute_dsmall()
print ("d_small from calc", d_small)
d_small_table = -(np.sin(theta/2))**3
print ("d_small from table", d_small)

#Case1 j=3/2, m=3/2, mp=-3/2
j, m, mp, theta_0, theta,phi= 3/2, 3/2, -3/2, np.pi/3, np.pi/2, np.pi/4

#Test rotaion matrix U_rot
U_calc_1 = U_rot(j, m, mp, theta_0, theta,phi)
U_table_1= 1j*((np.sin(theta_0/2)*np.sin(theta)*(np.exp(-1j*phi)))**3)
print("U_table_1", U_table_1)
print("U_calc_1", U_calc_1)



