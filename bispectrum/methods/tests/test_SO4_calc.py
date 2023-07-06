from bispectrum.methods.calc.INPUT_param import  neighbor_params
from bispectrum.methods.calc.SO4 import Bispectrum
from sympy import *
import timeit
from sympy.physics.quantum.cg import CG
import numpy as np
path = '/Users/duonghoang/Documents/GitHub/bispectrum_component//data/avgBL-Model.cif'
center_atom_id = 17
r_mu = 0.0779
R_cut = 0.25
input_file_path = path
output_directory = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data'
#read data from file
data = neighbor_params(center_atom_id, r_mu, R_cut, input_file_path, output_directory, file_type='cif')

#possible list of m values
B = Bispectrum(j=2, j1=3/2, j2=1/2, params=data)
keep_list, full_list = B.generate_m_values(3,1,2)
print (keep_list)
j1 = 1
m1 = 1
j2 = 1/2
m2 = 1/2
j = 3/2
m = 3/2
m1p =-1
m2p = 1/2
mp = 1/2

#test Clebsch Gordan Coefficient calculation
CG_coeff = B.clebsch_gordan(j1, j2, j, m1, m2, m)
print(CG_coeff)

CG_coeff_p = B.clebsch_gordan(j1, j2, j, m1p, m2p, mp)
print(CG_coeff_p)
#Compare sympy
cg = CG(j1,m1,j2,m2,j,m)
cg = cg.doit()

#test Coupling coefficient calculation
H_coeff = B.H(j1, j2, j, m1, m2, m, m1p, m2p, mp)
print(H_coeff)

#Test Wigner D
j, m, mp, theta_0, theta,phi= 1, -1, 0, np.pi, np.pi/2, 0
WD = B.wigner_D(j, m, mp, theta_0, theta,phi)
#Test U_rot
U = B.U_rot(j, m, mp, theta_0, theta,phi)
print (WD)
print(U)
#Test u
u = B.u_small(j, m, mp, params=data)
print(u)

#Test bispectrum
B = Bispectrum(j=2, j1=3/2, j2=1/2, params=data)
Bispectrum_SO4 = B.evaluate(j=2, j1=3/2, j2=1/2, params=data)
print(Bispectrum_SO4)

B = Bispectrum(j=2, j1=1/2, j2=3/2, params=data)
Bispectrum_SO4 = B.evaluate(j=2, j1=1/2, j2=3/2, params=data)
print(Bispectrum_SO4)