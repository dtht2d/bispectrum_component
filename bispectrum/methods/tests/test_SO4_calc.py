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
print ("test generate set", keep_list)


#test Clebsch Gordan Coefficient calculation
j1, j2, j, m1, m2, m = 2,1/2,5/2,1,1/2,3/2
CG_coeff = B.clebsch_gordan(j1, j2, j, m1, m2, m)
print("test CG calculationg from our function", CG_coeff)

#Compare sympy
cg = CG(j1,m1,j2,m2,j,m)
cg = cg.doit()
print("CG calculation from Sympy",cg)

#test Coupling coefficient calculation
m1p, m2p, mp = 1, 1/2, 3/2
H_coeff = B.H(j1, j2, j, m1, m2, m, m1p, m2p, mp)
print("test coupling coefficient calculation",H_coeff)

#Test Wigner D
j1, m1, m1p, theta_0, theta,phi= 3/2, 3/2, -3/2, np.pi/4, np.pi/3, np.pi/4
WD = B.wigner_D(j1, m1, m1p, theta_0, theta,phi)
print ("test Wigner_D calculation",WD)

#Test U_rot
U = B.U_rot(j1, m1, m1p, theta_0, theta,phi)
print("test U_rot calculation", U)

#Test u
u = B.u_small(j, m, mp, params=data)
print("test u calculation",u)


#Test bispectrum
B = Bispectrum(j=2, j1=3/2, j2=1/2, params=data)
Bispectrum_SO4 = B.evaluate(j=2, j1=3/2, j2=1/2, params=data)
print("Bispectrum coeeficient", Bispectrum_SO4)


