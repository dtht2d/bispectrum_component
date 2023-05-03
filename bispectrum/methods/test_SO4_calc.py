from bispectrum.methods.calc.INPUT_param import  neighbor_params
from bispectrum.methods.calc.SO4 import B
from sympy import *
import timeit
from sympy.physics.quantum.cg import CG
path = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif'
center_atom_id = 17
r_mu = 0.0779
R_cut = 0.25
input_file_path = path
output_directory = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data'
#read data from file
data = neighbor_params(center_atom_id, r_mu, R_cut, input_file_path, output_directory, file_type='cif')

#possible list of m values
B = B(j=1,j1=2,j2=3, params=data)
list_generate = B.generate_m_values(3,1,2)
print (list_generate)
j1, j2, j, m1, m2, m = 1/2, 1/2, 1, -1/2, -1/2, -1
CG_coeff = B.cg(j1, j2, j, m1, m2, m)
print(CG_coeff)

#Compare sympy
t0=timeit.default_timer()
cg = CG(j1,m1,j2,m2,j,m)
cg = cg.doit()
time for CG function from Sympy:", t1-t0, "seconds")