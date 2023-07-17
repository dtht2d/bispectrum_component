from bispectrum.methods.calc.INPUT_param import  neighbor_params
from bispectrum.methods.calc.SO4 import *
import pandas as pd
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
df= pd.DataFrame(data)
print(df)
#Test bispectrum
B = Bispectrum(j=5/2, j1=2, j2=1/2, params=data)
Bispectrum_SO4 = B.evaluate(j=5/2, j1=2, j2=1/2, params=data)
print("Bispectrum coefficient for atom ID17", Bispectrum_SO4)


#test generate set+Clebsch-Gordan coefficients
keep_set, full_set = B.generate_m_values(5/2,2,1/2)
keep_set_arr = np.array(keep_set)
print ("test generate set", keep_set_arr.shape)
df = pd.DataFrame(keep_set, columns=["m1", "m2", "m", "m1p", "m2p", "mp"])
unique_combinations = df[["m1", "m2", "m"]].drop_duplicates()
CG_unique =[]
j=5/2
j1=2
j2=1/2
for _, row in unique_combinations.iterrows():
    m, m1, m2 = row["m"], row["m1"], row["m2"]
    CB = B.clebsch_gordan(j1,j2,j,m1,m2,m)
    CG_unique.append(CB)
unique_combinations['CG'] = CG_unique
print(unique_combinations)




