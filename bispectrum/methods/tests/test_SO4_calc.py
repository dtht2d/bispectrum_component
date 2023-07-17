from bispectrum.methods.calc.INPUT_param import  neighbor_params
from bispectrum.methods.calc.SO4 import *
import pandas as pd
from sympy import *
import timeit
from sympy.physics.quantum.cg import CG
import numpy as np
import timeit

t0=timeit.default_timer()
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
j,j1,j2= 5/2,2,1/2
B = Bispectrum(j, j1, j2,data)
Bispectrum_SO4 = B.evaluate(j, j1, j2,data)
t1=timeit.default_timer()
print("Bispectrum coefficient for atom ID17", np.round(Bispectrum_SO4,2))
print("Execution time for bispectrum calculation", round(t1-t0,2), "seconds")


t2=timeit.default_timer()
#test generate set+Clebsch-Gordan coefficients
keep_set, full_set = B.generate_m_values(j,j1,j2)
keep_set_arr = np.array(keep_set)
print ("Test generate set", keep_set_arr.shape)
print (pd.DataFrame(keep_set))
df = pd.DataFrame(keep_set, columns=["m1", "m2", "m", "m1p", "m2p", "mp"])
unique_combinations = df[["m1", "m2", "m"]].drop_duplicates()
CG_unique =[]
for _, row in unique_combinations.iterrows():
    m, m1, m2 = row["m"], row["m1"], row["m2"]
    CB = B.clebsch_gordan(j1,j2,j,m1,m2,m)
    CG_unique.append(CB)
unique_combinations['CG'] = CG_unique
print(unique_combinations)
t3=timeit.default_timer()
print("Execution time for CB calculation")
print("for unique sets:", round(t3-t2,2), "seconds")




