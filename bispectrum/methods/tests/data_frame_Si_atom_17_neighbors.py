import numpy as np
import pandas as pd
import sympy as sp
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_d
from sympy.physics.quantum.spin import Rotation
from sympy import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from bispectrum.methods.plot.neighbors import plot_atoms
from bispectrum.methods.calc.INPUT_param import test_neighbor
import matplotlib.pyplot as plt
import itertools
path = "/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif"
dico = MMCIF2Dict(path)
df_cif = pd.DataFrame.from_dict(dico, orient='index')
x = df_cif.iloc[-3]
y = df_cif.iloc[-2]
z = df_cif.iloc[-1]
atom_type = df_cif.iloc[-4]
x_array = np.array(x[0],dtype=float)
y_array = np.array(y[0],dtype=float)
z_array = np.array(z[0],dtype=float)
atom_type_array = np.array(atom_type[0], dtype=str)
df = pd.DataFrame({"atom_type":atom_type_array,"X" : x_array, "Y":y_array, "Z": z_array})

#Estimate list of potentially atoms in the center cell
df_atoms = df[(df['X'].between(0.5,0.7,inclusive='both'))
                         & (df['Y'].between(0.5,0.7,inclusive='both'))
                         & (df['Z'].between(0.5,0.7, inclusive='both'))]
#print (df_atoms)


# id
x_i = df['X'].iloc[17]
y_i = df['Y'].iloc[17]
z_i = df['Z'].iloc[17]
#print(x_i,y_i,z_i)
X_array = df['X'].to_numpy()
Y_array = df['Y'].to_numpy()
Z_array = df['Z'].to_numpy()
X_k_array = X_array - x_i
Y_k_array = Y_array - y_i
Z_k_array = Z_array - z_i
r_ik= np.sqrt(np.square(X_k_array)+np.square(Y_k_array)+np.square(Z_k_array))
df['X_k'],df['Y_k'], df['Z_k'],df['r_ik']= X_k_array, Y_k_array, Z_k_array, r_ik

#Check to see if chosen center atom coordinate sets to (0,0,0)
print(df.iloc[17])

#INPUT values
atomic_radius = 1.46            #silicon atomic radius, unit: angstrom
cell_length = df.iloc[4]        #index row start from 0 _cell_length_a at row 5 index [4]

r_mu = 0.0779                   #scale atomic radius w.r.t cell length
R_cut = 0.25                    #scaled value w.r.t cell length (for Si-Si case)
df_ik = df[(df['r_ik'] + r_mu)<= (R_cut)].copy(deep=true)
print(df_ik[['X_k', 'Y_k', 'Z_k', 'r_ik']])
print (df_ik[['X', 'Y', 'Z']])
print(df_ik)


