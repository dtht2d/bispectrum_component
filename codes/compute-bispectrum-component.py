import numpy as np
import pandas as pd
import sympy as sp
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_d
from sympy.physics.quantum.spin import Rotation
from sympy import *
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from methods import *
import itertools
path = ".../data/avgBL-Model.cif"
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

#Choose a center atom i, in this example we choose atom 'Name'=17 from df_atoms dataframe
atom_i =df.iloc[17]

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

#ANGEL CONVERSION
#theta_0
r_ik_array = df_ik['r_ik'].to_numpy() #r_ik from selected neighbors
r_0_array = np.full((r_ik_array.shape),R_cut)
theta_0_array = np.pi*(np.divide(r_ik_array,r_0_array))
#theta
Z_k_abs_array = np.abs(df_ik['Z_k'].to_numpy())
theta_array = np.arccos(np.divide(Z_k_abs_array,r_ik_array))
#phi
X_k_array = df_ik['X_k'].to_numpy()
Y_k_array = df_ik['Y_k'].to_numpy()
phi_array = np.arctan(np.divide(Y_k_array, X_k_array))
#convert angle to positive value between [0,2pi]
phi_array_convert = np.mod(phi_array, 2*np.pi)
for angle in phi_array_convert:
    if (angle >=2*np.pi) and (angle < 0):
        raise ValueError('phi angle in between 0 and 2pi')
#replace NaN with 0: (code will have error for invalid value center atom values 0/0)_

df_ik['theta_0'] = theta_0_array
df_ik['theta_0'] = df_ik['theta_0'].replace(np.nan,0)
df_ik['theta'] = theta_array
df_ik['theta'] = df_ik['theta'].replace(np.nan,0)
df_ik['phi'] = phi_array_convert
df_ik['phi'] = df_ik['phi'].replace(np.nan,0)

#EXAMPLE
j,m,mp = 3,2,3
#array for weight coefficient w.r.t to atom type
w_ik_arr = np.full((r_ik_array.shape),1)
#delta function delta=1 if i and k has the same element type, if not delta =0
delta = np.full((r_ik_array.shape),0)
delta_arr = np.where(df_ik['atom_type']==df_ik['atom_type'].iloc[0],1,delta)
u_jmmp= getDensityFunction_u(3,2,3,w_ik_arr,delta_arr,r_ik_array,0,R_cut,theta_0_array,theta_array,phi_array)
print(u_jmmp)

#create input set [j1,j2,j,m1,m2,m,m1p,m2p,mp]
j = 3
j1 = 1
j2 = 2
m = np.linspace(-j, j, int(2 * j + 1)).tolist()
mp = np.linspace(-j, j, int(2 * j + 1)).tolist()
m1 = np.linspace(-j1, j1, int(2 * j1 + 1)).tolist()
m1p = np.linspace(-j1, j1, int(2 * j1 + 1)).tolist()
m2 = np.linspace(-j2, j2, int(2 * j2 + 1)).tolist()
m2p = np.linspace(-j2, j2, int(2 * j2 + 1)).tolist()
B_sum = 0
for i in itertools.product(m1,m2,m,m1p,m2p,mp):
    m1, m2, m, m1p, m2p, mp = i
    j,j1,j2=3,1,3
    H = getCoeffH(j1,j2,j,m1,m2,m,m1p,m2p,mp)
    if H==0:
        pass
    else:
        u_jmmp = getDensityFunction_u(j, m, mp, w_ik_arr, delta_arr, r_ik_array, 0,
                                      R_cut, theta_0_array, theta_array,phi_array)
        u1_j1m1m1p = getDensityFunction_u(j1, m1, m1p, w_ik_arr, delta_arr, r_ik_array, 0,
                                          R_cut, theta_0_array,theta_array, phi_array)
        u2_j2m2m2p = getDensityFunction_u(j2, m2, m2p, w_ik_arr, delta_arr, r_ik_array, 0,
                                          R_cut, theta_0_array,theta_array, phi_array)
        B_each = np.conj(u_jmmp) * H * (u1_j1m1m1p) * (u2_j2m2m2p)
        B = N(B_each)
        B_sum += B
print (B_sum)
