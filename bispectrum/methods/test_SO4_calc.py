import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import json
#DATA PREPARATION
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

center_atom_id = 17
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

#INPUT values
atomic_radius = 1.46            #silicon atomic radius, unit: angstrom
cell_length = df.iloc[4]        #index row start from 0 _cell_length_a at row 5 index [4]

r_mu = 0.0779                   #scale atomic radius w.r.t cell length
R_cut = 0.25                    #scaled value w.r.t cell length (for Si-Si case)
df_ik = df[(df['r_ik'] + r_mu)<= (R_cut)].copy(deep=True)
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
#array for weight coefficient w.r.t to atom type
w_ik_arr = np.full((r_ik_array.shape),1)
#delta function delta=1 if i and k has the same element type, if not delta =0
delta = np.full((r_ik_array.shape),0)
delta_arr = np.where(df_ik['atom_type']==df_ik['atom_type'].iloc[0],1,delta)
df_ik['w_ik']= w_ik_arr
df_ik['delta']= delta_arr
print (df_ik)
# save the DataFrame as a JSON file
neighbor_list_path = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data/atom-'+ \
                     str(center_atom_id)+'-neighbor-list.json'
#store data in a lightweight format
#df_ik.to_json(neighbor_list_path, orient='records')
# read the .json data from a file
with open(neighbor_list_path, 'r') as f:
    json_data = json.load(f)

# get the keys from the first item in the JSON data
keys = list(json_data[0].keys())

# create an empty NumPy array with the same number of columns as the JSON data
num_columns = len(keys)
data_array = np.empty((len(json_data), num_columns), dtype=object)

# populate the NumPy array with data from the JSON data
for i, row in enumerate(json_data):
    for j, key in enumerate(keys):
        value = row[key]
        if isinstance(value, str) and not value.isdigit():
            data_array[i][j] = value
        else:
            data_array[i][j] = float(value)
print (data_array)
#crete dictionary mapping column name
column_names = ['atom_type', 'X', 'Y', 'Z', 'X_k', 'Y_k', 'Z_k', 'r_ik', 'theta_0', 'theta', 'phi', 'w_ik', 'delta']
column_indices = {name: i for i, name in enumerate(column_names)}
print (column_indices)

#extract data from array
atom_type_array = data_array[:, 1]
r_ik_array = data_array[:, 7]
theta_0_array = data_array[:,8]
theta_array = data_array[:, 9]
phi_array = data_array[:, 10]
w_ik_array = data_array[:, 11]
delta_array = data_array[:, 12]

