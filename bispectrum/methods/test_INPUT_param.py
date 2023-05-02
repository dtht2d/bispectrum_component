import numpy as np
from itertools import product
import pandas as pd
from bispectrum.methods.calc.INPUT_param import neighbor_params
center_atom_id = 17
r_mu =0.0779
R_cut= 0.25
input_file = "/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif"
output_dir = "/Users/duonghoang/Documents/GitHub/bispectrum_component/data"
file_type = 'cif'
#input_data = neighbor_params(center_atom_id, r_mu, R_cut, input_file, output_dir, file_type='cif')

#print (input_data)
# Convert dictionary to DataFrame
#df = pd.DataFrame(input_data)
# Display DataFrame
#print(df)
j1,j2,j = 1,2,3
#genertate_m_values
m1_vals = np.linspace(-j1, j1, int(2*j1+1))
m2_vals = np.linspace(-j2, j2, int(2*j2+1))
m_vals = np.linspace(-j, j, int(2*j+1))
mp_vals = m_vals.copy()
m1, m2, m, m1p, m2p, mp = np.meshgrid(m1_vals, m2_vals, m_vals, m1_vals, m2_vals, mp_vals)
s = np.stack((m1.ravel(), m2.ravel(), m.ravel(), m1p.ravel(), m2p.ravel(), mp.ravel()), axis=1)
print (s.shape)


def generate_m_val(j1, j2, j):
    m1_vals = np.linspace(-j1, j1, int(2*j1+1))
    m2_vals = np.linspace(-j2, j2, int(2*j2+1))
    m_vals = np.linspace(-j, j, int(2*j+1))
    mp_vals = m_vals.copy()
    m1, m2, m, m1p, m2p, mp = np.meshgrid(m1_vals, m2_vals, m_vals, m1_vals, m2_vals, mp_vals)

    s = np.stack((m1.ravel(), m2.ravel(), m.ravel(), m1p.ravel(), m2p.ravel(), mp.ravel()), axis=1)
