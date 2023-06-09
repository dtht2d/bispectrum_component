"""
This module is to get neighbor atoms and  generate all possible parameters for bispectrum calculation
"""
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import numpy as np
import json

def neighbor_params(center_atom_id:int, r_mu, R_cut, input_file:str, output_dir:str, file_type:str):
    """
    This function is used to get input value from cif file and define neighbor atoms parameters
    """
    if file_type == "cif":
        # DATA PREPARATION
        dico = MMCIF2Dict(input_file)
        df_cif = pd.DataFrame.from_dict(dico, orient='index')
        x = df_cif.iloc[-3]
        y = df_cif.iloc[-2]
        z = df_cif.iloc[-1]
        atom_type = df_cif.iloc[-4]
        x_array = np.array(x[0], dtype=float)
        y_array = np.array(y[0], dtype=float)
        z_array = np.array(z[0], dtype=float)
        atom_type_array = np.array(atom_type[0], dtype=str)
        df = pd.DataFrame({"atom_type": atom_type_array, "X": x_array, "Y": y_array, "Z": z_array})
        # Estimate list of potentially atoms in the center cell
        id = center_atom_id
        # id
        x_i = df['X'].iloc[id]
        y_i = df['Y'].iloc[id]
        z_i = df['Z'].iloc[id]
        # print(x_i,y_i,z_i)
        X_array = df['X'].to_numpy()
        Y_array = df['Y'].to_numpy()
        Z_array = df['Z'].to_numpy()
        X_k_array = X_array - x_i
        Y_k_array = Y_array - y_i
        Z_k_array = Z_array - z_i
        r_ik = np.sqrt(np.square(X_k_array) + np.square(Y_k_array) + np.square(Z_k_array))
        df['X_k'], df['Y_k'], df['Z_k'], df['r_ik'] = X_k_array, Y_k_array, Z_k_array, r_ik
        # INPUT values
        cell_length = df.iloc[4]  # index row start from 0 _cell_length_a at row 5 index [4]
        #r_mu = 0.0779  # scale atomic radius w.r.t cell length
        #R_cut = 0.25  # scaled value w.r.t cell length (for Si-Si case)
        condition = (df['r_ik'] + r_mu <= R_cut) & (df['r_ik'] != 0)
        df_ik = df.loc[condition].copy(deep=True)
        # ANGEL CONVERSION
        # theta_0
        r_ik_array = df_ik['r_ik'].to_numpy()  # r_ik from selected neighbors
        r_0_array = np.full((r_ik_array.shape), R_cut)
        theta_0_array = np.pi * (np.divide(r_ik_array, r_0_array))
        # theta
        Z_k_abs_array = np.abs(df_ik['Z_k'].to_numpy())
        theta_array = np.arccos(np.divide(Z_k_abs_array, r_ik_array))
        # phi
        X_k_array = df_ik['X_k'].to_numpy()
        Y_k_array = df_ik['Y_k'].to_numpy()
        phi_array = np.arctan(np.divide(Y_k_array, X_k_array))
        # convert angle to positive value between [0,2pi]
        phi_array_convert = np.mod(phi_array, 2 * np.pi)
        for angle in phi_array_convert:
            if (angle >= 2 * np.pi) and (angle < 0):
                raise ValueError('phi angle in between 0 and 2pi')
        # replace NaN with 0: (code will have error for invalid value center atom values 0/0)_
        df_ik['theta_0'] = theta_0_array
        df_ik['theta_0'] = df_ik['theta_0'].replace(np.nan, 0)
        df_ik['theta'] = theta_array
        df_ik['theta'] = df_ik['theta'].replace(np.nan, 0)
        df_ik['phi'] = phi_array_convert
        df_ik['phi'] = df_ik['phi'].replace(np.nan, 0)
        # array for weight coefficient w.r.t to atom type
        w_ik_arr = np.full((r_ik_array.shape), 1)
        # delta function delta=1 if i and k has the same element type, if not delta =0
        delta = np.full((r_ik_array.shape), 0)
        delta_arr = np.where(df_ik['atom_type'] == df_ik['atom_type'].iloc[0], 1, delta)
        df_ik['w_ik'] = w_ik_arr
        df_ik['delta'] = delta_arr
    # save the DataFrame as a JSON file
        neighbor_list_path = output_dir + '-atom-' + str(center_atom_id) + '-neighbor-list.json'
        df_ik.to_json(neighbor_list_path, orient='records')
        """
        This function is to read in the neighbor list file
        """
        with open(neighbor_list_path, 'r') as f:
            json_data = json.load(f)

            # create a list of dictionaries from the JSON data
            data_list = [dict(row) for row in json_data]

            # extract data from the list of dictionaries
            r_ik_array = np.array([float(row['r_ik']) for row in data_list])
            theta_0_array = np.array([float(row['theta_0']) for row in data_list])
            theta_array = np.array([float(row['theta']) for row in data_list])
            phi_array = np.array([float(row['phi']) for row in data_list])
            w_ik_array = np.array([float(row['w_ik']) for row in data_list])
            delta_array = np.array([float(row['delta']) for row in data_list])
            r_cut_array = np.full((r_ik_array.shape), R_cut)
            # return a dictionary with the extracted data
    else:
        raise ValueError(f"Unsupported file type: {file_type}")
    return {
            'r_ik': r_ik_array,
            'theta_0': theta_0_array,
            'theta': theta_array,
            'phi': phi_array,
            'w_ik': w_ik_array,
            'delta': delta_array,
            'r_cut': r_cut_array
        }
def generate_m_values(j, j1, j2):
    """
    This function generates (m1, m2, m, m1p, m2p, mp) from input set (j1,j2,j)
    and only keeps sets that satisfy the conditions.
    Return: a one and only unique set (m1, m2, m, m1p, m2p, mp)
    """
    J=j1 + j2 + j
    #Condition 1
    if not ((abs(j1 - j2) <= j <= j1 + j2) or isinstance(J,int) or (J >=0)):
        raise ValueError("Condition (1) is not satisfied.")
    #Condition 2
    if not ((isinstance(val,int) or val >=0 for val in [j + j1 - j2,j - j1 + j2,j1 + j2 - j])):
        raise ValueError("(𝑗+𝑗1−𝑗2) and (𝑗−𝑗1+𝑗2) and (𝑗1+𝑗2−𝑗) and is non-negative integer")
    #Condition 3
    if j1 > J or j2 > J or j > J:
        raise ValueError("𝑗1,𝑗2,𝑗 not exceed a positive integer 𝐽=𝑗1+𝑗2+𝑗")
    # Generate m values
    m1_vals = np.linspace(-j1, j1, int(2 * j1 + 1))
    m2_vals = np.linspace(-j2, j2, int(2 * j2 + 1))
    m_vals = np.linspace(-j, j, int(2 * j + 1))
    mp_vals = m_vals.copy()
    m1p_vals = m1_vals.copy()
    m2p_vals = m2_vals.copy()
    m1, m2, m, m1p, m2p, mp = np.meshgrid(m1_vals, m2_vals, m_vals, m1p_vals, m2p_vals, mp_vals)
    s = np.stack((m1.ravel(), m2.ravel(), m.ravel(), m1p.ravel(), m2p.ravel(), mp.ravel()), axis=1)
    full_list = s
    keep_list = []
    for i in range(len(s)):
      m1_val, m2_val, m_val, m1p_val, m2p_val, mp_val = s[i]
      #Condition 4-8: Clebsch-Gordan calc for set (𝑗1,𝑗2,𝑗,𝑚1,𝑚2,𝑚), (𝑗1,𝑗2,𝑗,𝑚1p,𝑚2p,𝑚p)
      c4 = [isinstance(val, (int, float)) and val % 0.5 == 0 for val in [m1_val, m2_val, m_val, m1p_val, m2p_val, mp_val]]
      c5 = [isinstance(vals, int) and vals >= 0 for vals in [j1 + m1_val, j1 - m1_val, j2 + m2_val, j2 - m2_val, j + m_val, j - m_val]]
      c5_p = [isinstance(vals, int) and vals >= 0 for vals in [j1 + m1p_val, j1 - m1p_val, j2 + m2p_val, j2 - m2p_val, j + mp_val, j - mp_val]]
      c6 = [m1_val+m2_val == m_val and m1p_val+m2p_val == mp_val]
      c7 = [abs(val) <= limit for val, limit in [(m1_val, j1), (m2_val, j2), (m_val, j), (m1p_val, j1), (m2p_val, j2), (mp_val, j)]]
      c8 = [j2 + j + m1_val >= 0 and j1 - j2 - m_val >= 0 and isinstance(j2 + j + m1_val, int) and isinstance(j1 - j2 - m_val, int)]
      c8_p = [j2 + j + m1p_val >= 0 and j1 - j2 - mp_val >= 0 and isinstance(j2 + j + m1p_val, int) and isinstance(j1 - j2 - mp_val, int)]
      #Condition 9: Wigner-D calc
      c9 = [isinstance(vals, int) and vals >= 0 for vals in [mp_val-m_val, m1p_val-m1_val, m2p_val-m2_val]]
    if (c4) and (c5 and c5_p) and (c6) and (c7) and (c8 and c8_p) and (c9):
      keep_list.append(s[i])
    else:
      pass
    return keep_list, full_list