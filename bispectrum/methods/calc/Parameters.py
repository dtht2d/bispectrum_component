"""
This module is to generate all possible parameters for bispectrum calculation
"""
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import numpy as np
def generate_m_val(j1, j2, j):
    """
    This function generates (m1, m2, m, m1p, m2p, mp) from input set (j1,j2,j)
    and only keep set that satisfy the condition
    """
    # Generate m values
    m = np.linspace(-j, j, int(2 * j + 1))
    mp = m.copy()
    m1 = np.linspace(-j1, j1, int(2 * j1 + 1))
    m1p = m1.copy()
    m2 = np.linspace(-j2, j2, int(2 * j2 + 1))
    m2p = m2.copy()
    s = product(m1, m2, m, m1p, m2p, mp)
    keep_list = []
    for i in s:
        m1, m2, m, m1p, m2p, mp = i[0], i[1], i[2], i[3], i[4], i[5]
        # Check input parameter conditions
        # Condition 1 & 2 & 5
        if not (abs(j1 - j2) <= j <= j1 + j2 and m1 + m2 == m and j1 + j2 - j % 1 != 0.5):
            pass
        if not (abs(j1 - j2) <= j <= j1 + j2 and m1p + m2p == mp and
                j1 + j2 - j % 1 != 0.5):
            pass
        # Condition 3 & 6
        if not all(abs(x) <= y and (x % 0.5 == 0 or x % 1 == 0) for x, y in
                   zip([m1, m2, m], [j1, j2, j])):
            pass
        if not all(abs(x) <= y and (x % 0.5 == 0 or x % 1 == 0) for x, y in
                   zip([m1p, m2p, mp], [j1, j2, j])):
            pass
        # Condition 4
        J = (j1 + j2 + j)
        if J < (int(j1 + j2 + j)) and J < 0:
            raise ValueError("Invalid input parameters: j1, j2, j \
                              must not exceed a positive integer J")
        # Condition 7
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and x >= 0
                   and (x % 0.5 == 0 or x % 1 == 0) for x in [j1, j2, j]):
            raise ValueError("Invalid input parameters: j1, j2, j must be integer \
                             or half-integer non-negative numbers")
        # Condition 8
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and x >= 0
                   and x % 1 == 0 for x in [j1 + m1, j2 + m2, j + m, j1 + j2 + j]):
            pass
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and x >= 0
                   and x % 1 == 0 for x in [j1 + m1p, j2 + m2p, j + mp, j1 + j2 + j]):
            pass
        else:
            keep_list.append(i)
    return keep_list
def neighbor_atoms (input_file, center_atom_id, rcut):
    '''
    This function is to create a neighbor_list for a chosen center atoms
    Parameters:
        input_file (string) directory of the input file
        ceter_atom_id (integer): center atom ID
        rcut (float): cutoff radius depends on the size of the unit cell and atom type
                note:choosen rcut needs to divide by the true cell length
                    since atom coordinate (x,y,z)are fraction with cell dimension (1,1,1)
    Returns: neighbor atoms pandas frame with atom ID and x,y,z fractional distance from center atom
    '''
    path = "/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif"
    dico = MMCIF2Dict(path)
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

    # id
    x_i = df['X'].iloc[center_atom_id]
    y_i = df['Y'].iloc[center_atom_id]
    z_i = df['Z'].iloc[center_atom_id]
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
    df_ik = df[(df['r_ik']) <= (rcut)].copy(deep=True)

    # ANGEL CONVERSION
    # theta_0
    r_ik_array = df_ik['r_ik'].to_numpy()  # r_ik from selected neighbors
    r_0_array = np.full((r_ik_array.shape), rcut)
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
    return df_ik
