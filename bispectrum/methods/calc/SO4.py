import numpy as np
import json
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from typing import Dict
# create a lookup table for factorials
factorials = [1]
for i in range(1, 1001):
    factorials.append(factorials[-1] * i)
def fact(n):
    """
    This function is used to calculate factorial of a number by using look the lookup table method
    """
    if n < 0:
        return None
    elif n <= 1000:
        return factorials[n]
    else:
        return float('inf')  # factorial too large for lookup table
class Bispectrum:
    """
    Calculate bispectrum components
    """
    def __init__(self, j, j1, j2, params: Dict[str, np.ndarray] ):
        '''
            j: j index
            j1: j1 index
            j2: j2 index
            input_data: input data dictionary with extracted values:
            r_ik (array): dictance from center atom to n neighbor atom, dim = [n,]
            theta_0 (array): first angle of rotation [0, pi] , dim = [n,]
            theta (array): second angle of rotation [0, pi], dim = [n,]
            phi (array): third angle of rotation [0, 2pi], dim = [n,]
            w_ik (array): weight coefficient, dim = [n,]
            delta (array): delta function, dim = [n,]
            r_cut (array): cutoff distance, dim = [n,]
        '''
        self.j = j
        self.j1 = j1
        self.j2 = j2
        self.input_val = input_data #dictionary with values
        self.params = params['params']
        w_ik_array = self.params['w_ik']
        delta_array = self.params['delta']
        r_ik_array = self.params['r_ik']
        r_cut_array = self.params['r_cut']
        theta_0_array = self.params['theta_0']
        theta_array = self.params['theta']
        phi_array = self.params['phi']

        #Condition check
        if not (abs(j1 - j2) <= j <= j1 + j2 and j1 + j2 - j % 1 != 0.5):
            raise ValueError("Invalid input parameters: j1, j2, j must satisfy the triangle inequality.\ "
                             "j1+j2-j must not be a half-integer")
        J = (j1 + j2 + j)
        if J < (int(j1 + j2 + j)) and J < 0:
            raise ValueError("Invalid input parameters: j1, j2, j must not exceed a positive integer J")

    def u(self, j, m, mp):
      '''
      Parameter:
        j (scalar): angular momentum
        m (scalar): eigenvalue of angular momentum
        mp (scalar): eigenvalue of j along rotated axis
      Returns: expansion coefficients density function u_jm_mp
       '''
      w_ik_array = self.params['w_ik']
      delta_array = self.params['delta']
      r_ik_array = self.params['r_ik']
      r_cut_array = self.params['r_cut']
      theta_0_array = self.params['theta_0']
      theta_array = self.params['theta']
      phi_array = self.params['phi']

      # Calculate cutoff_function
      f_cut_arr = (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array, r_cut_array))) + 1)
      # Calculate rotational matrix U for all k=n neighbor atoms
      U_ik_array = np.array([Wigner_D(self.j, m, mp, theta_0, theta, phi).U_rot for theta_0, theta, phi in zip(theta_0_array, theta_array

class Bispectrum:
    """
    Calculate bispectrum components
    """
    def __init__(self, j, j1, j2, params: Dict[str, np.ndarray]):
        '''
            j: j index
            j1: j1 index
            j2: j2 index
            input_data: input data dictionary with extracted values:
            r_ik (array): dictance from center atom to n neighbor atom, dim = [n,]
            theta_0 (array): first angle of rotation [0, pi] , dim = [n,]
            theta (array): second angle of rotation [0, pi], dim = [n,]
            phi (array): third angle of rotation [0, 2pi], dim = [n,]
            w_ik (array): weight coefficient, dim = [n,]
            delta (array): delta function, dim = [n,]
            r_cut (array): cutoff distance, dim = [n,]
        '''
        self.j = j
        self.j1 = j1
        self.j2 = j2
        self.input_val = input_data
        #Condition check
        if not (abs(j1 - j2) <= j <= j1 + j2 and j1 + j2 - j % 1 != 0.5):
            raise ValueError("Invalid input parameters: j1, j2, j must satisfy the triangle inequality.\ "
                             "j1+j2-j must not be a half-integer")
        J = (j1 + j2 + j)
        if J < (int(j1 + j2 + j)) and J < 0:
            raise ValueError("Invalid input parameters: j1, j2, j must not exceed a positive integer J")
    def generate_m_val(self):
        """
        This function generates (m1, m2, m, m1p, m2p, mp) from input set (j1,j2,j)
        and only keep set that satisfy the condition
        """
        # Generate m values
        m = np.linspace(-self.j, self.j, int(2 * self.j + 1))
        mp = m.copy()
        m1 = np.linspace(-self.j1, self. j1, int(2 * self.j1 + 1))
        m1p = m1.copy()
        m2 = np.linspace(-self.j2, self.j2, int(2 * self.j2 + 1))
        m2p = m2.copy()
        s = product(m1, m2, m, m1p, m2p, mp)
        keep_list = []
        for i in s:
            m1, m2, m, m1p, m2p, mp = i[0], i[1], i[2], i[3], i[4], i[5]
            # Check input parameter conditions
            # Condition 1 & 2 & 5
            if not (abs(j1 - j2) <= j <= j1 + j2 and m1 + m2 == m
                    and j1 + j2 - j % 1 != 0.5):
                pass
            if not (abs(j1 - j2) <= j <= j1 + j2 and m1p + m2p == mp
                    and j1 + j2 - j % 1 != 0.5):
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








