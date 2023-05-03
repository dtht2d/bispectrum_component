import numpy as np
import json
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from bispectrum.methods.calc.ClebschGordan import Clebsch_Gordan, H_coeff
from bispectrum.methods.calc.WignerD import Wigner_D, U_rot, u
from typing import Dict

fact_cache = {}
def fact(n):
    if n < 0 or not np.isclose(n, int(n)):
        raise ValueError("Invalid input parameter: n must be a non-negative integer.")
    if n in fact_cache:
        return fact_cache[n]
    result = np.math.factorial(int(n))
    fact_cache[n] = result
    return result
class Clebsch_Gordan:
    """
    Definition:
    A Clebsch-Gordan coefficients are vector addition coefficients. They play an important role in decomposition of
    reducible representations of rotation. Let j1 and j2 with projections on m1 and m2 on the quantization axis.
    The coefficients represent the probability amplitude that j1 and j2 are coupled into a resultant angular momentum
    j with projection m.
    Args:
        j1 (scalar): angular momentum
        j2 (scalar): angular momentum
        j (scalar): angular momentum
        m1 (scalar): eigenvalue of angular momentum
        m2 (scalar): eigenvalue of angular momentum
        m (scalar): eigenvalue of angular momentum
    Returns: complex number, Clebsh Gordan function
    ==========================Reference==================================
    [5] Chapter 8 D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    [12] Chapter 3 Biedenharn, L., & Louck, J.D. ,
        Encyclopedia of Mathematics and its Applications (1981)
    """
    def __init__(self, j1, j2, j, m1, m2, m):
        self.j1 = j1
        self.j2 = j2
        self.j = j
        self.m1 = m1
        self.m2 = m2
        self.m = m
        self.J = j1 + j2 + j
        if self.m1 + self.m2 != self.m:
            return 0.0  # delta function fails
        prefactor = cmath.sqrt((2 * self.j + 1) * fact(self.j + self.j1 - self.j2) * fact(self.j-self.j1 + self.j2) \
                               * fact(self.j1 + self.j2 - self.j) / fact(self.j + self.j1 + self.j2 + 1))
        coefficient = cmath.sqrt(fact(self.j + self.m) * fact(self.j - self.m) / (fact(self.j1 + self.m1) \
                               * fact(self.j1 - self.m1) * fact(self.j2 + self.m2) * fact(self.j2 - self.m2)))
        sum = 0.0
        smin= max(0, int(self.m1-self.j1),int(self.j2-self.j1+self.m))
        smax= min(int(self.j2+self.j+self.m1),int(self.j-self.j1+self.j2),\
                   int(self.j+self.m))

        for s in range(smin,smax+1):
            den = fact(s) * fact(self.j - self.j1 + self.j2 - s) * fact(self.j + self.m - s) \
                  * fact(self.j1 - self.j2 - self.m + s)
            num = ((-1) ** (self.j2 + self.m2 + s))* fact(self.j2 + self.j + self.m1 - s) * fact(self.j1 - self.m1 + s)
            sum += num / den
        cg = prefactor * coefficient * sum
        return cg

def H_coeff(j1,j2,j,m1,m2,m,m1p,m2p,mp):
    '''
    This function calculate coupling coefficient H via computing
    the Clebsch-Gordan coefficient for cg(j1,m1,j2,m2,j,m)
    and cg(j1,m1p,j2,m2p,j,mp)
    Parameters:
        j1: angular momentum 1
        j2: angular momentum 2
        j: total angular momentum (j1+j2)
        m1: eigenvalue of angular momentum j1
        m2: eigenvalue of angular momentum j2
        m: eigenvalue of angular momentum j
        m1p: eigenvalue of j1 along rotated axis
        m2p: eigenvalue of j2 along rotated axis
        mp: eigenvalue of j along rotated axis
    Returns: Coupling coefficient H(j1,j2,j,m1,m2,m.m1p,m2p,mp)
    ======================Reference=========================
    [1] Thompson, Swiler, Trott, Foiles, Tucker,
        Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials (2015)
    [5] Chapter 8  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    '''
    CG = Clebsch_Gordan(j1,j2,j,m1,m2,m)
    cg = CG.cg()
    CGp = Clebsch_Gordan(j1,j2,j,m1p,m2p,mp)
    cg_p = CGp.cg()
    H = (cg)*(cg_p)
    return H
class B:
    """
    Calculate bispectrum- S(0)4 components
    """
    def __init__(self, j: float, j1: float, j2: float, params: Dict[str, np.ndarray]):
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
        self.params = params
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

    @staticmethod
    def generate_m_values(j, j1, j2):
        """
        This function generates (m1, m2, m, m1p, m2p, mp) from input set (j1,j2,j)
        and only keeps sets that satisfy the conditions.
        Return: a one and only unique set (m1, m2, m, m1p, m2p, mp)
        """
        J = j1 + j2 + j
        # Condition 1
        if not ((abs(j1 - j2) <= j <= j1 + j2) or isinstance(J, int) or (J >= 0) or (j1 + j2) == j):
            raise ValueError("Condition (1) is not satisfied.")
        # Condition 2
        if not ((isinstance(val, int) or val >= 0 for val in [j + j1 - j2, j - j1 + j2, j1 + j2 - j])):
            raise ValueError("(𝑗+𝑗1−𝑗2) and (𝑗−𝑗1+𝑗2) and (𝑗1+𝑗2−𝑗) and is non-negative integer")
        # Condition 3
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
        keep_list = []
        for i in range(len(s)):
            m1_val, m2_val, m_val, m1p_val, m2p_val, mp_val = s[i]
            # Condition 4-8: Clebsch-Gordan calc for set (𝑗1,𝑗2,𝑗,𝑚1,𝑚2,𝑚), (𝑗1,𝑗2,𝑗,𝑚1p,𝑚2p,𝑚p)
            c4 = [isinstance(val, (int, float)) and val % 0.5 == 0 for val in
                  [m1_val, m2_val, m_val, m1p_val, m2p_val, mp_val]]
            c5 = [isinstance(vals, int) and vals >= 0 for vals in
                  [j1 + m1_val, j1 - m1_val, j2 + m2_val, j2 - m2_val, j + m_val, j - m_val]]
            c5_p = [isinstance(vals, int) and vals >= 0 for vals in
                    [j1 + m1p_val, j1 - m1p_val, j2 + m2p_val, j2 - m2p_val, j + mp_val, j - mp_val]]
            c6 = [m1_val + m2_val == m_val and m1p_val + m2p_val == mp_val]
            c7 = [abs(val) <= limit for val, limit in
                  [(m1_val, j1), (m2_val, j2), (m_val, j), (m1p_val, j1), (m2p_val, j2), (mp_val, j)]]
            c8 = [j2 + j + m1_val >= 0 and j1 - j2 - m_val >= 0 and isinstance(j2 + j + m1_val, int) and isinstance(
                j1 - j2 - m_val, int)]
            c8_p = [
                j2 + j + m1p_val >= 0 and j1 - j2 - mp_val >= 0 and isinstance(j2 + j + m1p_val, int) and isinstance(
                    j1 - j2 - mp_val, int)]
            # Condition 9: Wigner-D calc
            c9 = [isinstance(vals, int) and vals >= 0 for vals in [mp_val - m_val, m1p_val - m1_val, m2p_val - m2_val]]
        if (c4) and (c5 and c5_p) and (c6) and (c7) and (c8 and c8_p) and (c9):
            keep_list.append(s[i])
        else:
            pass
        return keep_list
    #Clebsch-Gordan Coefficient
    @staticmethod
    def cg (j1, j2, j, m1, m2, m):
        if m1 + m2 != m:
            return 0.0  # delta function fails
        prefactor = np.sqrt((2 * j + 1) * fact(j + j1 - j2) * fact(j - j1 + j2) * fact(j1 + j2 - j) / fact(j + j1 + j2 + 1))
        coefficient = np.sqrt(fact(j + m) * fact(j - m) / (fact(j1 + m1) * fact(j1 - m1) * fact(j2 + m2) * fact(j2 - m2)))
        sum = 0.0
        smin = max(0, int(m1 - j1), int(j2 - j1 + m))
        smax = min(int(j2 + j + m1), int(j - j1 + j2), int(j +m))
        for s in range(smin, smax + 1):
            den = fact(s) * fact(j - j1 + j2 - s) * fact(j + m - s) * fact(j1 - j2 - m + s)
            num = ((-1) ** (j2 + m2 + s)) * fact(j2 + j + m1 - s) * fact(j1 - m1 + s)
            sum += num / den
        cg = prefactor * coefficient * sum
        return cg










