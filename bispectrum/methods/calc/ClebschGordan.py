import numpy as np
import cmath
def fact(n):
    """
    This function is used to calculate factorial of a number by using
    an iterative approach instead of recursive approach
    """
    return np.prod(np.arange(1, n + 1))
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
        #Condition 1 & 2 & 5
        if not (abs(j1 - j2) <= j <= j1 + j2 and m1 + m2 == m and j1 + j2 - j % 1 != 0.5):
            raise ValueError("Invalid input parameters: j1, j2, j, m1, m2, and m must satisfy the triangle inequality.\ "
                             "j1+j2-j must not be a half-integer")
        #Condition 3 & 6
        if not all(abs(x) <= y and (x % 0.5 == 0 or x % 1 == 0) for x, y in zip([m1, m2, m], [j1, j2, j])):
            raise ValueError("Invalid input parameters: |m1| <= j1, |m_2| <= j2, |m| <= j\ "
                             "m1, m2, m must be integer or half-integer (positive or negative) numbers")
        # Condition 4
        J =(j1+j2+j)
        if J < (int(j1+j2+j)) and J <0:
            raise ValueError("Invalid input parameters: j1, j2, j must not exceed a positive integer J")
        # Condition 7
        if not all(
                isinstance(x, (int, float, np.integer, np.floating)) and x >= 0 and (x % 0.5 == 0 or x % 1 == 0) for x
                in [j1, j2, j]):
            raise ValueError(
                "Invalid input parameters: j1, j2, j must be integer or half-integer non-negative numbers")
        # Condition 8
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and x >= 0 and x % 1 == 0 for x in
                   [j1 + m1, j2 + m2, j + m, j1 + j2 + j]):
            raise ValueError("Invalid input parameters: j1+m1,j2+m2,j+m,j1+j2+j must be integer non-negative numbers")
    def cg(self):
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
    return H, cg, cg_p