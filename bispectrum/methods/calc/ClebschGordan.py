import numpy as np
import cmath

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
    Returns: numbers
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
        #Conditions
        J = j1 + j2 + j
        # Condition 1
        if not ((abs(j1 - j2) <= j <= j1 + j2) or isinstance(J, int) or (J >= 0)):
            raise ValueError("|𝑗1−𝑗2| ≤𝑗≤𝑗1+𝑗2 and 𝑗1+𝑗2+𝑗 are non-negative integer and 𝑗1+𝑗2=𝑗")
        # Condition 2
        if not ((isinstance(val, int) or val >= 0 for val in [j + j1 - j2, j - j1 + j2, j1 + j2 - j])):
            raise ValueError("(𝑗+𝑗1−𝑗2) and (𝑗−𝑗1+𝑗2) and (𝑗1+𝑗2−𝑗) and is non-negative integer")
        # Condition 3
        if j1 > J or j2 > J or j > J:
            raise ValueError("𝑗1,𝑗2,𝑗 not exceed a positive integer 𝐽=𝑗1+𝑗2+𝑗")
        # Condition 4
        if not (isinstance(val, (int, float)) and val % 0.5 == 0 for val in [m1, m2, m]):
            raise ValueError("𝑚1,𝑚2,𝑚 must be integer or half-integer (positive or negative) numbers")
        # Condition 5
        if not (isinstance(vals, int) and vals >= 0 for vals in [j1 + m1, j1 - m1, j2 + m2, j2 - m2, j + m, j - m]):
            raise ValueError("𝑗1+𝑚1,𝑗1−𝑚1, 𝑗2+𝑚2,𝑗2−𝑚2 𝑗+𝑚,𝑗−𝑚 are non-negative integer")
        # Condition 6
        if not m1 + m2 == m:
            raise ValueError("𝑚1+𝑚2=𝑚 and 𝑚1′+𝑚2′=𝑚′")
        # Condition 7
        if not (abs(val) <= limit for val, limit in [(m1, j1), (m2, j2), (m, j)]):
            raise ValueError("|𝑚1|≤𝑗1, |𝑚2|≤𝑗2, |𝑚|≤𝑗")
        # Condition 8
        if not (j2 + j + m1 >= 0 and j1 - j2 - m >= 0 and isinstance(j2 + j + m1, int) and isinstance(j1 - j2 - m, int)):
            raise ValueError("𝑗2+𝑗+𝑚1≥0 and 𝑗1−𝑗2−𝑚≥0 and 𝑗2+𝑗+𝑚1 and 𝑗1−𝑗2−𝑚 are non-negative integer")
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
    [4] Meremianin,
        Multipole expansions in four-dimensional hyperspherical harmonics (2006
    [5] Chapter 8  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    '''
    CG = Clebsch_Gordan(j1,j2,j,m1,m2,m)
    cg = CG.cg()
    CGp = Clebsch_Gordan(j1,j2,j,m1p,m2p,mp)
    cg_p = CGp.cg()
    H = (cg)*(cg_p)
    return H