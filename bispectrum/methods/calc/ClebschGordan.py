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
        #Condition 1 & 2
        if j1 + j2 < j or np.abs(j1 - j2) > j or m1 + m2 != m or j1+j2!=j:
            raise ValueError("Invalid input parameters: j1, j2, j, m1, m2, and m must satisfy the triangle inequality.")
        # Condition 3
        if abs(m1) > j1 or abs(m2) > j2 or abs(m) > j:
            raise ValueError("Invalid input parameters: |m1| <= j1, |m_2| <= j2, |m| <= j")
        # Condition 4
        J =(j1+j2+j)
        if J < (int(j1+j2+j)) and J <0:
            raise ValueError("Invalid input parameters: j1, j2, j must not exceed a positive integer J")
        # Condition 5
        if j1 + j2 - j % 1 == 0.5:
            raise ValueError("Invalid input parameters: j1+j2-j must not be a half-integer")
        # Condition 6
        if not all(isinstance(x, (int, float, np.integer, np.floating)) and (x % 0.5 == 0 or x % 1 == 0) for x in
                   [m1, m2, m]):
            raise ValueError(
                "Invalid input parameters: m1, m2, m must be integer or half-integer (positive or negative) numbers")
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


    def cb(self):
        if self.m1 + self.m2 != self.m:
            return 0.0  # delta function fails
        prefactor = cmath.sqrt((2 * self.j + 1) * fact(self.j + self.j1 - self.j2) * fact(self.j-self.j1 + self.j2) \
                               * fact(self.j1 + self.j2 - self.j) / fact(self.j + self.j1 + self.j2 + 1))
        coefficient = cmath.sqrt(fact(self.j + self.m) * fact(self.j - self.m) / (fact(self.j1 + self.m1) \
                               * fact(self.j1 - self.m1) * fact(self.j2 + self.m2) * fact(self.j2 - self.m2)))
        sum = 0.0
        s_min= max(0, int(self.m1-self.j1,int(self.j2-self.j1+self.m))
        s_max= min(int(self.j1+self.j+self.m1),int(self.j-self.j1+self.j2),int(j+m))
        for s in range(s_min, s_max + 1):
            den = fact(s) * fact(self.j - self.j1 + self.j2 - s) * fact(self.j + self.m - s) \
                  * fact(self.j1 + self.j2 - self.j - self.m + self.s)
            num = ((-1) ** (self.j2 + self.m2 + s))* fact(self.j2 + self.j + self.m1 - s) * fact(self.j1 - self.m1 + s)
            sum += num / den
        cb = prefactor * coefficient * sum
        return cb