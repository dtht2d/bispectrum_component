import numpy as np
import cmath

def fact(n):
    """
    This function is used to calculate factorial of a number by using
    an iterative approach instead of recursive approach
    """
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result


class Wigner_D:
    """
    Args:
        j (scalar): angular momentum
        m (scalar): eigenvalue of angular momentum
        mp (scalar): eigenvalue of j along rotated axis
        theta_0 (scalar): first angle of rotation [0, pi]
        theta (scalar): second angle of rotation [0, pi]
        phi (scalar): third angle of rotation [0, 2*pi]
    Returns: complex number, Wigner D function
    ==========================Reference==================================
    [5] Chapter 4.3-(p.76,eq.1)  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    """
    def __init__(self, j, m, mp, theta_0, theta, phi):
        if j < 0 or not np.isclose(j, int(j)) or (j % 1 == 0.5 and (m % 1 != 0 or mp % 1 != 0)):
            raise ValueError("Invalid input parameters: j must be a non-negative integer or half-integer, "
                             "m and mp must be between -j and j.")
        if theta_0 < 0 or theta_0 > 2*np.pi or theta < 0 or theta > np.pi or phi < 0 or phi > 2 * np.pi:
            raise ValueError(
                "Invalid input parameters: theta_0, theta, and phi must be within [0, pi] and [0, 2pi], respectively.")
        self.j = j
        self.m = m
        self.mp = mp
        self.theta_0 = theta_0
        self.theta = theta
        self.phi = phi
    @property
    def compute_dsmall(self):
        """
        This method is used to calculate the Wigner d small- real function involving trigonometric functions
        ==========================Reference==================================
        [5] Chapter 4.3.1-(p.76,eq.4)  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Returns: Wigner d small (theta) - real function
        """
        kmax = max(0, self.m - self.mp)
        kmin = min(self.j + self.m, self.j - self.mp)
        term1 = np.sqrt(fact(self.j + self.m) * fact(self.j - self.m) * fact(self.j + self.mp) * fact(self.j - self.mp))
        sum = 0
        for k in range(kmax, kmin + 1):
            numerator = (-1) ** k * (cmath.cos(self.theta / 2)) ** (2 * self.j - 2 * k + self.m - self.mp) * \
                        (cmath.sin(self.theta / 2)) ** (2 * k - self.m + self.mp)
            denominator = fact(k) * fact(self.j + self.m - k) * fact(self.j - self.mp - k) * fact(self.mp - self.m + k)
            sum += numerator / denominator
        return sum*term1

    def wigner_D(self):
        # term1 = np.exp(-1j * self.m * self.theta_0)
        term1 = np.cos(self.m * self.theta_0) - 1j * (np.sin(self.m * self.theta_0))
        term2 = self.compute_dsmall
        # term3 = np.exp(-1j * self.mp * self.phi)
        term3 = np.cos(self.mp * self.phi) - 1j * (np.sin(self.mp * self.phi))
        result = term1 * term2 * term3
        return result




