
"""
Function to calculate the Wigner D function, Rotational Matrix U, and expansion coefficients density function u
"""
import numpy as np
fact_cache = {}

def fact(n):
    if n < 0 or not np.isclose(n, int(n)):
        raise ValueError("Invalid input parameter: n must be a non-negative integer.")
    if n in fact_cache:
        return fact_cache[n]
    result = np.math.factorial(int(n))
    fact_cache[n] = result
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
        #Conditions for j, m, mp, theta_0, theta, phi
        #𝑗+𝑚, 𝑗−𝑚, 𝑗+𝑚′, 𝑗−𝑚′ are non-negative integers
        if not (isinstance(vals, int) and vals >= 0 for vals in [j + m, j - m, j + mp, j - mp]):
            raise ValueError("Invalid input parameters: j+m, j-m, j+m', and j-m' must be non-negative integers")
        #|𝑚|≤𝑗, |𝑚′|≤𝑗
        if not ([abs(val) <= limit for val, limit in [(m, j), (mp, j)]]):
          raise ValueError("Invalid input parameters: |m| and |mp| must be less than or equal to j")
        if not ((j >= 0) or ((j % 1 == 0.0) or (j % 1 == 0.5))):
          raise ValueError("Invalid input parameters: j must be a non-negative integer or half-integer")
        if theta_0 < 0 or theta_0 > 2*np.pi or abs(theta) > np.pi or abs(phi) > 2 * np.pi:
          raise ValueError("Invalid input parameters: theta_0, theta, and phi must be within [0, pi] and [0, 2pi], respectively.")
        self.j = j
        self.m = m
        self.mp = mp
        self.theta_0 = theta_0
        self.theta = theta
        self.phi = phi
    def compute_dsmall(self):
        """
        This method is used to calculate the Wigner d small- real function involving trigonometric functions
        ==========================Reference==================================
        [5] Chapter 4.3.1-(p.76,eq.4)  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Returns: Wigner d - real function
        """
        kmax = max(0, self.m - self.mp)
        kmin = min(self.j + self.m, self.j - self.mp)
        term1 = np.sqrt(fact(self.j + self.m) * fact(self.j - self.m) * fact(self.j + self.mp) * fact(self.j - self.mp))
        sum = 0
        for k in range(int(kmax), int(kmin) + 1):
            numerator = (-1) ** k * (np.cos(self.theta / 2)) ** (2 * self.j - 2 * k + self.m - self.mp) * \
                        (np.sin(self.theta / 2)) ** (2 * k - self.m + self.mp)
            denominator = fact(k) * fact(self.j + self.m - k) * fact(self.j - self.mp - k) * fact(self.mp - self.m + k)
            sum += numerator / denominator
        return sum*term1
    def wigner_D(self):
        #term1 = np.exp(-1j * self.m * self.theta_0)
        term1 = np.cos(self.m * self.theta_0) - 1j*(np.sin(self.m * self.theta_0))
        term2 = self.compute_dsmall()
        #term3 = np.exp(-1j * self.mp * self.phi)
        term3 = np.cos(self.mp * self.phi) -1j*(np.sin(self.mp * self.phi))

        result = term1 * term2 * term3
        return result
def U_rot(j, m, mp, theta_0, theta, phi):
    '''
    Parameters:
      j: integer/half integer number
         angular momentum
      m: integer/half integer number
         Eigenvalue of angular momentum along rotated axis
      mp: eigenvalue of j along rotated axis
      mpp:
      theta_0: fist angle of rotation [0,pi]
      theta: second angle of rotation [0,pi]
      phi: third angle of rotation [0,2pi]
           rotational od a coordinate system through an angle theta_0
          about an axis n(theta,phi)
      Returns: Rotational matrix U
      ==========================Reference==================================
      [5] Chapter 4.5.2 (a) Eq. (3)  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
              Quantum Theory of Angular Momentum (1988)
    '''
    mpp_vals = np.linspace(-j, j, int(2*j+1))
    U = 0
    for mpp in mpp_vals:
        WD_1 = Wigner_D(j, m, mpp, phi, theta,-phi)
        term1 = WD_1.wigner_D()
        term2 = np.cos(mpp * theta_0) - 1j*(np.sin(mpp * theta_0))
        WD_2 = Wigner_D(j, mpp, mp, phi, -theta, -phi)
        term3 = WD_2.wigner_D()
        Um_mp = term1 * term2 * term3
        U += Um_mp
    return U
def u(j, m, mp, params):
    '''
    Parameter:
        j (scalar): angular momentum
        m (scalar): eigenvalue of angular momentum
        mp (scalar): eigenvalue of j along rotated axis
        params (dict): a dictionary containing the following keys, its values:
            - w_ik (array): the coefficients that are dimensionless weights that are chosen to distinguish atoms
              of different types, while the central atom is arbitrarily assigned a unit weight, dimensin (1,k)
            - delta (array): the Dirac delta function, indicates only neighbor atom of element the same as center atom
              contribute to partial density,  dimension (1,k)
            - r_ik (array): distance from center atom to neighbor atom, dimension (1,k), k is number of neighbor atoms
              in cutoff radius, array exclude center atom as well
            - r_cut (array): cutoff radius
            - theta_0: array for theta_0 angel (fist angle of rotation [0,pi])
              of neighbor atoms in reference frame of center atom, dimension (k+1,)
            - theta: array for theta angel ( second angle of rotation [0,pi])
              of neighbor atoms in reference frame of center atom, dimension (k+1,)
            - phi: array for phi angel (third angle of rotation [0,2pi])
              of neighbor atoms in reference frame of center atom, dimension (k+1,)
    Returns: expansion coefficients density function u_jm_mp
    '''
    w_ik_array = params['w_ik']
    delta_array = params['delta']
    r_ik_array = params['r_ik']
    r_cut_array = params['r_cut']
    theta_0_array = params['theta_0']
    theta_array = params['theta']
    phi_array = params['phi']
    # Calculate cutoff_function
    mask = r_ik_array >= r_cut_array
    # Set elements of f_cut_arr to 0 where the mask is True
    f_cut_arr = np.where(mask, 0, (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array, r_cut_array))) + 1))
    # Calculate rotational matrix U for all k=n neighbor atoms
    U_ik_array = np.array([U_rot(j, m, mp, theta_0, theta, phi) for theta_0, theta, phi in zip(theta_0_array, theta_array, phi_array)], dtype='complex')
    # Compute u_jmmp
    u_jmmp = np.dot((f_cut_arr * U_ik_array), (w_ik_array * delta_array))
    return u_jmmp


