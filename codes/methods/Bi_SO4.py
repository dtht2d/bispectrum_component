
class Wigner_D (j,m,mp,theta_0,theta,phi):
    '''
    Args:
        j (scalar): angular momentum
        m (scalar): eigenvalue of angular momentum
        mp (scalar): eigenvalue of j along rotated axis
        theta_0 (scalar): fist angle of rotation [0,pi]
        theta (scalar): second angle of rotation [0,pi]
        phi (scalar): third angle of rotation [0,2pi]
    Returns: complex number, wigner D function
    ==========================Reference==================================
    [5] Chapter 4-(p.76,eq.1&3)  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    '''
    def __init__(self,j,m,mp,theta_0,theta,phi):
        self.j = j
        self.m = m
        self.mp = mp
        self.theta_0 = theta_0
        self.theta = theta
        self.phi = phi
    def fact(seln):
        if n == 0:
            return 1
        else:
            return n * self.factorial(n-1)