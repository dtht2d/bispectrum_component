import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_d
from sympy.physics.quantum.spin import Rotation
from sympy import *

def getCoeffH(j1,j2,j,m1,m2,m,m1p,m2p,mp):
    '''
    This function calculate coupling coefficient H via computing
    the Clebsch-Gordan coefficient for cg(j1,m1,j2,m2,j,m)
    and cg_p(j1,m1p,j2,m2p,j,mp)
    Parameters:
        j1 (scalar): angular momentum 1
        j2 (scalar): angular momentum 2
        j (scalar): total angular momentum (j1+j2)
        m1 (scalar): eigenvalue of angular momentum j1
        m2 (scalar): eigenvalue of angular momentum j2
        m (scalar): eigenvalue of angular momentum j
        m1p (scalar): eigenvalue of j1 along rotated axis
        m2p (scalar): eigenvalue of j2 along rotated axis
        mp (scalar): eigenvalue of j along rotated axis
    Returns:
        Scalar, Coupling coefficient H(j1,j2,j,m1,m2,m, m1p,m2p,mp)
    ====================Reference=========================
    [1] Thompson, Swiler, Trott, Foiles, Tucker,
        Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials (2015)
    [5] Chapter 8  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    '''
    cg = CG(j1,m1,j2,m2,j,m)
    cg = cg.doit()
    cg_p = CG(j1,m1p,j2,m2p,j,mp)
    cg_p = cg_p.doit()
    H_coeff = (cg)*(cg_p)
    H = N(H_coeff)
    return H
def getRotationalMatrixU(j,m,mp,theta_0,theta,phi):
    '''
    Parameters:
        j (scalar): integer/half integer number
           Total angular momentum
        m (scalar): integer/half integer number
           eigenvalue of angular momentum along rotated axis
        mp (scalar): eigenvalue of j along rotated axis
        mpp (scalar):
        theta_0 (scalar): fist angle of rotation [0,pi]
        theta (scalar): second angle of rotation [0,pi]
        phi (scalar): third angle of rotation [0,2pi]
        rotational od a coordinate system through an angle theta_0
        about an axis n(theta,phi)
    Returns: complec number, a single rotational matrix U function
    ==========================Reference==================================
    [5] Chapter 4  D.A. Varshalovich, A.N. Moskalev, V.K Khersonskii,
        Quantum Theory of Angular Momentum (1988)
    '''
    mpp_list = np.linspace(-j, j, int(2 * j + 1))
    U = 0
    for mpp in mpp_list:
        rot1 = Rotation.D(j, m, mpp, phi, theta, - phi)
        rot2 = Rotation.D(j, mpp, mp, phi, -theta, -phi)
        Dj_m_mpp = rot1.doit()
        Dj_mpp_mp = rot2.doit()
        Um_mp = Dj_m_mpp * (exp(-I * mpp * theta_0)) * Dj_mpp_mp
        U += N(Um_mp) #N() evaluated U expression
    return U
def getCutoffFunction (r_ik,r_min0, R_cut):
    '''
    Parameter:
        r_ik (ndarray): shape (k+1,), distance from center atom  to k neighbor atoms
        r_min0 (scalar): parameter in distance to angle conversion(distance unit)
        R_cut (scalar): cut off radius
    Returns:
        ndarray, dimension [k,]: cut of function f_cut(r_ik) where r_ik < R_cut
    ==========================References==================================
    [3] A. Bartok, M. C. Payne, K. Risi, G. Csanyi, Gaussian approximation
        potentials: the accuracy of quantum mechanics, without the electrons (2010)
    [7] LAMMPS, Compute sna/atom command,
        https://docs.lammps.org/compute_sna_atom.html
    '''
    R_cut_array = np.full((r_ik_array.shape), R_cut) #R_cut_array shape (k,)
    f_cut_arr = (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array - r_min0, R_cut_array - r_min0))) + 1)
    return f_cut_arr
def getDensityFunction_u(j,m,mp,w_ik_arr, delta_arr,r_ik_array, r_min0, R_cut, theta_0_array, theta_array, phi_array):
    '''
    Args:
        j (scalar): angular momentum
        m (scalar): eigenvalue of angular momentum
        mp (scalar): eigenvalue of j along rotated axis
        w_ik_arr (ndarray): shape(k+1,) the coefficients that are dimensionless weights that are chosen
                            to distinguish atoms of different types, while the central atom is arbitrarily
                            assigned a unit weight
        delta_arr (ndarray): shape(k+1,) the Dirac delta function, indicates only neighbor atom of element
                             the same as center atom contribute to partial density
        r_ik_array (ndarray): shape(k+1) distance from center atom to neighbor atom, k is number of neighbor atoms
                              in cutoff radius, array include center atom as well
        r_min0 (scalar): parameter in distance to angle conversion (distance units), choose
        R_cut (scalar): cutoff radius
        theta_0_array (ndarray): shape(k+1,) theta_0 angel (fist angle of rotation [0,pi])
                                 of neighbor atoms in reference frame of center atom
        theta_array (ndarray): shape (k+1,) theta angel ( second angle of rotation [0,pi])
                               of neighbor atoms in reference frame of center atom
        phi_array (ndarray): shape(k+1,) phi angel (third angle of rotation [0,2pi])
                             of neighbor atoms in reference frame of center atom
    Returns: complex number, expansion coefficients density function u_jmmp
    '''
    R_cut_array = np.full((r_ik_array.shape), R_cut)
    f_cut_arr = (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array - r_min0, R_cut_array - r_min0))) + 1)
    U_jmmp = []
    for theta_0, theta, phi in zip(theta_0_array, theta_array, phi_array):
        U = getRotationalMatrixU(j,m,mp,theta_0,theta,phi) #U is a complex number
        U_jmmp.append(U)
    U_jmmp_array = np.array(U_jmmp, dtype='complex')
    U_jmmp_arr = np.where(np.isnan(U_jmmp_array), 0, U_jmmp_array)
    u_jmmp = np.dot((f_cut_arr * U_jmmp_arr), (w_ik_arr * delta_arr))
    return u_jmmp

