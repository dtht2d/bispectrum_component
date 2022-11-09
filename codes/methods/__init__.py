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
        j1: angular momentum 1
        j2: angular momentum 2
        j: total angular momentum (j1+j2)
        m1: eigenvalue of angular momentum j1
        m2: eigenvalue of angular momentum j2
        m: eigenvalue of angular momentum j
        m1p: eigenvalue of j1 along rotated axis
        m2p: eigenvalue of j2 along rotated axis
        mp: eigenvalue of j along rotated axis
    Returns:
        Coupling coefficient H(j1,j2,j,m1,m2,m, m1p,m2p,mp)
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
        j: integer/half integer number
           Total angular momentum
        m: integer/half integer number
           Eigenvalue of angular momentum along rotated axis
        mp:
        mpp:
        theta_0: fist angle of rotation [0,pi]
        theta: second angle of rotation [0,pi]
        phi: third angle of rotation [0,2pi]
        rotational od a coordinate system through an angle theta_0
        about an axis n(theta,phi)
    Returns: A single rotational matrix U function
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
        Um_mp = nsimplify(Um_mp)
        U += N(Um_mp) #N() evaluated U expression
    return U
def getCutoffFunction (r_ik,r_min0, R_cut):
    '''
    Parameter:
        r_ik: 1D array dimension [k,], distance from center atom  to k neighbor atoms
        r_min0: parameter in distance to angle conversion(distance unit)
        R_cut: cut off radius
    Returns:
        array, dimension [k,]: cut of function f_cut(r_ik) where r_ik < R_cut
    ==========================References==================================
    [3] A. Bartok, M. C. Payne, K. Risi, G. Csanyi, Gaussian approximation
        potentials: the accuracy of quantum mechanics, without the electrons (2010)
    [7] LAMMPS, Compute sna/atom command,
        https://docs.lammps.org/compute_sna_atom.html
    '''
    R_cut_array = np.full((r_ik_array.shape), R_cut)
    f_cut = (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array - r_min0, R_cut_array - r_min0))) + 1)
    return f_cut
def getDensityFunction_u(j,m,mp,w_ik_arr, delta_arr,r_ik_array, r_min0, R_cut, theta_0_array, theta_array, phi_array):
    '''
    Args:
        j: angular momentum
        m: eigenvalue of angular momentum
        mp: eigenvalue of j along rotated axis
        w_ik_array: array for the coefficients that are dimensionless weights that are chosen to distinguish atoms
                    of different types, while the central atom is arbitrarily assigned a unit weight, dimensin (1,k)
        delta_array: array for the Dirac delta function, indicates only neighbor atom of element the same as center atom
                    contribute to partial density,  dimension (1,k)
        r_ik_array: array for distance from center atom to neighbor atom, dimension (1,k+1), k is number of neighbor atoms
                    in cutoff radius, array include center atom as well
        r_min0: number, parameter in distance to angle conversion (distance units), choose
        R_cut: number, cutoff radius
        theta_0_array: array for theta_0 angel (fist angle of rotation [0,pi])
                    of neighbor atoms in reference frame of center atom, dimension (1, k+1)
        theta_array: array for theta angel ( second angle of rotation [0,pi])
                    of neighbor atoms in reference frame of center atom, dimension (1, k+1)
        phi_array: array for phi angel (third angle of rotation [0,2pi])
                    of neighbor atoms in reference frame of center atom, dimension (1, k+1)
    Returns: expansion coefficients density function u_jmmp
    '''
    R_cut_array = np.full((r_ik_array.shape), R_cut)
    f_cut_arr = (1 / 2) * (np.cos(np.pi * (np.divide(r_ik_array - r_min0, R_cut_array - r_min0))) + 1)
    U_jmmp = []
    for theta_0, theta, phi in zip(theta_0_array, theta_array, phi_array):
        mpp_list = np.linspace(-j, j, int(2 * j + 1))
        U = 0
        for mpp in mpp_list:
            rot1 = Rotation.D(j, m, mpp, phi, theta, - phi)
            rot2 = Rotation.D(j, mpp, mp, phi, -theta, -phi)
            Dj_m_mpp = rot1.doit()
            Dj_mpp_mp = rot2.doit()
            Um_mp = Dj_m_mpp * (exp(-I * mpp * theta_0)) * Dj_mpp_mp
            Um_mp = nsimplify(Um_mp)
            U += N(Um_mp)
        U_jmmp.append(U)
    U_jmmp_array = np.array(U_jmmp, dtype='complex')
    U_jmmp_arr = np.where(np.isnan(U_jmmp_array), 0, U_jmmp_array)
    u_jmmp = np.dot((f_cut_arr * U_jmmp_arr), (w_ik_arr * delta_arr))
    return u_jmmp
