"""
Inprogress
"""
import numpy as np
import cmath
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd

class Bispectrum:
    """
    Calculate bispectrum components
    """
    def __init__(self, j1, j2, j,rcut,input_file):
        self.j = j
        self.m = m
        self.mp = mp

    def neighbor_ID(input_file,center_atom_id, rcut):
        """
        This function is to create a neighbor_list for a chosen center atoms
        Parameters:
            input_file (string) directory of the input file
            center_atom_id (integer): center atom ID
            rcut (float): cutoff radius depends on the size of the unit cell and atom type
                    note:choosen rcut needs to divide by the true cell length
                        since atom coordinate (x,y,z)are fraction with cell dimension (1,1,1),
        Returns: neighbor atoms data frame with atom ID and distance from center
                    atom, theta_0, theta, phi
        """
        path = input_file
        dico = MMCIF2Dict(path)
        df_cif = pd.DataFrame.from_dict(dico, orient='index')
        x = df_cif.iloc[-3]
        y = df_cif.iloc[-2]
        z = df_cif.iloc[-1]







