"""
This file is used to study the cases example of SO(4) bispectrum calculation.
A compound where the central atom has only two symmetrical neighbors of the same atom type is known as a
diatomic molecule. In diatomic molecules, the central atom is bonded to two identical atoms.
"""
import numpy as np
import matplotlib.pyplot as plt
from bispectrum.methods.calc.SO4 import Bispectrum
from bispectrum.methods.calc.INPUT_param import test_neighbor
from bispectrum.methods.plot.neighbors import plot_atoms
#Case 1-Special case: 1 center atom in an empty cell
center_atom = (0.5, 0.5, 0.5)
cell_length = 1.0
R_cut = 0.2
j,j1,j2= 5/2, 3/2, 1
center_atom = [0.5, 0.5, 0.5]
cell_length = 1.0
R_cut = 0.2
j,j1,j2= 5/2, 3/2, 1
data_1 = {
        'w_ik': np.array([1]),
        'delta': np.array([0]),
        'r_ik': np.array([0]),
        'r_cut': np.array([0.2]),
        'theta_0': np.array([0]),
        'theta': np.array([0]),
        'phi': np.array([0])
    }
#bispectrum calculation
B_1= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_1)
B_data_1 = B_1.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_1)
print(B_data_1)
#plot
case_1_plot="case_1_plot.png"
plot_atoms(center_atom, None, R_cut, cell_length, B_data_1)