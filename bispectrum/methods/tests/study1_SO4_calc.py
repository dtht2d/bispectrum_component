"""
This file is used to study the cases example of SO(4) bispectrum calculation.
"""
import numpy as np
import matplotlib.pyplot as plt
from bispectrum.methods.calc.SO4 import Bispectrum
from bispectrum.methods.calc.INPUT_param import test_neighbor
from bispectrum.methods.plot.neighbors import plot_atoms

#Parameters stay constant for all cases
center_atom = (0.5, 0.5, 0.5)
cell_length = 1.0
R_cut = 0.2
j,j1,j2= 5/2, 3/2, 1

#Case 1-Special case: 1 center atom in an empty cell

data_1 = {
        'w_ik': np.array([1]),
        'delta': np.array([0]),
        'r_ik': np.array([0]),
        'r_cut': np.array([0.2]),
        'theta_0': np.array([0]),
        'theta': np.array([0]),
        'phi': np.array([0])
    }
B_1= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_1)
B_data_1 = B_1.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_1)
print(B_data_1)
case_1_plot=save_path = "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_1.png"
plot_atoms(center_atom, None, R_cut, cell_length, B_data_1, save_path=case_1_plot)

#Case 2-[Same cell at (1)] with 1 neighbor atom at (0.6, 0.5, 0.5)
center_atom = [0.5, 0.5, 0.5]
neighbor_atoms_2 = [(0.6, 0.5, 0.5)]
data_2 = test_neighbor(center_atom, neighbor_atoms_2, R_cut)
B_2= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_2)
B_data_2 = B_2.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_2)
print(B_data_2)
case_2_plot=save_path = "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_2.png"
plot_atoms(center_atom, neighbor_atoms_2, R_cut, cell_length, B= B_data_2, save_path=case_2_plot)

#Case 3-[Same cell at (1)] with 1 neighbor atoms at (0.6, 0.5, 0.5) and (0.5, 0.4, 0.5)
neighbor_atoms_3 = [(0.5,0.4,0.5)]
data_3 = test_neighbor(center_atom, neighbor_atoms_3, R_cut)
B_3= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_3)
B_data_3 = B_3.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_3)
print(B_data_3)
case_3_plot=save_path = "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_3.png"
plot_atoms(center_atom, neighbor_atoms_3, R_cut, cell_length, B= B_data_3, save_path=case_3_plot)

#Case 4-[Same cell at (1)] with 1 neighbor atom at (0.4, 0.5, 0.5)
neighbor_atoms_4 = [(0.4,0.5,0.5)]
data_4 = test_neighbor(center_atom, neighbor_atoms_4, R_cut)
B_4= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_4)
B_data_4 = B_4.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_4)
print(B_data_4)
case_4_plot=save_path = "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_4.png"
plot_atoms(center_atom, neighbor_atoms_4, R_cut, cell_length, B= B_data_4, save_path=case_4_plot)

#Case 5-[Same cell at (1)] with 1 neighbor atom at (0.5, 0.6, 0.5)
neighbor_atoms_5 = [(0.5,0.6,0.5)]
data_5 = test_neighbor(center_atom, neighbor_atoms_5, R_cut)
B_5= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_5)
B_data_5 = B_5.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_5)
print(B_data_5)
case_5_plot= "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_5.png"
plot_atoms(center_atom, neighbor_atoms_5, R_cut, cell_length, B= B_data_5, save_path=case_5_plot)

#Case 6-[Same cell at (1)] 2 neighbor atoms  combine case 2+3
neighbor_atoms_23 = [[0.6, 0.5, 0.5],[0.5,0.4,0.5]]
data_23 = test_neighbor(center_atom, neighbor_atoms_23, R_cut)
B_23= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_23)
B_data_23 = B_23.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_23)
print(B_data_23)
case_23_plot= "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_23.png"
plot_atoms(center_atom, neighbor_atoms_23, R_cut, cell_length, B= B_data_23, save_path=case_23_plot)

#Case 7-[Same cell at (1)] 2 neighbor atoms  combine case 4+5
neighbor_atoms_45 = [[0.4,0.5,0.5],[0.5,0.6,0.5]]
data_45 = test_neighbor(center_atom, neighbor_atoms_23, R_cut)
B_45= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_45)
B_data_45 = B_5.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_45)
print(B_data_45)
case_45_plot= "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_45.png"
plot_atoms(center_atom, neighbor_atoms_45, R_cut, cell_length, B= B_data_45, save_path=case_45_plot)


#Case 8-[Same cell at (1)] 3 neighbor atoms  combine case 2+3+4
neighbor_atoms_234 = [[0.6, 0.5, 0.5],[0.5,0.4,0.5], [0.4,0.5,0.5]]
data_234 = test_neighbor(center_atom, neighbor_atoms_234, R_cut)
B_234= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_234)
B_data_234 = B_234.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_234)
print(B_data_234)
case_234_plot= "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_234.png"
plot_atoms(center_atom, neighbor_atoms_234, R_cut, cell_length, B= B_data_234, save_path=case_234_plot)

#Case 9-[Same cell at (1)] 3 neighbor atoms  combine case 2+3+5
neighbor_atoms_235 = [[0.6, 0.5, 0.5],[0.5,0.4,0.5], [0.5, 0.6, 0.5]]
data_235 = test_neighbor(center_atom, neighbor_atoms_235, R_cut)
B_235= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_235)
B_data_235 = B_235.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_235)
print(B_data_235)
case_235_plot= "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_235.png"
plot_atoms(center_atom, neighbor_atoms_235, R_cut, cell_length, B= B_data_235, save_path=case_235_plot)


#Case 10-[Same cell at (1)] 3 neighbor atoms  combine case 3+4+5
neighbor_atoms_345 = [[0.5,0.4,0.5],[0.4, 0.5, 0.5], [0.5, 0.6, 0.5]]
data_345 = test_neighbor(center_atom, neighbor_atoms_345, R_cut)
B_345= Bispectrum(j=5/2, j1=3/2, j2=2/2, params=data_345)
B_data_345 = B_345.evaluate(j=5/2, j1=3/2, j2=2/2, params=data_345)
print(B_data_345)
case_345_plot= "/Users/duonghoang/Documents/GitHub/bispectrum_component/plots/plot_case_345.png"
plot_atoms(center_atom, neighbor_atoms_345, R_cut, cell_length, B= B_data_345, save_path=case_345_plot)
