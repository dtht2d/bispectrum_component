"""
Same-cell: Plot the neighbor atoms of a center atom in 1D.
"""
import matplotlib.pyplot as plt
import numpy as np

def plot_atoms(center_atom, neighbor_atoms, R_cut, cell_length, B, save_path=None):
    # Plot
    fig, ax = plt.subplots()
    ax.set_xlim([0, cell_length])
    ax.set_ylim([0, cell_length])
    ax.set_aspect('equal')

    ax.scatter(center_atom[0], center_atom[1], color='red', label='Center Atom')

    if neighbor_atoms is not None:
        for neighbor_atom in neighbor_atoms:
            ax.scatter(neighbor_atom[0], neighbor_atom[1], color='blue', label='Neighbor Atom')

    circle = plt.Circle((center_atom[0], center_atom[1]), R_cut, fill=False, color='green')
    ax.add_artist(circle)

    ax.annotate(r'$R_{\mathrm{cut}}$', xy=(center_atom[0] + R_cut, center_atom[1]), xytext=(center_atom[0] + R_cut, center_atom[1] + 0.1),
                 arrowprops=dict(facecolor='blue', arrowstyle='->'))
    ax.set_title(r'$B_{{j_1,j_2,j}}={}$'.format(np.round(B, 2)))
    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()
