from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import numpy as np
from bispectrum.methods.calc.Parameters import neighbor_atoms
input_file = "/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif"
center_atom_id = 17
rcut = 0.15
df_ik = neighbor_atoms(input_file, center_atom_id, rcut)
print(df_ik)