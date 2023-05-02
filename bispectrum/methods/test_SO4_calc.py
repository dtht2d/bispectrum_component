
from bispectrum.methods.calc.SO4 import B, generate_m_val
from bispectrum.methods.calc.INPUT_param import get_INPUT_value
path = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif'
center_atom_id = 17
r_mu = 0.0779
R_cut = 0.25
input_file_path = path
output_directory = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data'
#read data from file
data = get_INPUT_value(center_atom_id, r_mu, R_cut, input_file_path, output_directory, file_type='cif')

#possible list of m values
B = B(j,j1,j2, params=data)
list_generate = B.generate_m_val(1,2,3)
print (list_generate)