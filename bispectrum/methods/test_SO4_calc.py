
from bispectrum.methods.calc.SO4 import Bispectrum, get_INPUT_value
path = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif'
center_atom_id = 17
r_mu = 0.0779
R_cut = 0.25
input_file_path = path
output_directory = '/Users/duonghoang/Documents/GitHub/bispectrum_component/data'
#read data from file
data = get_INPUT_value(center_atom_id, r_mu, R_cut, input_file_path, output_directory, file_type='cif')
print (data)
print (data.keys())
#crete dictionary mapping column name
dic_ID = {name: i for i, name in enumerate(data.keys())}
print (dic_ID)
print (type(data['phi']))