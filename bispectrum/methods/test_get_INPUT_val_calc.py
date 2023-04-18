from bispectrum.methods.calc.get_INPUT_val import get_INPUT_value
import pandas as pd
center_atom_id = 17
r_mu =0.0779
R_cut= 0.25
input_file = "/Users/duonghoang/Documents/GitHub/bispectrum_component/data/avgBL-Model.cif"
output_dir = "/Users/duonghoang/Documents/GitHub/bispectrum_component/data"
file_type ='cif'
input_data = get_INPUT_value(center_atom_id, r_mu, R_cut, input_file, output_dir, file_type)

print (input_data)
# Convert dictionary to DataFrame
df = pd.DataFrame(input_data)
# Display DataFrame
print(df)