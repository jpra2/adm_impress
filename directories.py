import os
import yaml
from impress.preprocessor import directories as direc_impress

# entities_lv0 = ['nodes', 'edges', 'faces', 'volumes']
entities_lv0 = direc_impress.entities_lv0
entities_lv0_0 = direc_impress.entities_lv0_0

names_data_loaded_lv0 = ['read_permeability', 'file_name_permeability', 'Permeability']

names_data_loaded_lv2 = ['type', 'value', 'p0', 'p1']

types_region_data_loaded = ['all', 'box']

variables_impress = {'permeability': 'permeability', 'poro': 'poro', 'k_harm': 'k_harm',
                     'area': 'area', 'dist_cent': 'dist_cent', 'u_normal': 'u_normal'}

impress = 'impress'
tcc = 'tcc'
flying = 'flying'
adm = 'adm'
input_file_0 = 'inputs0.yml'
ext_h5m = '.h5m'
name_out_file = 'output'
states = ['0', '1', '2']
state_file = 'state.npy'

parent_dir = os.path.dirname(os.path.abspath(__file__)) # adm_impress_dir
adm_impress_dir = parent_dir
# path_ant = os.getcwd()

path_adm = os.path.join(adm_impress_dir, adm)
path_impress = direc_impress.impress_path

path_local_variables = direc_impress.path_local_variables
path_local_info_data = direc_impress.path_local_info_data

path_flying_geral = os.path.join(adm_impress_dir, flying)
output_file = os.path.join(flying, name_out_file)
state_path = os.path.join(flying, state_file)
last_file_name = 'last_file_name.npy'
path_local_last_file_name = os.path.join(flying, last_file_name)

names_outfiles_steps = [output_file+'0-all.h5m', output_file+'1-all.h5m']

names_outfiles_variables_steps = [os.path.join(flying, 'variables0.npz'),
                                  os.path.join(flying, 'variables1.npz')]

with open(input_file_0, 'r') as f:
    data_loaded = yaml.safe_load(f)
