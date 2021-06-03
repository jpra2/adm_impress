import os
import yaml
# from impress.preprocessor import directories as direc_impress
from .preprocess import directories_impress as direc_impress

entities_lv0 = ['nodes', 'edges', 'faces', 'volumes']
# entities_lv0 = direc_impress.entities_lv0
entities_lv0_0 = direc_impress.entities_lv0_0

names_data_loaded_lv0 = ['read_permeability', 'file_name_permeability_and porosity', 'Permeability', 'Crs',
                         'Saturation', 'set_permeability', 'set_porosity', 'Porosity']

names_data_loaded_lv2 = ['type', 'value', 'p0', 'p1']

types_region_data_loaded = ['all', 'box', 'ring']
types_region_for_saturation = ['all', 'wells']
types_presc = ['dirichlet', 'neumann']
types_wells = ['injector', 'producer']

# variables_impress = {'permeability': 'permeability', 'poro': 'poro', 'k_harm': 'k_harm',
#                      'area': 'area', 'dist_cent': 'dist_cent', 'u_normal': 'u_normal'}

variables_impress = {'poro': 'phi'}

impress = 'impress'
tcc = 'tcc'
flying = 'flying'
adm = 'adm'
# input_file_0 = 'inputs0.yml'
input_file_finescale_preprocess = 'input_cards/FIM/finescale_preprocess.yml'
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

names_outfiles_steps = [output_file+'0-all.h5m',
                        output_file+'1-all.h5m',
                        output_file+'2-all.h5m',
                        output_file+'3-all.h5m',
                        output_file+'4-all.h5m']

names_outfiles_variables_steps = [os.path.join(flying, 'variables0.npz'),
                                  os.path.join(flying, 'variables1.npz'),
                                  os.path.join(flying, 'variables2.npz'),
                                  os.path.join(flying, 'variables3.npz')
                                  ]
names_datas_contour = os.path.join(flying, 'datas_contour.npz')
file_adm_mesh_def = 'initial_mesh_def.yml'
file_adm_mesh_def = os.path.join('input_cards', file_adm_mesh_def)

with open(input_file_finescale_preprocess, 'r') as f:
    data_loaded = yaml.safe_load(f)

with open(file_adm_mesh_def, 'r') as f:
    file_adm_mesh_def = yaml.safe_load(f)


with open('input_cards/variable_input.yml', 'r') as f:
    variables_loaded = yaml.safe_load(f)

name_load = os.path.join(flying, 'load.npy')
name_hist = os.path.join(flying, 'current_biphasic_results.npy')
