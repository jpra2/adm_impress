import os
import yaml
# from impress.preprocessor import directories as direc_impress
from .preprocess import directories_impress as direc_impress

entities_lv0 = ['nodes', 'edges', 'faces', 'volumes']
# entities_lv0 = direc_impress.entities_lv0
entities_lv0_0 = direc_impress.entities_lv0_0

names_data_loaded_lv0 = ['read_permeability', 'file_name_permeability', 'Permeability', 'Crs',
                         'Saturation', 'set_permeability', 'set_porosity', 'Porosity']

names_data_loaded_lv2 = ['type', 'value', 'p0', 'p1']

types_region_data_loaded = ['all', 'box']
types_region_for_saturation = ['all', 'wells', 'box']
types_presc = ['dirichlet', 'neumann']
types_wells = ['injector', 'producer']
types_values = ['molar', 'volumetric']

types_simulation = ['compositional', 'biphasic', 'monophasic']

# variables_impress = {'permeability': 'permeability', 'poro': 'poro', 'k_harm': 'k_harm',
#                      'area': 'area', 'dist_cent': 'dist_cent', 'u_normal': 'u_normal'}

variables_impress = dict()

impress = 'impress'
tcc = 'tcc'
flying = 'flying'
adm = 'adm'
# input_file_0 = 'input_cards/inputs0.yml'
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

#with open('input_cards/inputs_compositional.yml', 'r') as f:
with open ('input_cards/input_file_name.yml','r') as f:
    names_files_load = yaml.safe_load(f)
    name_input_file_load = names_files_load['name_file']
    name_variable_inputs_file_load = names_files_load['variable_file']
    simulation_type = names_files_load['simulation_type']

with open(name_input_file_load, 'r') as f:
    data_loaded = yaml.safe_load(f)

with open(name_variable_inputs_file_load, 'r') as f:
    variables_loaded = yaml.safe_load(f)

name_load = os.path.join(flying, 'load.npy')
name_hist = os.path.join(flying, 'current_compositional_results.npy')
