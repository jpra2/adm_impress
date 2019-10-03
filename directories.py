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
                     'area': 'area', 'dist_cent': 'dist_cent'}

impress = 'impress'
tcc = 'tcc'
flying = 'flying'
adm = 'adm'
input_file_0 = 'inputs0.yml'

parent_dir = os.path.dirname(os.path.abspath(__file__)) # adm_impress_dir
adm_impress_dir = parent_dir
path_ant = os.getcwd()

path_adm = os.path.join(adm_impress_dir, adm)
path_impress = direc_impress.impress_path

### usar durante a simulacao tcc path
path_flying_geral = os.path.join(adm_impress_dir, flying)

with open(input_file_0, 'r') as f:
    data_loaded = yaml.safe_load(f)
