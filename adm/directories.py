import os
import yaml

entities_lv0 = ['nodes', 'edges', 'faces', 'volumes']

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path_ant = os.getcwd()

impress = 'impress'
tcc = 'tcc'
flying = 'flying'
adm = 'adm'
input_file_0 = 'inputs0.yml'

path_adm = os.path.join(parent_dir, adm)
path_impress = os.path.join(parent_dir, impress)

### usar durante a simulacao tcc path
flying_geral = os.path.join(parent_dir, flying)

with open(input_file_0, 'r') as f:
    data_loaded = yaml.load(f)
