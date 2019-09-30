import os

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

impress = 'impress'
tcc = 'tcc'
flying = 'flying'
adm = 'adm'
input_file_1 = 'inputs.yaml'

path_adm = os.path.join(parent_dir, adm)
path_impress = os.path.join(parent_dir, impress)

### usar durante a simulacao tcc path
flying_geral = os.path.join(parent_dir, flying)
