import os

impress_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'impress')

entities_lv0 = ['nodes', 'edges', 'faces', 'volumes']
entities_lv0_0 = ['internal_faces', 'vols_viz_faces', 'vols_viz_internal_faces',
                  'vols_viz_boundary_faces', 'boundary_faces']
names_datas = ['data_size', 'data_format', 'entity', 'level']
data_formats = ['float', 'int', 'bool']
name_variables = 'variables.npz'
name_info_data = 'info_data.txt'

flying = 'flying'
path_local_variables = flying + '/' + name_variables
path_local_info_data = flying + '/' + name_info_data
