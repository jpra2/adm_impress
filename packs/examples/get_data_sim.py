from packs.manager.mesh_data import MeshData
import os
from packs import defpaths
import scipy.io  as sio
import numpy as np
import matplotlib.pyplot as plt
from packs.manager.meshmanager import create_initial_3D_mesh_prperties, MeshProperty, load_mesh_properties
from packs.examples.load_3d_mesh import first_func, second_func

#######################################################
mesh_name = '9x9x1_sim.h5m'
mesh_path = os.path.join(defpaths.mesh, mesh_name)
mesh_properties_name = '9x9x1_sim'

first_func(mesh_path, mesh_properties_name)
second_func(mesh_properties_name)

mesh_properties = load_mesh_properties(mesh_properties_name)

to_export_name = 'test1_sim'

mesh_data = MeshData(dim=3, mesh_path=mesh_path)

nblocks = len(mesh_properties.volumes)
ids = np.arange(nblocks)

mesh_data.create_tag('id', data_type='int')
mesh_data.insert_tag_data('id', ids, 'volumes')

permeability = np.ones(nblocks)*100
mesh_data.create_tag('permeability')
mesh_data.insert_tag_data('permeability', permeability, 'volumes')

porosity = np.ones(nblocks)*0.28
porosity[mesh_properties.volumes_centroids[:, 0] > 720] = 0.3
mesh_data.create_tag('porosity')
mesh_data.insert_tag_data('porosity', porosity, 'volumes')


mesh_data.export_all_elements_type_to_vtk(to_export_name, 'volumes')



