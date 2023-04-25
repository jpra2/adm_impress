from packs.manager.mesh_data import MeshData
import os
from packs import defpaths
import scipy.io  as sio
import numpy as np

mesh_name = '20x1x1.h5m'
mesh_path = os.path.join(defpaths.mesh, mesh_name)

data_path = os.path.join(defpaths.flying, 'resp_p1.mat')
mdata = sio.loadmat(data_path)

name_export = 'p1_'

mesh_data = MeshData(dim=3, mesh_path=mesh_path)

tags = ['pressure', 'saturation']
for tag in tags:
    mesh_data.create_tag(tag)

n_pressure = len(mdata['all_pressures'])
n_volumes = len(mdata['all_pressures'][0])
for i in range(n_pressure):
    pressure = mdata['all_pressures'][i]
    saturation = mdata['all_saturations'][i]
    mesh_data.insert_tag_data('pressure', pressure, 'volumes', np.arange(n_volumes))
    mesh_data.insert_tag_data('saturation', saturation, 'volumes', np.arange(n_volumes))
    to_export_name = name_export + str(i)
    mesh_data.export_all_elements_type_to_vtk(to_export_name, 'volumes')


import pdb; pdb.set_trace()

print(mesh_data)
