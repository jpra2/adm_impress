from packs.manager.mesh_data import MeshData
import os
from packs import defpaths
import scipy.io  as sio
import numpy as np
import matplotlib.pyplot as plt

mesh_name = '20x1x1.h5m'
mesh_path = os.path.join(defpaths.mesh, mesh_name)

data_path = os.path.join(defpaths.flying, 'resp_p1.mat')
mdata = sio.loadmat(data_path)

folder_to_export = 'p1_mat'
name_export = os.path.join(folder_to_export, 'p1_')

mesh_data = MeshData(dim=3, mesh_path=mesh_path)

tags = ['pressure', 'saturation']
for tag in tags:
    mesh_data.create_tag(tag)

n_pressure = len(mdata['all_pressures'])
n_volumes = len(mdata['all_pressures'][0])
volumes = np.arange(n_volumes)
for i in range(n_pressure):
    pressure = mdata['all_pressures'][i]
    saturation = mdata['all_saturations'][i]
    mesh_data.insert_tag_data('pressure', pressure, 'volumes', volumes)
    mesh_data.insert_tag_data('saturation', saturation, 'volumes', volumes)
    to_export_name = name_export + str(i)
    mesh_data.export_all_elements_type_to_vtk(to_export_name, 'volumes')


vpis = mdata['all_vpi']
oil_cumulative = np.absolute(mdata['cumulative_oil_prod'])
wor_ratio = mdata['all_wor_ratio']
qo_flux = mdata['all_qo_flux']




import pdb; pdb.set_trace()

print(mesh_data)
