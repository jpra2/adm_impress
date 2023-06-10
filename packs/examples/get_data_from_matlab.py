from packs.manager.mesh_data import MeshData
import os
from packs import defpaths
import scipy.io  as sio
import numpy as np
import matplotlib.pyplot as plt

#######################################################
mesh_name = '20x1x1.h5m'
mesh_path = os.path.join(defpaths.mesh, mesh_name)

data_path = os.path.join(defpaths.flying, 'resp_p1.mat')
mdata = sio.loadmat(data_path)

folder_to_export = 'p1_mat'
problem_name = '_p1_'
vpi_name = 'vpi_'
name_export = os.path.join(folder_to_export, 'p1_')

import pdb; pdb.set_trace()
#######################################################


#######################################################
mesh_name = '21x21x1.h5m'
mesh_path = os.path.join(defpaths.mesh, mesh_name)

data_path = os.path.join(defpaths.flying, 'resp_p6.mat')
mdata = sio.loadmat(data_path)


folder_to_export = 'p6_mat'
problem_name = '_p6_'
vpi_name = 'vpi_'
name_export = os.path.join(folder_to_export, 'p6_')
#######################################################

mesh_data = MeshData(dim=3, mesh_path=mesh_path)

def test_value_sat():
    mesh_data.create_tag('sat_test')
    sat_path = os.path.join(defpaths.flying, 'sat.mat')
    data = sio.loadmat(sat_path)
    sat = data['S'].flatten()
    n = len(sat)
    import pdb; pdb.set_trace()
    vols = np.arange(n)
    mesh_data.insert_tag_data('sat_test', sat, 'volumes', vols)
    to_export_name = os.path.join(defpaths.flying, 'satu_test')
    mesh_data.export_all_elements_type_to_vtk(to_export_name, 'volumes')

# test_value_sat()
# import pdb; pdb.set_trace()

tags = ['pressure', 'saturation']
for tag in tags:
    mesh_data.create_tag(tag)

n_pressure = len(mdata['all_pressures'])
n_volumes = len(mdata['all_pressures'][0])
volumes = np.arange(n_volumes)

permeabilityx = mdata['perm'][:,1,1].astype(float)
mesh_data.create_tag('permeabilityx', data_size=1)
mesh_data.insert_tag_data('permeabilityx', permeabilityx, 'volumes', volumes)

all_vpis = np.around(mdata['all_vpi'].flatten(), decimals=3)
for i in range(n_pressure):
    pressure = mdata['all_pressures'][i]
    saturation = mdata['all_saturations'][i]
    mesh_data.insert_tag_data('pressure', pressure, 'volumes', volumes)
    mesh_data.insert_tag_data('saturation', saturation, 'volumes', volumes)
    name = problem_name + str(i)
    to_export_name = os.path.join(folder_to_export, name)
    mesh_data.export_all_elements_type_to_vtk(to_export_name, 'volumes')


# vpis = mdata['all_vpi']
# oil_cumulative = np.absolute(mdata['cumulative_oil_prod'])
# wor_ratio = mdata['all_wor_ratio']
# qo_flux = mdata['all_qo_flux']

import pdb; pdb.set_trace()

print(mesh_data)