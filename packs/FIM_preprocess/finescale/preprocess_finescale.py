from packs.running.initial_mesh_properties import initial_mesh
# from packs.FIM_preprocess.finescale.finescale_mesh import MeshManager
from pymoab import core, types, rng, topo_util

import numpy as np
from ... import directories_FIM as direc
def get_finescale_data():
    M=construct_finescale_object(direc.data_loaded['mesh_name'])

    centroids=np.array([[0.5,0.5,1.0],[1.5, 1.5, 1.5]])
    ks, phis = get_permeability_and_phi_spe10(centroids)

@profile
def construct_finescale_object(mesh_file):
    # mb=core.Core()
    # mb.load_file(mesh_file)
    M, elements_lv0, data_impress, wells = initial_mesh()
    Ts=M.k_harm[M.faces.internal].T[0]

    adjs=M.faces.bridge_adjacencies(M.faces.internal_elements[:],2,3)
    import pdb; pdb.set_trace()
    # volumes=mb.get_entities_by_dimension(0,3)
    # nodes=mb.get_entities_by_dimension(0, 0)
    # mtu = topo_util.MeshTopoUtil(mb)
    # mtu.construct_aentities(nodes)
    # faces=mb.get_entities_by_dimension(0, 2)
    # GID_0_tag = mb.tag_get_handle("GID_0", 1, types.MB_TYPE_INTEGER, types.MB_TAG_DENSE, True)
    # centroids=np.array([mtu.get_average_position([v]) for v in volumes])
    # mb.tag_set_data(GID_0_tag, np.array(volumes), range(len(volumes)))
    # adjs=[mtu.get_bridge_adjacencies(face,3, 3) for face in faces]
    # adjs=np.array([ad for ad in adjs if len(ad)==2])
    # areas=get_area(mb, faces)
    # ks, phis = get_permeabilities_and_porosities(centroids)
    return adjs

def get_area(mb, faces):
    nodes=[mb.get_adjacencies(face, 0) for face in faces if len(mb.get_adjacencies(face,3, 3))==2]
    return nodes

def get_permeabilities_and_porosities(centroids):
    load_permeability=direc.data_loaded['read_permeability']
    load_porosity=direc.data_loaded['read_porosity']
    if load_permeability and load_porosity:
        ks, phis = get_permeability_and_phi_spe10(centroids)
    return ks, phis

def get_permeability_and_phi_spe10(centroids):
    data_spe10 = np.load(direc.data_loaded['file_name_permeability_and_porosity'])
    ks = data_spe10['perms']
    phi = data_spe10['phi']
    phi = phi.flatten()

    nx = 60
    ny = 220
    nz = 85

    ijk0 = np.array([centroids[:, 0]//1.0, centroids[:, 1]//1.0, centroids[:, 2]//2.0])
    ee = ijk0[0] + ijk0[1]*nx + ijk0[2]*nx*ny
    ee = ee.astype(np.int32)
    return ks[ee], phi[ee]
