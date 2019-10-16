from pymoab import core, types, rng, topo_util
import numpy as np

def get_box_dep0(all_centroids, limites):

    '''
    all_centroids->coordenadas dos centroides do conjunto
    limites-> diagonal que define os volumes objetivo (numpy array com duas coordenadas)
    Retorna os indices cujo centroide está dentro de limites
    '''

    inds0 = np.where(all_centroids[:,0] > limites[0,0])[0]
    inds1 = np.where(all_centroids[:,1] > limites[0,1])[0]
    inds2 = np.where(all_centroids[:,2] > limites[0,2])[0]
    c1 = set(inds0) & set(inds1) & set(inds2)
    inds0 = np.where(all_centroids[:,0] < limites[1,0])[0]
    inds1 = np.where(all_centroids[:,1] < limites[1,1])[0]
    inds2 = np.where(all_centroids[:,2] < limites[1,2])[0]
    c2 = set(inds0) & set(inds1) & set(inds2)
    inds_vols = list(c1 & c2)
    return inds_vols

def get_box(all_centroids, limites):
    ids = np.arange(len(all_centroids))

    '''
    all_centroids->coordenadas dos centroides do conjunto
    limites-> diagonal que define os volumes objetivo (numpy array com duas coordenadas)
    Retorna os indices cujo centroide está dentro de limites
    '''

    inds_vols = ids[
            (all_centroids[:,0] > limites[0,0]) & (all_centroids[:,1] > limites[0,1]) & \
            (all_centroids[:,2] > limites[0,2]) & (all_centroids[:,0] < limites[1,0]) & \
            (all_centroids[:,1] < limites[1,1]) & (all_centroids[:,2] < limites[1,2])]
    return inds_vols

def getting_tag(mb, name, n, t1, t2, create, entitie, tipo, tags, tags_to_infos):
    types_data = ['handle', 'integer', 'array', 'double']
    entities = ['nodes', 'edges', 'faces', 'volumes', 'root_set', 'intern_faces', 'boundary_faces', 'vols_viz_face',
                'coarse_volumes_lv1', 'coarse_volumes_lv2']

    assert tipo in types_data, f'tipo nao listado: {tipo}'

    if entitie not in entities:
        raise NameError(f'\nA entidade {entitie} nao esta na lista\n')

    tag = mb.tag_get_handle(name, n, t1, t2, create)
    tag_to_infos = dict(zip(['entitie', 'type', 'n'], [entitie, tipo, n]))
    tags[name] = tag
    tags_to_infos[name] = tag_to_infos

def Min_Max(e, M1):
    verts = M1.mb.get_connectivity(int(e))
    coords = M1.mb.get_coords(verts).reshape(len(verts), 3)
    xmax, ymax, zmax = coords.max(axis=0)
    xmin, ymin, zmin = coords.min(axis=0)
    return([xmin,xmax,ymin,ymax,zmin,zmax])
