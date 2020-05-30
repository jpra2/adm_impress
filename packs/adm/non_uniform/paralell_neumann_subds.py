import time
import numpy as np
from scipy.sparse import csc_matrix

def get_subdomains_paralell(neumann_subds, transmissibility_faces, pms, comm=None):
    print("starting subdomain preprocess")
    tini=time.time()
    pms_flux_faces = np.zeros(len(transmissibility_faces))
    n_volumes = len(pms)
    list_of_subdomains = []
    for neumann_subd in neumann_subds:
        intersect_faces=neumann_subd.intersect_faces
        intern_local_faces = neumann_subd.intern_local_faces
        intern_boundary_volumes=neumann_subd.intern_boundary_volumes
        ind_diric = neumann_subd.ind_diric
        adjs_intersect_faces = neumann_subd.adjs_intersect_faces
        adj_intern_local_faces = neumann_subd.adj_intern_local_faces

        ind_neum = neumann_subd.intern_boundary_volumes
        lines = neumann_subd.lines
        cols = neumann_subd.cols
        ind_diric_local = neumann_subd.ind_diric_local

        transmissibility = transmissibility_faces[intern_local_faces]
        val_diric = pms[ind_diric]
        v0 = adjs_intersect_faces
        pms0 = pms[v0[:,0]]
        pms1 = pms[v0[:,1]]
        t0 = transmissibility_faces[intersect_faces]
        pms_flux_faces_local = get_flux_faces(pms1, pms0, t0)
        pms_flux_faces[intersect_faces] = pms_flux_faces_local

        map_vol = np.repeat(-1,adjs_intersect_faces.max()+1)
        map_vol[intern_boundary_volumes] = np.arange(len(intern_boundary_volumes))
        ls=np.concatenate([v0[:, 0], v0[:, 1]])
        ds=np.concatenate([pms_flux_faces_local, -pms_flux_faces_local])
        local_ids=map_vol[ls]
        presc_flux_intern_boundary_volumes=np.bincount(local_ids[local_ids>-1], weights=ds[local_ids>-1])
        val_neum = list(presc_flux_intern_boundary_volumes)
        list_of_subdomains.append(Subdomain(ind_diric, ind_neum, val_diric, val_neum, intern_local_faces, adj_intern_local_faces, transmissibility, lines, cols, ind_diric_local))
    print('preprocessing process finished after {} seconds'.format(time.time()-tini))
    comm.send([list_of_subdomains, pms_flux_faces])


def get_flux_faces(p1, p0, t0, flux_grav_faces=None):

    if flux_grav_faces != None:
        flux = -((p1 - p0) * t0 - flux_grav_faces)
    else:
        flux = -((p1 - p0) * t0)

    return flux

class Subdomain():

    def __init__(self, ind_diric, ind_neum, val_diric, val_neum,
        intern_faces, adjs_intern_faces, transmissibility, lines, cols, ind_diric_local):
        self.transmissibilities=transmissibility
        self.T_local = self.get_local_matrix_with_boundary_condition(lines, cols, transmissibility, ind_diric_local)
        volumes=np.unique(adjs_intern_faces)
        self.volumes = volumes
        self.ids_local = np.arange(len(volumes))
        self.ind_diric = ind_diric
        self.ind_neum = ind_neum
        self.val_diric = val_diric
        self.val_neum = val_neum
        self.intern_faces = intern_faces
        self.adjs_intern_faces = adjs_intern_faces
        self.map_gid_in_lid = np.repeat(-1, volumes.max()+1)
        self.map_gid_in_lid[volumes] = self.ids_local

    def get_local_matrix_with_boundary_condition(self,lines, cols, transmissibility, ind_diric_local):
        d=transmissibility.copy()
        data=np.concatenate([d,d,-d,-d])
        data[lines==ind_diric_local]=0
        data[(lines==ind_diric_local) & (cols==ind_diric_local)]=0.5
        n_vols=lines.max()+1
        T2=csc_matrix((data,(lines, cols)),shape=(n_vols,n_vols))
        return T2
