import numpy as np
import scipy.sparse as sp

class masterNeumanNonNested:
    def __init__(self,M, data_impress, elements_lv0, ml_data, n_levels, wells,neumann_subds):
        self.data_impress = data_impress
        self.elements_lv0 = elements_lv0
        self.ml_data = ml_data
        self.n_levels = n_levels
        self.wells = wells
        # self.mesh = M
        self.neumann_subds = neumann_subds

    def run(self):
        list_of_subdomains, self.global_ms_flux_faces = self.get_subdomains_vec()
        faces_intern, ms_flux_faces, vols, pcorr = LocalSolution(list_of_subdomains).run()

        self.global_ms_flux_faces[np.concatenate(faces_intern)] = np.concatenate(ms_flux_faces)
        global_pcorr = np.zeros(len(self.data_impress['GID_0']))
        global_pcorr[np.concatenate(vols)] = np.concatenate(pcorr)
        return self.global_ms_flux_faces.copy(), global_pcorr

    def get_subdomains_vec(self):
        list_of_subdomains = []
        pms_flux_faces = np.zeros(len(self.elements_lv0['faces']))
        pms = self.data_impress['pms']
        transmissibilities=self.data_impress['transmissibility']
        coarse_ids = self.ml_data['coarse_primal_id_level_1']
        gids_level = np.unique(coarse_ids)
        for gidc in gids_level:
            intersect_faces=self.neumann_subds.intersect_faces[gidc]
            intern_local_faces = self.neumann_subds.intern_local_faces[gidc]
            intern_boundary_volumes=self.neumann_subds.intern_boundary_volumes[gidc]
            ind_diric = self.neumann_subds.ind_diric[gidc]
            adjs_intersect_faces = self.neumann_subds.adjs_intersect_faces[gidc]
            transmissibility = transmissibilities[intern_local_faces]
            ind_neum = self.neumann_subds.intern_boundary_volumes[gidc]
            ind_diric_local=self.neumann_subds.ind_diric_local[gidc]
            T=self.neumann_subds.T[gidc]
            b=self.neumann_subds.b[gidc]
            volumes=self.neumann_subds.volumes[gidc]
            adj0=self.neumann_subds.adj0[gidc]
            adj1=self.neumann_subds.adj1[gidc]
            ind_diric_local=self.neumann_subds.ind_diric_local[gidc]
            ind_neum_local=self.neumann_subds.ind_neum_local[gidc]
            l=self.neumann_subds.l[gidc]
            c=self.neumann_subds.c[gidc]
            local_bound_ids=self.neumann_subds.local_bound_ids[gidc]

            val_diric = pms[ind_diric]
            v0 = adjs_intersect_faces.T

            pms0 = pms[v0[0]]
            pms1 = pms[v0[1]]
            t0 = self.data_impress['transmissibility'][intersect_faces]
            pms_flux_faces_local = get_flux_faces(pms1, pms0, t0)
            pms_flux_faces[intersect_faces] = pms_flux_faces_local

            ds=np.concatenate([pms_flux_faces_local, -pms_flux_faces_local])

            val_neum=np.bincount(local_bound_ids[local_bound_ids>-1], weights=ds[local_bound_ids>-1])

            list_of_subdomains.append(Subdomain(val_diric, val_neum, intern_local_faces, transmissibility, ind_diric_local, T, volumes, adj0, adj1, ind_neum_local, b, l, c))

        return list_of_subdomains, pms_flux_faces

def get_flux_faces(p1, p0, t0, flux_grav_faces=None):

    if flux_grav_faces != None:
        flux = -((p1 - p0) * t0 - flux_grav_faces)
    else:
        flux = -((p1 - p0) * t0)

    return flux

class Subdomain():
    def __init__(self, val_diric, val_neum, intern_faces, transmissibility, ind_diric_local, T, volumes, adj0, adj1, ind_neum_local, b, l, c):
        self.transmissibilities=transmissibility
        self.T_local = self.get_local_matrix_with_boundary_condition(transmissibility, ind_diric_local, T, l, c)
        self.volumes = volumes
        self.intern_faces = intern_faces
        self.adj0=adj0
        self.adj1=adj1
        b[ind_diric_local] = val_diric
        b[ind_neum_local] = val_neum
        self.b=b

    def get_local_matrix_with_boundary_condition(self, transmissibility, ind_diric_local, T, l, c):
        d=transmissibility
        T[l,c]=d
        T[c,l]=d
        np.fill_diagonal(T,np.bincount(np.concatenate([l,c]),weights=-np.concatenate([d,d])))
        T[ind_diric_local]=0
        T[ind_diric_local,ind_diric_local]=1
        return T

class LocalSolution:
    def __init__(self, subdomains):
        self.subdomains = subdomains

    def run(self):
        faces_intern=[]
        ms_flux_faces=[]
        vols=[]
        pcorr=[]
        for subd in self.subdomains:
            volumes = subd.volumes
            T_local = subd.T_local
            intern_faces = subd.intern_faces
            transmissibilities = subd.transmissibilities
            adj0=subd.adj0
            adj1=subd.adj1
            b=subd.b
            x=np.linalg.solve(T_local,b)
            p0 = x[adj0]
            p1 = x[adj1]
            ms_flux = -transmissibilities*(p1 - p0) #get_flux_faces(p1, p0, t0)
            faces_intern.append(intern_faces)
            ms_flux_faces.append(ms_flux)
            vols.append(volumes)
            pcorr.append(x)
        faces_intern=np.concatenate([faces_intern])
        ms_flux_faces=np.concatenate([ms_flux_faces])
        vols=np.concatenate([vols])
        pcorr=np.concatenate([pcorr])
        return faces_intern, ms_flux_faces, vols, pcorr
