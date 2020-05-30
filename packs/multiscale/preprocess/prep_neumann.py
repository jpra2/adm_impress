import numpy as np

class NeumannSubdomains:
    def __init__(self,elements_lv0,ml_data, data_impress):
        self.ind_diric=[]
        self.ind_neum=[]
        self.intern_local_faces=[]
        self.adj_intern_local_faces=[]
        self.adjs_intersect_faces=[]
        self.intern_boundary_volumes=[]
        self.intersect_faces=[]
        self.lines=[]
        self.cols=[]
        self.ind_diric_local=[]
        self.neumann_subds=[]
        self.create_neumann_subdomains(elements_lv0, ml_data, data_impress)

    def create_neumann_subdomains(self,elements_lv0, ml_data, data_impress):
        remaped_internal_faces = elements_lv0['remaped_internal_faces']
        neig_internal_faces = elements_lv0['neig_internal_faces']
        gid0 = data_impress['GID_0']
        all_gids_coarse = data_impress['GID_1']
        all_intern_boundary_volumes = ml_data['internal_boundary_fine_volumes_level_1']
        all_intersect_faces = ml_data['coarse_intersect_faces_level_1']
        all_intern_faces = ml_data['coarse_internal_faces_level_1']
        all_faces = ml_data['coarse_faces_level_1']
        all_fine_vertex = ml_data['fine_vertex_coarse_volumes_level_1']
        coarse_ids = ml_data['coarse_primal_id_level_1']
        gids_level = np.unique(all_gids_coarse)
        for gidc in gids_level:
            intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
            adjs_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
            intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
            adj_intern_local_faces = neig_internal_faces[remaped_internal_faces[intern_local_faces]]
            intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
            vertex = all_fine_vertex[coarse_ids==gidc]

            self.intersect_faces.append(intersect_faces)
            self.ind_diric.append(vertex)
            self.ind_neum.append(intern_boundary_volumes)
            self.intern_local_faces.append(intern_local_faces)
            self.adj_intern_local_faces.append(adj_intern_local_faces)
            self.adjs_intersect_faces.append(adjs_intersect_faces)
            self.intern_boundary_volumes.append(intern_boundary_volumes)

            adjs=adj_intern_local_faces
            volumes=np.unique(adjs)
            map_volumes=np.repeat(-1,adjs.max()+1)
            map_volumes[volumes]=range(len(volumes))
            ind_diric_local=map_volumes[vertex]
            l=map_volumes[adjs[:,0]]
            c=map_volumes[adjs[:,1]]
            lines=np.concatenate([l,c,l,c])
            cols=np.concatenate([c,l,l,c])
            self.lines.append(lines)
            self.cols.append(cols)
            self.ind_diric_local.append(ind_diric_local)
            self.neumann_subds.append(PrimalSubdomain(elements_lv0, ml_data, data_impress, gidc))

class PrimalSubdomain:
    def __init__(self,elements_lv0, ml_data, data_impress, gidc):
        self.create_primal_subdomain(elements_lv0, ml_data, data_impress, gidc)

    def create_primal_subdomain(self,elements_lv0, ml_data, data_impress, gidc):
        remaped_internal_faces = elements_lv0['remaped_internal_faces']
        neig_internal_faces = elements_lv0['neig_internal_faces']
        gid0 = data_impress['GID_0']
        all_gids_coarse = data_impress['GID_1']
        all_intern_boundary_volumes = ml_data['internal_boundary_fine_volumes_level_1']
        all_intersect_faces = ml_data['coarse_intersect_faces_level_1']
        all_intern_faces = ml_data['coarse_internal_faces_level_1']
        all_faces = ml_data['coarse_faces_level_1']
        all_fine_vertex = ml_data['fine_vertex_coarse_volumes_level_1']
        coarse_ids = ml_data['coarse_primal_id_level_1']
        gids_level = np.unique(all_gids_coarse)

        intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
        adjs_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
        intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
        adj_intern_local_faces = neig_internal_faces[remaped_internal_faces[intern_local_faces]]
        intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
        vertex = all_fine_vertex[coarse_ids==gidc]

        self.intersect_faces=intersect_faces
        self.ind_diric=vertex
        self.ind_neum=intern_boundary_volumes
        self.intern_local_faces=intern_local_faces
        self.adj_intern_local_faces=adj_intern_local_faces
        self.adjs_intersect_faces=adjs_intersect_faces
        self.intern_boundary_volumes=intern_boundary_volumes

        adjs=adj_intern_local_faces
        volumes=np.unique(adjs)
        map_volumes=np.repeat(-1,adjs.max()+1)
        map_volumes[volumes]=range(len(volumes))
        ind_diric_local=map_volumes[vertex]
        l=map_volumes[adjs[:,0]]
        c=map_volumes[adjs[:,1]]
        lines=np.concatenate([l,c,l,c])
        cols=np.concatenate([c,l,l,c])
        self.lines=lines
        self.cols=cols
        self.ind_diric_local=ind_diric_local
