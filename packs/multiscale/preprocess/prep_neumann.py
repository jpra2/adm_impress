import numpy as np

class NeumannSubdomains:
    def __init__(self,elements_lv0,ml_data, data_impress, wells):
        self.wells=wells['all_wells']
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
        self.T=[]
        self.neumann_subds=[]
        self.lines_d=[]
        self.lines_off=[]
        self.cols_d=[]
        self.cols_off=[]
        self.volumes=[]
        self.adj0=[]
        self.adj1=[]
        self.ind_diric_local=[]
        self.ind_neum_local=[]
        self.b=[]
        self.l=[]
        self.c=[]
        self.local_bound_ids=[]
        self.local_vertex=[]

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
            # import pdb; pdb.set_trace()
            pocos=np.intersect1d(self.wells,gid0[all_gids_coarse==gidc])
            if len(pocos)>0:
                vertex=pocos
            else:
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
            self.volumes.append(volumes)
            map_gid_in_lid = np.repeat(-1, volumes.max()+1)
            map_gid_in_lid[volumes] = np.arange(len(volumes))
            self.adj0.append(map_gid_in_lid[adj_intern_local_faces[:,0]])
            self.adj1.append(map_gid_in_lid[adj_intern_local_faces[:,1]])
            self.ind_diric_local.append(map_gid_in_lid[vertex])
            self.ind_neum_local.append(map_gid_in_lid[intern_boundary_volumes])


            self.T.append(np.zeros((len(volumes),len(volumes))))
            self.b.append(np.zeros(len(volumes)))
            map_volumes=np.repeat(-1,adjs.max()+1)
            map_volumes[volumes]=range(len(volumes))

            ind_diric_local=map_volumes[vertex]

            l=map_volumes[adjs[:,0]]
            c=map_volumes[adjs[:,1]]
            if (map_volumes[volumes]!=map_gid_in_lid[volumes]).sum():
                print("ffkffldifere")
                import pdb; pdb.set_trace()
            self.local_vertex.append(map_volumes[all_fine_vertex[coarse_ids==gidc]])
            self.l.append(l)
            self.c.append(c)
            lines=np.concatenate([l,c,l,c])
            cols=np.concatenate([c,l,l,c])
            self.lines.append(lines)
            self.cols.append(cols)
            self.lines_d.append(np.concatenate([l,c]))
            self.cols_d.append(np.concatenate([l,c]))
            self.lines_off.append(np.concatenate([l,c]))
            self.cols_off.append(np.concatenate([c,l]))

            data_impress['initial_diag'][volumes]=map_volumes[volumes]
            data_impress['raz_flux_tag'][volumes[map_volumes[volumes]==ind_diric_local]]=100
            v0=adjs_intersect_faces.T
            map_vol = np.repeat(-1,adjs_intersect_faces.max()+1)
            map_vol[intern_boundary_volumes] = np.arange(len(intern_boundary_volumes))

            ls=np.concatenate([v0[0], v0[1]])
            self.local_bound_ids.append(map_vol[ls])
            data_impress['val_diric'][vertex]=1
            data_impress['val_neum'][intern_boundary_volumes]=1

            self.neumann_subds.append(PrimalSubdomain(elements_lv0, ml_data, data_impress, gidc,all_fine_vertex[coarse_ids==gidc]))

class PrimalSubdomain:
    def __init__(self,elements_lv0, ml_data, data_impress, gidc, local_vertex):
        self.create_primal_subdomain(elements_lv0, ml_data, data_impress, gidc, local_vertex)

    def create_primal_subdomain(self,elements_lv0, ml_data, data_impress, gidc, local_vertex):
        remaped_internal_faces = elements_lv0['remaped_internal_faces']
        neig_internal_faces = elements_lv0['neig_internal_faces']
        # gid0 = data_impress['GID_0']
        all_gids_coarse = data_impress['GID_1']
        all_intern_boundary_volumes = ml_data['internal_boundary_fine_volumes_level_1']
        all_intersect_faces = ml_data['coarse_intersect_faces_level_1']
        all_intern_faces = ml_data['coarse_internal_faces_level_1']
        # all_faces = ml_data['coarse_faces_level_1']
        all_fine_vertex = ml_data['fine_vertex_coarse_volumes_level_1']
        coarse_ids = ml_data['coarse_primal_id_level_1']
        # gids_level = np.unique(all_gids_coarse)

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
        self.local_vertex=local_vertex
        adjs=adj_intern_local_faces
        volumes=np.unique(adjs)
        map_volumes=np.repeat(-1,adjs.max()+1)
        map_volumes[volumes]=range(len(volumes))
        self.map_volumes = map_volumes
        ind_diric_local=map_volumes[vertex]
        l=map_volumes[adjs[:,0]]
        c=map_volumes[adjs[:,1]]
        lines=np.concatenate([l,c,l,c])
        cols=np.concatenate([c,l,l,c])
        self.lines=lines
        self.cols=cols
        self.gid1=gidc
        self.volumes=volumes
        self.ind_diric_local=ind_diric_local
