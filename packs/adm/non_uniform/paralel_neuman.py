import numpy as np
import scipy.sparse as sp

class Subdomain:

    def __init__(self, volumes, ind_diric, ind_neum, val_diric, val_neum, intern_faces, T_global):
        pass







class masterNeumanNonNested:

    def __init__(self, data_impress, elements_lv0, ml_data, n_levels, T_without, wells):
        self.data_impress = data_impress
        self.elements_lv0 = elements_lv0
        self.ml_data = ml_data
        self.n_levels = n_levels

    def get_subdomains(self):

        '''
            ordem de envio:

                volumes: global id dos volumes locais
                ind_diric: indice dos volumes com pressao prescrita
                ind_neum: indice dos volumes com vazao prescrita
                val_diric: valores de pressao prescrita
                val_neum: valores de vazao prescrita
                local_transm: transmissibilidade local

                all_faces: todas faces do coarse volume
                intern_faces: faces internas do coarse volume
                intersect_faces: faces na interseccao
        '''

        list_of_subdomains = []
        levels = self.data_impress['LEVEL']
        pms = self.data_impress['pms']
        remaped_internal_faces = self.elements_lv0['remaped_internal_faces']
        neig_internal_faces = self.elements_lv0['neig_internal_faces']
        n_volumes = len(levels)


        for level in range(1, self.n_levels):
            str_level = str(level)
            set_level = set([level])

            all_gids_coarse = self.data_impress['GID_'+ str_level]
            # all_local_ids_coarse = self.data_impress['COARSE_LOCAL_ID_'+ str_level]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+ str_level]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+ str_level]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+ str_level]
            all_faces = self.ml_data['coarse_faces_level_'+ str_level]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+ str_level]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+ str_level]
            gids_level = np.unique(all_gids_coarse[levels==level])

            import pdb; pdb.set_trace()

            for gidc in gids_level:

                intersect_faces = all_intersect_faces[coarse_ids==gidc][0] # faces na interseccao
                intern_local_faces = all_intern_faces[coarse_ids==gidc][0] # faces internas
                faces = all_faces[coarse_ids==gidc][0] # faces do volume
                intern_boundary_volumes = all_intern_boundary_volumes[coarse_ids==gidc][0] # volumes internos no contorno
                vertex = all_fine_vertex[coarse_ids==gidc]
                pressure_vertex = pms[vertex]
                volumes = self.elements_lv0['volumes'][all_gids_coarse==gidc]
                level_volumes = levels[volumes]

                if set_level ^ set(level_volumes):
                    print('tratamento aqui')
                    import pdb; pdb.set_trace()
                    pass
                else:

                    adjs_intersect_faces = neig_internal_faces[remaped_internal_faces[intersect_faces]]
                    v0 = adjs_intersect_faces
                    volumes_for_pms_flux = np.unique(np.concatenate([adjs_intersect_faces]))
                    pms0 = pms[adjs_intersect_faces[:,0]]
                    pms1 = pms[adjs_intersect_faces[:,1]]
                    t0 = self.data_impress['transmissibility'][intersect_faces]
                    pms_flux_faces = get_flux_faces(pms1, pms0, t0)

                    lines = np.concatenate([v0[:, 0], v0[:, 1]])
                    cols = np.repeat(0, len(lines))
                    data = np.concatenate([pms_flux_faces, -pms_flux_faces])
                    flux_pms_volumes = sp.csc_matrix((data, (lines, cols)), shape=(n_volumes, 1)).toarray().flatten()
                    presc_flux_intern_boundary_volumes = flux_pms_volumes[intern_boundary_volumes]

                    ind_neum = intern_boundary_volumes
                    val_neum = presc_flux_intern_boundary_volumes
                    ind_diric = vertex
                    val_diric = pms[vertex]















        import pdb; pdb.set_trace()
        return list_of_subdomains





    def run(self):
        list_of_subdomains = self.get_subdomains()


def get_flux_faces(p1, p0, t0, flux_grav_faces=None):

    if flux_grav_faces:
        flux = -((p1 - p0) * t0 - flux_grav_faces)
    else:
        flux = -((p1 - p0) * t0 - flux_grav_faces)

    return flux
