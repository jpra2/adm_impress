import numpy as np

class Subdomain:

    def __init__(self, data_impress, elements_lv0, ml_data, level):
        coarse_ids = np.unique(data_impress['GID_'+str(level)])






class masterNeumanNonNested:

    def __init__(self, data_impress, elements_lv0, ml_data, n_levels):
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

        import pdb; pdb.set_trace()

        for level in range(1, self.n_levels):
            str_level = str(level)
            set_level = set([level])

            all_gids_coarse = self.data_impress['GID_'+ str_level]
            # all_local_ids_coarse = self.data_impress['COARSE_LOCAL_ID_'+ str_level]
            all_intern_boundary_volumes = self.ml_data['internal_boundary_fine_volumes_level_'+ str_level]
            all_intersect_faces = self.ml_data['coarse_intersect_faces_level_'+ str_level]
            all_intern_faces = self.ml_data['coarse_internal_faces_level_'+ st_level]
            all_faces = self.ml_data['coarse_faces_level_'+ str_level]
            all_fine_vertex = self.ml_data['fine_vertex_coarse_volumes_level_'+ str_level]
            coarse_ids = self.ml_data['coarse_primal_id_level_'+ str_level]
            gids_level = np.unique(all_gids_coarse[levels==level])

            for cid in gids_level:
                volumes = self.elements_lv0['volumes'][all_gids_coarse==cid]
                level_volumes = levels[volumes]


                if set_level ^ set(level_volumes)  :





    def run(self):
        list_of_subdomains = self.get_subdomains()
