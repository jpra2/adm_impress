
from .. import directories as direc


import numpy as np

class Preprocess0:

    def set_permeability(self, M):
        data_loaded = direc.data_loaded

        read = data_loaded['read_permeability']
        if read:
            set_permeability_spe10(M)
        else:
            self.set_permeability_regions(M)

    def set_permeability_regions(self, M):
        pass


def set_permeability_spe10(M):
    ks = np.load(direc.data_loaded['file_name_permeability'])['perms']
    phi = np.load(direc.data_loaded['file_name_permeability'])['phi']
    phi = phi.flatten()

    nx = 60
    ny = 220
    nz = 85
    perms = []
    phis = []

    k = 1.0  #para converter a unidade de permeabilidade
    centroids=M.data.centroids[direc.entities_lv0[3]]
    cont=0
    for v in self.all_volumes:
        permeabilidade = np.zeros(9)
        centroid = centroids[cont]
        cont+=1
        ijk = np.array([centroid[0]//20.0, centroid[1]//10.0, centroid[2]//2.0])
        e = int(ijk[0] + ijk[1]*nx + ijk[2]*nx*ny)
        # perm = ks[e]*k
        # fi = phi[e]
        permeabilidade[0] = ks[e][0]
        permeabilidade[4] = ks[e][1]
        permeabilidade[8] = ks[e][2]
        permeabilidade *= k
        perms.append(permeabilidade)
        phis.append(phi[e])
        self.mb.tag_set_data(self.perm_tag, v, permeabilidade)

    self.mb.tag_set_data(self.phi_tag, self.all_volumes, phis)
