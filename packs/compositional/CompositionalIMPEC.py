from ..directories import data_loaded
from ..data_class.data_manage import dataManager
from ..utils import relative_permeability
from .. import directories as direc
import numpy as np

class CompositionalTPFA(DataManager):
    def __init__(self,data_name: str='CompositionalTPFA.npz'):
        super().__init__(data_impress, load=load)
        # self.M = M
        # self.elements_lv0 = elements_lv0
        # self.relative_permeability = getattr(relative_permeability, self.compositional_data['relative_permeability'])
        # load = data_loaded['load_compositional_data']
        # self.mesh_name = os.path.join(direc.flying, 'compositional_')
        # self.all_compositional_results = self.get_empty_current_compositional_results()
        # self.solver = SolverSp()

    def update_transmissibility(self,M,data_impress,elements_lv0,fluid_properties):
        v0 = elements_lv0['neig_internal_faces']
        internal_faces = elements_lv0['internal_faces']
        pretransmissibility_faces = M.data[M.data.variables_impress['pretransmissibility']]
        pretransmissibility_internal_faces = transmissibility_faces[internal_faces]
        mesh_grid_blocks = len(elements_lv0['volumes']) #número de blocos da malha -  ver como puxar isso

        #z - molar fraction of the component in each phase - ex: z = [[zc1.liq, zc2.liq],[zc1.vap,zc2.vap]
        molar_fractions = np.zeros([fluid_properties.Nc,2,mesh_grid_blocks])
        mass_densities = np.zeros([1,2,mesh_grid_blocks])

        molar_fractions[:,0,:] = fluid_properties.x; molar_fractions[:,1,:] = fluid_properties.y
        mass_densities[0,0,:] = [fluid_properties.rho_L]; mass_densities[1,1,:] = [fluid_properties.rho_V]

        #mobilities = relative_permeabilities/phase_viscosities #2 column vectors - ver o formato das outras coisas pra alterar o dela e vetorizar

        t0 = (molar_fractions*mass_densities*mobilities).sum(axis=1).sum(axis=1)
        t0 = t0*pretransmissibility_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1], v0[:, 0], v0[:, 1]]).flatten()
        cols = np.array([v0[:, 1], v0[:, 0], v0[:, 0], v0[:, 1]]).flatten()
        data = np.array([t0, t0, -t0, -t0]).flatten()

        T = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, self.n_volumes))

        # Falta incluir o termo da derivada de V por P que soma na diagonal principal
        #e o termo da derivada de V por Nk que multiplica a transmissibilidade t0 (ele
        #entra no segundo somatório)
