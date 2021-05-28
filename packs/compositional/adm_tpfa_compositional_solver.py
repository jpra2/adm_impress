from packs.compositional.pressure_solver import TPFASolver
from packs.utils.test_functions import test_kwargs_keys, test_instance
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from ..utils import constants as ctes
from packs.multiscale.ms_utils.multiscale_solver import multilevel_pressure_solver
import scipy.sparse as sp
import numpy as np



class AdmTpfaCompositionalSolver(TPFASolver):
    
    _kwargs_keys = {
        'get_pressure': [
            'multilevel_data',
            'multilevel_operators'
        ]
    }
    
    def get_pressure(self, M, wells, fprop, delta_t, **kwargs):
        params = kwargs.get('params')

        # test_kwargs_keys(
        #     AdmTpfaCompositionalSolver._kwargs_keys['get_pressure'],
        #     kwargs.keys()
        # )
        
        T, T_noCC = self.update_transmissibility(M, wells, fprop, delta_t)
        D = self.update_independent_terms(M, fprop, wells, delta_t)
        mlo = kwargs.get('multilevel_operators')
        test_instance(mlo, MultilevelOperators)
        # mlo.run(T_noCC, np.zeros(len(D)), np.zeros(len(D)))
        mlo.run(T_noCC, D, np.zeros(len(D)), return_correction_matrix=False)
        n_levels = 2
        #####
        ##Test with AMS operators
        #####
        prolongation_list = []
        restriction_list = []
        for level in range(1, n_levels):
            prolongation_list.append(mlo[mlo.prolongation + str(level)])
            restriction_list.append(mlo[mlo.restriction + str(level)])

        solution = multilevel_pressure_solver(
            T,
            D,
            prolongation_list,
            restriction_list
        )

        self.P = self.update_pressure(T, D, fprop) # OP*padm
        error = np.absolute(self.P - solution) / self.P
        data_impress = M.data
        data_impress['pressure'][:] = self.P
        data_impress['ms_pressure'][:] = solution
        data_impress['pressure_error'][:] = error
        data_impress.update_variables_to_mesh()
        m1 = M.core.mb.create_meshset()
        M.core.mb.add_entities(m1, M.core.all_volumes)
        M.core.mb.write_file('results/test_comp_1.vtk', [m1])
        import pdb; pdb.set_trace()

        Ft_internal_faces = self.update_total_flux_internal_faces(M, fprop) # pressao local
        self.update_flux_wells(fprop, wells, delta_t)
        # params['dVtdP'] = AdmTpfaCompositionalSolver.dVtP
        # params['dVtdNk'] = AdmTpfaCompositionalSolver.dVtk
        return self.P, Ft_internal_faces, self.q
    
    def update_transmissibility(self, M, wells, fprop, delta_t):
        self.t0_internal_faces_prod = fprop.xkj_internal_faces * \
                                      fprop.Csi_j_internal_faces * \
                                      fprop.mobilities_internal_faces

        ''' Transmissibility '''
        t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = t0 * ctes.pretransmissibility_internal_faces
        # T = np.zeros([ctes.n_volumes, ctes.n_volumes])
        T = sp.lil_matrix((ctes.n_volumes, ctes.n_volumes))

        # Look for a way of doing this not using a loop!!!
        for i in range(ctes.n_components):
            lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
            data = np.array([-t0[i,:], -t0[i,:], +t0[i,:], +t0[i,:]]).flatten()

            # Ta = (sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))).toarray()
            Ta = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))
            # T += Ta * self.dVtk[i,:, np.newaxis]
            T += Ta.multiply(self.dVtk[i,:, np.newaxis])

        # T = T * delta_t
        T *= delta_t
        ''' Transmissibility diagonal term '''
        # diag = np.diag((ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        # T += diag

        # diagT = T.diagonal().flatten()
        # diagT += ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP
        # T.setdiag(diagT)

        T.setdiag(T.diagonal().flatten() + (ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))

        # self.T_noCC = np.copy(T) #Transmissibility without contour conditions
        self.T_noCC = T.tocsc()

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1

        return T.tocsc(), self.T_noCC.copy()