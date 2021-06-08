import os.path

from packs.compositional.IMPEC.pressure_solver import TPFASolver
from packs.utils.test_functions import test_instance
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.utils import constants as ctes
from packs.multiscale.ms_utils.multiscale_functions import multilevel_pressure_solver, print_mesh_volumes_data, update_local_problem
import scipy.sparse as sp
import numpy as np
from packs.adm.non_uniform import monotonic_adm_subds



class AdmTpfaCompositionalSolver(TPFASolver):
    
    _kwargs_keys = {
        'get_pressure': [
            'multilevel_data',
            'multilevel_operators'
        ]
    }
    
    def get_pressure(self, M, wells, fprop, delta_t, **kwargs):
        adm_method = kwargs.get('adm_method')
        params = kwargs.get('params')
        neumann_subds = kwargs.get('neumann_subds')
        data_impress= kwargs.get('data_impress')
        elements_lv0 = kwargs.get('elements_lv0')

        # test_kwargs_keys(
        #     AdmTpfaCompositionalSolver._kwargs_keys['get_pressure'],
        #     kwargs.keys()
        # )
        
        T, T_noCC = self.update_transmissibility(M, wells, fprop, delta_t)
        D = self.update_independent_terms(M, fprop, wells, delta_t)
        mlo: MultilevelOperators = kwargs.get('multilevel_operators')
        test_instance(mlo, MultilevelOperators)
        # mlo.run(T_noCC, np.zeros(len(D)), np.zeros(len(D)))
        mlo.run(T_noCC, D, np.zeros(len(D)), return_correction_matrix=False)
        n_levels = 2
        data_impress['transmissibility'][elements_lv0['internal_faces']] = ctes.pretransmissibility_internal_faces

        preprocessed_primal_objects, critical_groups = monotonic_adm_subds.get_preprossed_monotonic_primal_objects(
            data_impress, elements_lv0, mlo.get_prolongation_by_level(1), neumann_subds.neumann_subds, phiK_raz_lim=3)

        try:
            l_groups = np.concatenate([np.repeat(i, len(critical_groups[i])) for i in range(len(critical_groups))])
            groups_c = np.concatenate(critical_groups)
        except:
            l_groups = np.array([0])
            groups_c = np.array([0])

        volumes, netasp_array = monotonic_adm_subds.get_monotonizing_volumes(preprocessed_primal_objects, data_impress['transmissibility'])
        maxs = np.zeros(len(np.unique(volumes)))
        np.maximum.at(maxs, volumes, netasp_array)
        data_impress['nfp'][np.unique(volumes)] = maxs


        neta_lim_finescale = 1
        vols_orig = monotonic_adm_subds.get_monotonizing_level(l_groups, groups_c, critical_groups, data_impress,
                                                               elements_lv0, volumes, netasp_array, neta_lim_finescale)        # adm_method.set_level_wells_3()
        data_impress['LEVEL'][:] = 1
        adm_method.set_level_wells_only()

        if len(vols_orig) > 0:
            adm_method.set_monotonizing_level(vols_orig)
        # adm_method.set_saturation_level_simple(delta_sat_max)

        gid_0 = data_impress['GID_0'][data_impress['LEVEL'] == 0]
        gid_1 = data_impress['GID_0'][data_impress['LEVEL'] == 1]
        adm_method.set_adm_mesh_non_nested(v0=gid_0, v1=gid_1, pare=True)

        for level in range(1, n_levels):
            adm_method.organize_ops_adm(mlo, level)

        # import pdb; pdb.set_trace()

        #####
        ##Test with AMS operators
        #####
        prolongation_list = []
        restriction_list = []
        for level in range(1, n_levels):
            # prolongation_list.append(mlo[mlo.prolongation + str(level)])
            # restriction_list.append(mlo[mlo.restriction + str(level)])
            prolongation_list.append(adm_method['adm_prolongation_level_' + str(level)])
            restriction_list.append(adm_method['adm_restriction_level_' + str(level)])

        solution = multilevel_pressure_solver(
            T,
            D,
            prolongation_list,
            restriction_list
        )

        update_local_problems(neumann_subds.neumann_subds, T_noCC, fprop)



        import pdb; pdb.set_trace()

        self.P = self.update_pressure(T, D, fprop) # OP*padm
        error = np.absolute(self.P - solution) / self.P
        data_impress = M.data
        data_impress['pressure'][:] = self.P
        data_impress['ms_pressure'][:] = solution
        data_impress['pressure_error'][:] = error
        data_impress.update_variables_to_mesh()
        # m1 = M.core.mb.create_meshset()
        # M.    core.mb.add_entities(m1, M.core.all_volumes)
        # M.core.mb.write_file('results/test_comp_1.vtk', [m1])
        # insert_prolongation_operator_impress(
        #     M,
        #     prolongation_list[0],
        #     1,
        #     data_impress['GID_0'],
        #     data_impress['GID_0'],
        # )
        #
        # insert_restriction_operator_impress(
        #     M,
        #     restriction_list[0],
        #     1,
        #     data_impress['GID_0'],
        #     data_impress['GID_0'],
        # )

        print_mesh_volumes_data(
            M,
            os.path.join('results', 'prolongation_level_1.vtk')
        )



        import pdb; pdb.set_trace()

        Ft_internal_faces = self.update_total_flux_internal_faces(M, fprop) # pressao local
        self.update_flux_wells(fprop, wells, delta_t)
        params['dVtdP'] = AdmTpfaCompositionalSolver.dVtP
        params['dVtdNk'] = AdmTpfaCompositionalSolver.dVtk
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