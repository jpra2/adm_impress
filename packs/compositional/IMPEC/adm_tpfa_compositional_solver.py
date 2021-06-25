import os.path

from packs.compositional.IMPEC.pressure_solver import TPFASolver
from packs.utils.test_functions import test_instance
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.utils import constants as ctes
from packs.multiscale.ms_utils.multiscale_functions import multilevel_pressure_solver, print_mesh_volumes_data, update_local_problem
import scipy.sparse as sp
import numpy as np
from packs.adm.non_uniform import monotonic_adm_subds
from packs.multiscale.neuman_local_problems.master_local_solver import MasterLocalSolver
from packs.multiscale.ms_utils.multiscale_functions import update_local_transmissibility

def update_local_parameters(dt, fprop, **kwargs):
    params = kwargs.get('params')
    kwargs.update({
        'dVtdP': params['dVtdP'],
        'dVtdk': params['dVtdNk'],
        'P': fprop.P,
        'xkj_internal_faces': fprop.xkj_internal_faces,
        'Csi_j_internal_faces': fprop.Csi_j_internal_faces,
        'mobilities_internal_faces': fprop.mobilities_internal_faces,
        'Pcap': fprop.Pcap,
        'Vp': fprop.Vp,
        'Vt': fprop.Vt,
        'delta_t': dt,
        'rho_j': fprop.rho_j,
        'rho_j_internal_faces': fprop.rho_j_internal_faces
    })
    return kwargs




class AdmTpfaCompositionalSolver(TPFASolver):
    
    _kwargs_keys = {
        'get_pressure': [
            'multilevel_data',
            'multilevel_operators'
        ]
    }
    
    def get_pressure(self, M, wells, fprop, delta_t, Pold, **kwargs):

        adm_method = kwargs.get('adm_method')
        params = kwargs.get('params')
        neumann_subds = kwargs.get('neumann_subds')
        data_impress= kwargs.get('data_impress')
        elements_lv0 = kwargs.get('elements_lv0')
        ml_data = kwargs.get('multilevel_data')
        all_coarse_intersect_faces = np.unique(np.concatenate(ml_data['coarse_intersect_faces_level_1']))
        mlo: MultilevelOperators = kwargs.get('multilevel_operators')

        # test_kwargs_keys(
        #     AdmTpfaCompositionalSolver._kwargs_keys['get_pressure'],
        #     kwargs.keys()
        # )
        
        T, T_noCC = self.update_transmissibility(M, wells, fprop, delta_t, **kwargs)
        D = self.update_independent_terms(M, fprop, Pold, wells, delta_t)
        # D2 = GlobalIMPECPressureSolver.mount_independent_term(
        #     ctes.Vbulk,
        #     ctes.porosity,
        #     ctes.Cf,
        #     self.dVtP,
        #     fprop.P,
        #     ctes.n_volumes,
        #     ctes.n_components,
        #     3,
        #     ctes.v0,
        #     self.dVtk,
        #     ctes.z,
        #     fprop.xkj_internal_faces,
        #     fprop.Csi_j_internal_faces,
        #     fprop.mobilities_internal_faces,
        #     ctes.pretransmissibility_internal_faces,
        #     fprop.Pcap,
        #     fprop.Vp,
        #     fprop.Vt,
        #     wells['ws_q'],
        #     wells['values_q'],
        #     delta_t,
        #     ctes.g,
        #     wells['ws_p'],
        #     wells['values_p'],
        #     ctes.bhp_ind,
        #     fprop.rho_j,
        #     fprop.rho_j_internal_faces
        # )
        # print(np.allclose(D, D2))
        #
        # import pdb; pdb.set_trace()

        # test_instance(mlo, MultilevelOperators)
        # mlo.run(T_noCC, np.zeros(len(D)), np.zeros(len(D)))
        mlo.run(T_noCC, D, np.zeros(len(D)), return_correction_matrix=False)
        n_levels = 2
        transm_int_fac = np.array(T_noCC[ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        # data_impress['transmissibility'][elements_lv0['internal_faces']] = ctes.pretransmissibility_internal_faces
        data_impress['transmissibility'][elements_lv0['internal_faces']] = transm_int_fac
        data_impress['transmissibility'][elements_lv0['boundary_faces']] = 0

        ######################
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
        ##################

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

        self.P = self.update_pressure(T, D) # OP*padm
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

        # print_mesh_volumes_data(
        #     M,
        #     os.path.join('results', 'prolongation_level_1.vtk')
        # )

        Ft_internal_faces_orig = self.update_total_flux_internal_faces(fprop, self.P) # pressao local
        Ft_internal_faces_adm = self.update_total_flux_internal_faces(fprop, solution) # pressao local
        Ft_internal_faces = np.zeros(Ft_internal_faces_adm.shape)
        Ft_internal_faces[:, elements_lv0['remaped_internal_faces'][all_coarse_intersect_faces]] = Ft_internal_faces_adm[:, elements_lv0['remaped_internal_faces'][all_coarse_intersect_faces]]

        kwargs = update_local_parameters(delta_t, fprop, **kwargs)
        update_local_problem(
            neumann_subds.neumann_subds,
            T_noCC,
            params['diagonal_term'],
            self.P,
            Ft_internal_faces_orig,
            elements_lv0['remaped_internal_faces'],
            elements_lv0['volumes'],
            elements_lv0['neig_internal_faces'],
            all_coarse_intersect_faces,
            **kwargs
        )
        master = MasterLocalSolver(neumann_subds.neumann_subds, ctes.n_volumes)
        local_solution = master.run()
        del master
        # import pdb; pdb.set_trace()

        error2 = np.absolute(self.P - local_solution) / self.P
        data_impress['verif_po'][:] = local_solution
        data_impress['verif_rest'][:] = error2
        # data_impress['flux_volumes'][:] = error2

        ft_internal_faces_local_solution = self.update_total_flux_internal_faces(fprop, local_solution)
        other_faces = np.setdiff1d(elements_lv0['internal_faces'], all_coarse_intersect_faces)
        Ft_internal_faces[:, elements_lv0['remaped_internal_faces'][other_faces]] = ft_internal_faces_local_solution[:, elements_lv0['remaped_internal_faces'][other_faces]]

        data_impress.update_variables_to_mesh()
        print_mesh_volumes_data(
            M,
            os.path.join('results', 'prolongation_level_1.vtk')
        )
        import pdb; pdb.set_trace()

        self.update_flux_wells(fprop, wells, delta_t, solution)

        return self.P, Ft_internal_faces, self.q
    
    def update_transmissibility(self, M, wells, fprop, delta_t, **kwargs):
        params = kwargs.get('params')
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

        # T.setdiag(T.diagonal().flatten() + (ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        diagonal_term = ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP
        params['diagonal_term'] = diagonal_term
        T.setdiag(T.diagonal().flatten() + diagonal_term)

        # self.T_noCC = np.copy(T) #Transmissibility without contour conditions
        self.T_noCC = T.tocsc()

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1

        return T.tocsc(), self.T_noCC.copy()

    def update_total_flux_internal_faces(self, fprop, pressure):
        # Pot_hid = self.P + fprop.Pcap
        Pot_hid = pressure + fprop.Pcap
        Pot_hidj = Pot_hid[:,ctes.v0[:,0]]
        Pot_hidj_up = Pot_hid[:,ctes.v0[:,1]]
        z = ctes.z[ctes.v0[:,0]]
        z_up = ctes.z[ctes.v0[:,1]]
        Ft_internal_faces = - np.sum(fprop.mobilities_internal_faces
            * ctes.pretransmissibility_internal_faces * ((Pot_hidj_up - Pot_hidj) -
            ctes.g * fprop.rho_j_internal_faces * (z_up - z)), axis = 1)
        return Ft_internal_faces

    def update_flux_wells(self, fprop, wells, delta_t, pressure):
        wp = wells['ws_p']

        if len(wp) >= 1:
            well_term = (self.T_noCC[wp, :] @ pressure - self.pressure_term[wp] +
                         self.volume_term[wp]) / delta_t + self.capillary_term[wp] + \
                         self.gravity_term[wp]
            mob_ratio = fprop.mobilities[:, :, wp] / \
                        np.sum(fprop.mobilities[:, :, wp], axis=1)
            self.q[:, wp] = np.sum(fprop.xkj[:, :, wp] * mob_ratio *
                                   fprop.Csi_j[:, :, wp] * well_term, axis=1)
            fprop.q_phase = mob_ratio * well_term