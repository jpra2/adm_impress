import os.path
import pdb

import scipy

from packs.compositional.IMPEC.pressure_solver import TPFASolver
from packs.data_class import elements_lv0
from packs.utils.test_functions import test_instance
from packs.multiscale.multilevel.multilevel_operators import MultilevelOperators
from packs.utils import constants as ctes
from packs.multiscale.ms_utils.multiscale_functions import multilevel_pressure_solver, print_mesh_volumes_data, update_local_problem
import scipy.sparse as sp
import numpy as np
from packs.adm.non_uniform import monotonic_adm_subds
from packs.multiscale.neuman_local_problems.master_local_solver import MasterLocalSolver
from packs.multiscale.ms_utils.multiscale_functions import update_local_transmissibility
from packs.adm.non_uniform.adm_method_non_nested import AdmNonNested
from packs.multiscale.preprocess.dual_domains import DualSubdomain
from collections.abc import Sequence
from packs.solvers.solvers_trilinos.solvers_tril import solverTril
from packs.solvers.solvers_scipy.solver_sp import SolverSp
import time
from packs.multiscale.ms_solvers import TamsSolverFV
from packs.compositional.IMPEC.global_pressure_solver import GlobalIMPECPressureSolver as Gips
from packs.multiscale.operators.prolongation.AMS.paralell2.paralel_ams_new_2 import MasterLocalOperator

def update_local_parameters(dt, fprop, params):

    param_local_list = ['Vbulk', 'porosity', 'Cf', 'dVtdP', 'P', 'n_components', 'n_phases', 'dVtdk', 'z_centroids', 'xkj_internal_faces', 'Csi_j_internal_faces', 'mobilities_internal_faces', 'pretransmissibility_internal_faces', 'Pcap', 'Vp', 'Vt', 'delta_t', 'g', 'rho_j', 'rho_j_internal_faces', 'n_volumes']

    params.update({
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

    local_params = {i: params[i] for i in param_local_list}
    return local_params


def set_level0_negative_composition(fprop, params):
        
        test1 = fprop.Nk < 0
        if np.any(test1):
            for compositions in fprop.Nk:
                test = compositions < 0
                if np.any(test):
                    vols = np.unique(np.concatenate(params['volumes_to_volumes'][test]))
                    params['level0_negative_composition'][vols] = True

            # data_impress['LEVEL'][
            #     params['level0_negative_composition']
            # ] = 0

class AdmTpfaCompositionalSolver(TPFASolver):

    _kwargs_keys = {
        'get_pressure': [
            'multilevel_data',
            'multilevel_operators'
        ]
    }

    def get_pressure_adm_dep0(self, M, wells, fprop, Pold, delta_t, **kwargs):

        adm_method: AdmNonNested = kwargs.get('adm_method')
        params = kwargs.get('params')
        neumann_subds = kwargs.get('neumann_subds')
        data_impress= kwargs.get('data_impress')
        elements_lv0 = kwargs.get('elements_lv0')
        ml_data = kwargs.get('multilevel_data')
        all_coarse_intersect_faces = np.unique(np.concatenate(ml_data['coarse_intersect_faces_level_1']))
        mlo: MultilevelOperators = kwargs.get('multilevel_operators')
        dual_subdomains: Sequence[DualSubdomain] = kwargs.get('dual_subdomains')
        global_vector_update = kwargs.get('global_vector_update')
        OP_AMS = kwargs.get('OP_AMS')

        # import pdb; pdb.set_tracedescription = 'case18_adm_6k_5000_' # cr = 25()
        # centroids = M.volumes.center[:]
        # dx = centroids[:, 0][elements_lv0['neig_internal_faces'][:,1]] - centroids[:, 0][elements_lv0['neig_internal_faces'][:,0]]
        # dx = dx.min()
        # test = (centroids[:,0] >= centroids[:,0].min() - 0.1) & (centroids[:,0] < centroids.min() + 7*dx + 0.1)
        # vols = elements_lv0['volumes'][test]

        k = 1
        T, T_noCC, T_advec = self.update_transmissibility_adm(M, wells, fprop, delta_t, **kwargs)
        # T *= k
        # T_noCC *= k
        # T_advec *= k
        D = self.update_independent_terms(M, fprop, Pold, wells, delta_t)
        # D *= k

        # OP_AMS, cfs  = mlo.run_paralel_2(
        #     T_noCC,
        #     D,
        #     dual_subdomains,
        #     global_vector_update,
        #     params['diagonal_term'],
        #     OP_AMS,
        #     1
        # )
        # import pdb; pdb.set_trace()

        master_local_operator = kwargs.get('master_local_operator')
        OP_AMS, cfs  = mlo.run_paralel_2(
            T_advec,
            D,
            dual_subdomains,
            global_vector_update,
            # params['diagonal_term'],
            np.zeros(D.shape[0]),
            OP_AMS,
            1,
            master_local_operator
        )
        # import pdb; pdb.set_trace()

        n_levels = 2
        transm_int_fac = np.array(T_noCC[ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        # data_impress['transmissibility'][elements_lv0['internal_faces']] = ctes.pretransmissibility_internal_faces
        data_impress['transmissibility'][elements_lv0['internal_faces']] = transm_int_fac
        data_impress['transmissibility'][elements_lv0['boundary_faces']] = 0

        # #####################
        # preprocessed_primal_objects, critical_groups = monotonic_adm_subds.get_preprossed_monotonic_primal_objects(
        #     data_impress, elements_lv0, mlo.get_prolongation_by_level(1), neumann_subds.neumann_subds, phiK_raz_lim=3)

        # try:
        #     l_groups = np.concatenate([np.repeat(i, len(critical_groups[i])) for i in range(len(critical_groups))])
        #     groups_c = np.concatenate(critical_groups)
        # except:
        #     l_groups = np.array([0])
        #     groups_c = np.array([0])

        # volumes, netasp_array = monotonic_adm_subds.get_monotonizing_volumes(preprocessed_primal_objects, data_impress['transmissibility'])
        # maxs = np.zeros(len(np.unique(volumes)))
        # np.maximum.at(maxs, volumes, netasp_array)
        # data_impress['nfp'][np.unique(volumes)] = maxs

        # neta_lim_finescale = 1
        # vols_orig = monotonic_adm_subds.get_monotonizing_level(l_groups, groups_c, critical_groups, data_impress,
        #                                                        elements_lv0, volumes, netasp_array, neta_lim_finescale)        # adm_method.set_level_wells_3()
        # ########################

        # data_impress['LEVEL'][:] = 1
        # # data_impress['LEVEL'][vols] = 0
        # adm_method.set_level_wells_only()

        # if len(vols_orig) > 0:
        #     adm_method.set_monotonizing_level(vols_orig)
        # # adm_method.set_saturation_level_simple(delta_sat_max)
        # # import pdb; pdb.set_trace()

        # data_impress['LEVEL'][:] = 1
        # self.set_level0_by_composition(data_impress['LEVEL'], fprop.Csi_j, ctes.n_components, 0.1, elements_lv0['neig_internal_faces'], ctes.n_volumes)
        # self.set_level0_wells(data_impress['LEVEL'], adm_method.all_wells_ids, elements_lv0['volumes_face_volumes'], ctes.n_volumes)
        self.set_level0_wells_v2(data_impress['LEVEL'], adm_method.all_wells_ids, ctes.n_volumes, data_impress['GID_1']) # lvel 0 in all coarse ids

        gid_0 = data_impress['GID_0'][data_impress['LEVEL'] == 0]
        gid_1 = data_impress['GID_0'][data_impress['LEVEL'] == 1]
        adm_method.set_adm_mesh_non_nested(v0=gid_0, v1=gid_1, pare=True)
        # data_impress['LEVEL'][vols] = 0

        cfs[data_impress['LEVEL'] == 0] = 0

        for level in range(1, n_levels):
            adm_method.organize_ops_adm(mlo, level)
        ##################

        # import pdb; pdb.set_trace()

        #####
        ##Test with AMS operators
        #####
        prolongation_list = []
        restriction_list = []
        correction_function_list = []
        for level in range(1, n_levels):
            # prolongation_list.append(mlo[mlo.prolongation + str(level)])
            # restriction_list.append(mlo[mlo.restriction + str(level)])
            prolongation_list.append(adm_method['adm_prolongation_level_' + str(level)])
            restriction_list.append(adm_method['adm_restriction_level_' + str(level)])
            correction_function_list.append(np.zeros(adm_method['adm_prolongation_level_' + str(level)].shape[0]))

        # solution, n_active_volumes = multilevel_pressure_solver(
        #     T,
        #     D,
        #     prolongation_list,
        #     restriction_list,
        #     correction_function_list
        # )

        # ##########
        # ## iterate in fine scale
        # # trilinos_solver: solverTril = kwargs.get('trilinos_solver')
        # # solution = trilinos_solver.solve_linear_problem(T, D, x=solution, tolerance=1e-12)
        # scipy_solver: SolverSp = kwargs.get('scipy_solver')
        # tolerance = kwargs.get('tolerance')
        # # solution = scipy_solver.gmres_solver(T, D, x0=solution, tol=tolerance)
        # solution = scipy_solver.conjugate_gradient_solver(T, D, x0=solution, tol=tolerance)
        # ###############

        #####################
        ## Tams solver
        wells_producer = kwargs.get('wells_producer')
        solution, eps, iterations = TamsSolverFV.richardson_solver(
            T,
            D,
            fprop.P,
            restriction_list[0],
            prolongation_list[0],
            res_tol=1e-10,
            x_tol=1e-10,
            max_it = 1000
            # wells_producer = wells_producer
        )
        n_active_volumes = prolongation_list[0].shape[1]
        ##################################
        # from scipy.sparse.linalg import spsolve
        # solution2 = spsolve(T, D)
        # print('WELLS PRESSURE')
        # print(solution[wells_producer])
        # print(solution2[wells_producer])
        # print()
        # import pdb; pdb.set_trace()

        ###########
        # T_noCC /= k
        # T /= k
        # T_advec /= k
        # D /= k

        params.update({
            'active_volumes': n_active_volumes
        })
        Pnew = solution
        # import pdb; pdb.set_trace()

        # # if iterations == 100:
        # loop = kwargs.get('loop')
        # data_impress['pressure'][:] = solution
        # data_impress['So'][:] = fprop.So
        # data_impress['Sg'][:] = fprop.Sg
        # data_impress.update_variables_to_mesh()
        # m = M.core.mb.create_meshset()
        # M.core.mb.add_entities(m, M.core.all_volumes)
        # M.core.mb.write_file('results/test_primals_' + str(loop) + '.vtk', [m])
        # import pdb; pdb.set_trace()


        # P = self.update_pressure(T, D) # OP*padm
        # self.P = self.update_pressure(T, D) # OP*padm
        # error = np.absolute(self.P - solution) / self.P
        # data_impress = M.data
        # data_impress['pressure'][:] = self.P
        # data_impress['ms_pressure'][:] = solution
        # data_impress['pressure_error'][:] = error
        # import pdb; pdb.set_trace()
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

        # Ft_internal_faces_orig = self.update_total_flux_internal_faces(fprop, P) # pressao local

        ############################################
        ## calculo do fluxo
        Ft_internal_faces_adm = self.update_total_flux_internal_faces(M, fprop, solution) # pressao local

        Ft_internal_faces = np.zeros(Ft_internal_faces_adm.shape)
        Ft_internal_faces[:, elements_lv0['remaped_internal_faces'][all_coarse_intersect_faces]] = Ft_internal_faces_adm[:, elements_lv0['remaped_internal_faces'][all_coarse_intersect_faces]]

        kwargs = update_local_parameters(delta_t, fprop, **kwargs)
        # update_local_problem(
        #     neumann_subds.neumann_subds,
        #     T_noCC,
        #     params['diagonal_term'],
        #     self.P,
        #     Ft_internal_faces_orig,
        #     elements_lv0['remaped_internal_faces'],
        #     elements_lv0['volumes'],
        #     elements_lv0['neig_internal_faces'],
        #     all_coarse_intersect_faces,
        #     **kwargs
        # )

        update_local_problem(
            neumann_subds.neumann_subds,
            T_noCC,
            params['diagonal_term'],
            solution,
            Ft_internal_faces_adm,
            elements_lv0['remaped_internal_faces'],
            elements_lv0['volumes'],
            elements_lv0['neig_internal_faces'],
            all_coarse_intersect_faces,
            **kwargs
        )

        # master = MasterLocalSolver(neumann_subds.neumann_subds, ctes.n_volumes)
        # local_solution = master.run()
        # local_solution = master.run_serial()
        # del master
        master_neumann: MasterLocalSolver = kwargs.get('master_neumann')
        local_solution = master_neumann.run()
        # local_solution = master_neumann.run_serial()

        # error2 = np.absolute(self.P - local_solution) / self.P
        # data_impress['verif_po'][:] = local_solution
        # data_impress['verif_rest'][:] = error2
        # data_impress['flux_volumes'][:] = error2

        ft_internal_faces_local_solution = self.update_total_flux_internal_faces(M, fprop, local_solution)
        # from packs.compositional.IMPEC.global_pressure_solver import GlobalIMPECPressureSolver as Gips
        # global_flux = Gips.update_flux(
        #     ft_internal_faces_local_solution,
        #     kwargs['rho_j_internal_faces'],
        #     kwargs['mobilities_internal_faces'],
        #     kwargs['Pcap'],
        #     elements_lv0['neig_internal_faces'],
        #     kwargs['z_centroids'],
        #     kwargs['pretransmissibility_internal_faces'],
        #     kwargs['g'],
        #     kwargs['xkj_internal_faces'],
        #     kwargs['Csi_j_internal_faces'],
        #     kwargs['n_components'],
        #     kwargs['n_volumes']
        # )
        # global_flux2 = Gips.update_flux(
        #     Ft_internal_faces_orig,
        #     kwargs['rho_j_internal_faces'],
        #     kwargs['mobilities_internal_faces'],
        #     kwargs['Pcap'],
        #     elements_lv0['neig_internal_faces'],
        #     kwargs['z_centroids'],
        #     kwargs['pretransmissibility_internal_faces'],
        #     kwargs['g'],
        #     kwargs['xkj_internal_faces'],
        #     kwargs['Csi_j_internal_faces'],
        #     kwargs['n_components'],
        #     kwargs['n_volumes']
        # )

        # # import pdb; pdb.set_trace()

        # data_impress['flux_volumes'][:] = global_flux[0]
        other_faces = np.setdiff1d(elements_lv0['internal_faces'], all_coarse_intersect_faces)
        Ft_internal_faces[:, elements_lv0['remaped_internal_faces'][other_faces]] = ft_internal_faces_local_solution[:, elements_lv0['remaped_internal_faces'][other_faces]]



        # Ft_internal_faces = Ft_internal_faces_adm
        #########################################################

        # data_impress.update_variables_to_mesh()
        # print_mesh_volumes_data(
        #     M,
        #     os.path.join('results', 'test_flux_3.vtk')
        # )

        self.update_flux_wells(fprop, Pnew, wells, delta_t)

        # return self.P, Ft_internal_faces, self.q
        # import pdb; pdb.set_trace()
        # return self.P, Ft_internal_faces_orig, self.q
        # return solution, Ft_internal_faces, self.q

        return Pnew, Ft_internal_faces, self.q

    def update_transmissibility_adm(self, M, wells, fprop, delta_t, params, **kwargs):

        self.t0_internal_faces_prod = fprop.xkj_internal_faces * \
                                      fprop.Csi_j_internal_faces * \
                                      fprop.mobilities_internal_faces

        ''' Transmissibility '''
        t0 = (self.t0_internal_faces_prod).sum(axis = 1)
        t0 = t0 * ctes.pretransmissibility_internal_faces
        # T = np.zeros([ctes.n_volumes, ctes.n_volumes])
        T = sp.csr_matrix((ctes.n_volumes, ctes.n_volumes))

        ####################
        #### new matrix for prolongation operator
        t02 = fprop.mobilities_internal_faces.sum(axis=1) * ctes.pretransmissibility_internal_faces
        lines = np.array([ctes.v0[:, 0], ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        cols = np.array([ctes.v0[:, 1], ctes.v0[:, 0], ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        data = np.array([-t02[:], -t02[:], +t02[:], +t02[:]]).flatten()
        T_advec = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))
        del t02
        #############################

        # Look for a way of doing this not using a loop!!!
        for i in range(ctes.n_components):
            data = np.array([-t0[i,:], -t0[i,:], +t0[i,:], +t0[i,:]]).flatten()

            # Ta = (sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))).toarray()
            Ta = sp.csc_matrix((data, (lines, cols)), shape = (ctes.n_volumes, ctes.n_volumes))
            # T += Ta * self.dVtk[i,:, np.newaxis]
            T += Ta.multiply(self.dVtk[i,:, np.newaxis])

        # T = T * delta_t
        try:
            T *= delta_t
        except ValueError:
            T *= delta_t[0]
        except BaseException as e:
            print(e)
            import pdb; pdb.set_trace()
        # T_advec = T.copy()
        ''' Transmissibility diagonal term '''
        # diag = np.diag((ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        # T += diag

        # diagT = T.diagonal().flatten()
        # diagT += ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP
        # T.setdiag(diaparams['diagonal_term'] = diagonal_termgT)

        # T.setdiag(T.diagonal().flatten() + (ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP))
        diagonal_term = ctes.Vbulk * ctes.porosity * ctes.Cf - self.dVtP
        params.update({
            'diagonal_term': diagonal_term
        })
        T.setdiag(T.diagonal().flatten() + diagonal_term)

        # self.T_noCC = np.copy(T) #Transmissibility without contour conditions
        self.T_noCC = T.tocsc()

        ''' Includding contour conditions '''
        T[wells['ws_p'],:] = 0
        T[wells['ws_p'], wells['ws_p']] = 1

        return T.tocsc(), self.T_noCC.copy(), T_advec

    def set_level0_by_composition(self, level_vector, compositions, n_components, delta_c, adj_internal_faces, n_volumes):


        """define level by maximum value of composition variation
        Args:
            composition ([type]): array with composition of volumes
            n_components ([type]): number of components
            delta_c ([type]): maximum value of composition variation
            adj_internal_faces: internal faces adjacencies
            n_volumes: number of volumes
        """
        level0 = np.full(n_volumes, False, dtype=bool)

        for component in compositions:
            import pdb; pdb.set_trace()
            d_composition = component[adj_internal_faces]
            d_composition = np.absolute(d_composition[:, 1] - d_composition[:, 0])
            test = d_composition >= delta_c
            vols = np.unique(np.concatenate(adj_internal_faces[test]))
            level0[vols] = True

        level_vector[level0] = 0

    def set_level0_wells(self, level_vector, wells_ids, volumes_adjacencies_by_face, n_volumes):
        level0 = np.full(n_volumes, False, dtype=bool)
        adjs = np.unique(np.concatenate(volumes_adjacencies_by_face[wells_ids]))
        level0[wells_ids] = True
        level0[adjs] = True

        level_vector[level0] = 0

    def set_level0_wells_v2(self, level_vector, wells_ids, n_volumes, coarse_gid):
        level0 = np.full(n_volumes, False, dtype=bool)
        coarse_gid_wells = np.unique(coarse_gid[wells_ids])
        for cid in coarse_gid_wells:
            level0[coarse_gid == cid] = True

        level_vector[level0] = 0

    def get_pressure_finescale(self, M, wells, fprop, Pold, delta_t, params, **kwargs):
        # params = kwargs.get('params')
        # k = 1
        T = self.update_transmissibility(M, wells, fprop, delta_t)
        D = self.update_independent_terms(M, fprop, Pold, wells, delta_t)
        # T *= k
        # D *= k
        # tolerance = kwargs.get('tolerance')
        # Pnew = self.update_pressure(T, D)
        scipy_solver: SolverSp = kwargs.get('scipy_solver')
        # import pdb; pdb.set_trace()
        # # solution = scipy_solver.gmres_solver(T, D, x0=params['pressure'], tol=tolerance)
        # solution = scipy_solver.conjugate_gradient_solver(T, D, x0=params['pressure'], tol=tolerance)
        # solution = scipy_solver.LinearCG(T, D, x0=params['pressure'], tol=tolerance)
        solution = scipy_solver.direct_solver(T, D)
        Pnew = solution
        Ft_internal_faces = self.update_total_flux_internal_faces(M, fprop, Pnew)
        self.update_flux_wells(fprop, Pnew, wells, delta_t)

        # # if iterations == 100:
        # loop = kwargs.get('loop')
        # data_impress= kwargs.get('data_impress')
        # data_impress['pressure'][:] = solution
        # data_impress['So'][:] = fprop.So
        # data_impress['Sg'][:] = fprop.Sg
        # data_impress.update_variables_to_mesh()
        # m = M.core.mb.create_meshset()
        # M.core.mb.add_entities(m, M.core.all_volumes)
        # M.core.mb.write_file('results/test_primals_' + str(loop) + '.vtk', [m])
        # import pdb; pdb.set_trace()

        return Pnew, Ft_internal_faces, self.q

    def get_pressure(self, M, wells, fprop, Pold, delta_t, params, **kwargs):
        return self.get_pressure_adm(M, wells, fprop, Pold, delta_t, params, **kwargs)
        # return self.get_pressure_finescale(M, wells, fprop, Pold, delta_t, params, **kwargs)

    def get_pressure_adm(self, M, wells, fprop, Pold, delta_t, params, **kwargs):
        adm_solver = kwargs.get('adm_solver')
        if adm_solver == 'tams':
            return self.get_pressure_adm_tams(M, wells, fprop, Pold, delta_t, params, **kwargs)
        elif adm_solver == 'iterative-finescale':
            return self.get_pressure_adm_ms_direct(M, wells, fprop, Pold, delta_t, **kwargs)
        else:
            raise ValueError('adm_solver error')

    def get_pressure_adm_tams(self, M, wells, fprop, Pold, delta_t, params, **kwargs):
        adm_method: AdmNonNested = params.get('adm_method')
        neumann_subds = params.get('neumann_subds')
        data_impress= params.get('data_impress')
        elements_lv0 = kwargs.get('elements_lv0')
        ml_data = params.get('multilevel_data')
        # all_coarse_intersect_faces = np.unique(np.concatenate(ml_data['coarse_intersect_faces_level_1']))
        all_coarse_intersect_faces = params['all_coarse_intersect_faces_level_1']
        mlo: MultilevelOperators = params.get('multilevel_operators')
        dual_subdomains: Sequence[DualSubdomain] = params.get('dual_subdomains')
        global_vector_update = params.get('global_vector_update')
        OP_AMS = params.get('OP_AMS')
        master_local_operator = params.get('master_local_operator') ## prolongation operator functions
        master_neumann: MasterLocalSolver = params.get('master_neumann') ## neumann problem functions

        T, T_noCC, T_advec = self.update_transmissibility_adm(M, wells, fprop, delta_t, params, **kwargs)
        D = self.update_independent_terms(M, fprop, Pold, wells, delta_t)

        t0 = time.time()
        OP_AMS, cfs = self.get_basis_functions(
            OP_AMS,
            mlo,
            T_advec,
            D,
            dual_subdomains,
            global_vector_update,
            np.zeros(D.shape[0]),
            master_local_operator
        )
        t1 = time.time()
        print('#############################')
        print(f'\n OP_AMS time: {t1 - t0} \n')
        print('#############################')

        # OP_AMS, cfs  = mlo.run_paralel_2(
        #     T_advec,
        #     D,
        #     dual_subdomains,
        #     global_vector_update,
        #     # params['diagonal_term'],
        #     np.zeros(D.shape[0]),
        #     OP_AMS,
        #     1,
        #     master_local_operator
        # )

        n_levels = 2
        transm_int_fac = np.array(T_noCC[ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        # data_impress['transmissibility'][elements_lv0['internal_faces']] = ctes.pretransmissibility_internal_faces
        data_impress['transmissibility'][elements_lv0['internal_faces']] = transm_int_fac
        data_impress['transmissibility'][elements_lv0['boundary_faces']] = 0

        self.set_level0_wells_v2(data_impress['LEVEL'], adm_method.all_wells_ids, ctes.n_volumes, data_impress['GID_1']) # lvel 0 in all coarse ids
        # set_level0_negative_composition(fprop, data_impress, params)
        data_impress['LEVEL'][
                params['level0_negative_composition']
            ] = 0   
        
        gid_0 = data_impress['GID_0'][data_impress['LEVEL'] == 0]
        gid_1 = data_impress['GID_0'][data_impress['LEVEL'] == 1]
        
        adm_method.set_adm_mesh_non_nested(v0=gid_0, v1=gid_1, pare=True)

        cfs[data_impress['LEVEL'] == 0] = 0

        for level in range(1, n_levels):
            adm_method.organize_ops_adm(mlo, level)

        prolongation_list = []
        restriction_list = []
        correction_function_list = []
        for level in range(1, n_levels):
            prolongation_list.append(adm_method['adm_prolongation_level_' + str(level)])
            restriction_list.append(adm_method['adm_restriction_level_' + str(level)])
            correction_function_list.append(np.zeros(adm_method['adm_prolongation_level_' + str(level)].shape[0]))

        #####################
        ## Tams solver
        t0 = time.time()
        wells_producer = kwargs.get('wells_producer')
        solution, eps, iterations = TamsSolverFV.richardson_solver(
            T,
            D,
            fprop.P,
            restriction_list[0],
            prolongation_list[0],
            res_tol=1e-15,
            x_tol=1e-10,
            # wells_producer = wells_producer
            **kwargs
        )
        t1 = time.time()
        print('##################################')
        print(f'\nTAMS time: {t1 - t0}\n')
        print('##################################')
        n_active_volumes = prolongation_list[0].shape[1]
        ##################################
        # import pdb; pdb.set_trace()

        params.update({
            'active_volumes': n_active_volumes,
            'tams_itcounter': iterations
        })
        Pnew = solution

        ## updating params variables for neumann problems
        local_params = update_local_parameters(delta_t, fprop, params)
        local_params.update({
            'remaped_internal_faces': elements_lv0['remaped_internal_faces']
        })

        ############################################
        ## calculo do fluxo multiescala
        t0 = time.time()
        Ft_internal_faces_adm = self.update_total_flux_internal_faces(M, fprop, solution) # pressao local

        Ft_internal_faces = np.zeros(Ft_internal_faces_adm.shape)
        Ft_internal_faces[:, elements_lv0['remaped_internal_faces'][all_coarse_intersect_faces]] = Ft_internal_faces_adm[:, elements_lv0['remaped_internal_faces'][all_coarse_intersect_faces]]

        t2 = time.time()
        # local_solution = self.get_local_solution_serial(Ft_internal_faces_adm, all_coarse_intersect_faces, neumann_subds, T_noCC, params, solution, master_neumann, **kwargs)
        local_solution = solution
        t3 = time.time()
        print('######################################')
        print(f'\nLocal solution time: {t3 - t2} \n')
        print('######################################')

        # local_solution2 = self.get_local_solution_paralell(
        #     Ft_internal_faces_adm,
        #     all_coarse_intersect_faces,
        #     neumann_subds,
        #     T_noCC,
        #     params['diagonal_term'],
        #     local_params,
        #     solution,
        #     master_neumann,
        #     **kwargs
        # )

        ft_internal_faces_local_solution = self.update_total_flux_internal_faces(M, fprop, local_solution)
        # other_faces = np.setdiff1d(elements_lv0['internal_faces'], all_coarse_intersect_faces)
        other_faces = params['other_faces_level_1']
        Ft_internal_faces[:, elements_lv0['remaped_internal_faces'][other_faces]] = ft_internal_faces_local_solution[:, elements_lv0['remaped_internal_faces'][other_faces]]
        t1 = time.time()
        print('######################################')
        print(f'\nFlux time: {t1 - t0} \n')
        print('######################################')
        #########################################################

        # Ft_internal_faces = self.update_total_flux_internal_faces(M, fprop, solution)
        params['internal_faces_velocity'][:] = Ft_internal_faces

        # err_vel = np.absolute(Ft_internal_faces - Ft_internal_faces_adm)
        # print(f'erro maximo na velocidade {err_vel.max()}')
        # print(f'norma do erro na velocidade {np.linalg.norm(err_vel)}')
        # print()

        self.update_flux_wells(fprop, Pnew, wells, delta_t)
        # import pdb; pdb.set_trace()
        return Pnew, Ft_internal_faces, self.q
    
    def get_pressure_adm_iterative_finescale(self, M, wells, fprop, Pold, delta_t, **kwargs):
        adm_method: AdmNonNested = kwargs.get('adm_method')
        params = kwargs.get('params')
        neumann_subds = kwargs.get('neumann_subds')
        data_impress= kwargs.get('data_impress')
        elements_lv0 = kwargs.get('elements_lv0')
        ml_data = kwargs.get('multilevel_data')
        all_coarse_intersect_faces = np.unique(np.concatenate(ml_data['coarse_intersect_faces_level_1']))
        mlo: MultilevelOperators = kwargs.get('multilevel_operators')
        dual_subdomains: Sequence[DualSubdomain] = kwargs.get('dual_subdomains')
        global_vector_update = kwargs.get('global_vector_update')
        OP_AMS = kwargs.get('OP_AMS')

        T, T_noCC, T_advec = self.update_transmissibility_adm(M, wells, fprop, delta_t, **kwargs)

        D = self.update_independent_terms(M, fprop, Pold, wells, delta_t)

        master_local_operator = kwargs.get('master_local_operator')
        OP_AMS, cfs  = mlo.run_paralel_2(
            T_advec,
            D,
            dual_subdomains,
            global_vector_update,
            # params['diagonal_term'],
            np.zeros(D.shape[0]),
            OP_AMS,
            1,
            master_local_operator
        )

        n_levels = 2
        transm_int_fac = np.array(T_noCC[ctes.v0[:, 0], ctes.v0[:, 1]]).flatten()
        # data_impress['transmissibility'][elements_lv0['internal_faces']] = ctes.pretransmissibility_internal_faces
        data_impress['transmissibility'][elements_lv0['internal_faces']] = transm_int_fac
        data_impress['transmissibility'][elements_lv0['boundary_faces']] = 0

        self.set_level0_wells_v2(data_impress['LEVEL'], adm_method.all_wells_ids, ctes.n_volumes, data_impress['GID_1']) # lvel 0 in all coarse ids

        gid_0 = data_impress['GID_0'][data_impress['LEVEL'] == 0]
        gid_1 = data_impress['GID_0'][data_impress['LEVEL'] == 1]
        adm_method.set_adm_mesh_non_nested(v0=gid_0, v1=gid_1, pare=True)
        # data_impress['LEVEL'][vols] = 0

        cfs[data_impress['LEVEL'] == 0] = 0

        for level in range(1, n_levels):
            adm_method.organize_ops_adm(mlo, level)

        prolongation_list = []
        restriction_list = []
        correction_function_list = []
        for level in range(1, n_levels):
            # prolongation_list.append(mlo[mlo.prolongation + str(level)])
            # restriction_list.append(mlo[mlo.restriction + str(level)])
            prolongation_list.append(adm_method['adm_prolongation_level_' + str(level)])
            restriction_list.append(adm_method['adm_restriction_level_' + str(level)])
            correction_function_list.append(np.zeros(adm_method['adm_prolongation_level_' + str(level)].shape[0]))

        solution, n_active_volumes = multilevel_pressure_solver(
            T,
            D,
            prolongation_list,
            restriction_list,
            correction_function_list
        )

        scipy_solver: SolverSp = kwargs.get('scipy_solver')
        tolerance = kwargs.get('tolerance')
        solution = scipy_solver.conjugate_gradient_solver(T, D, x0=solution, tol=tolerance)

        params.update({
            'active_volumes': n_active_volumes
        })
        Pnew = solution

        ##########################
        ## calculo do fluxo
        Ft_internal_faces_adm = self.update_total_flux_internal_faces(M, fprop, solution)
        Ft_internal_faces = Ft_internal_faces_adm
        ##########################

        self.update_flux_wells(fprop, Pnew, wells, delta_t)

        return Pnew, Ft_internal_faces, self.q

    def get_local_solution_serial(self, Ft_internal_faces_adm, all_coarse_intersect_faces, neumann_subds, T_noCC, params, solution, master_neumann: MasterLocalSolver, **kwargs):

        elements_lv0 = kwargs.get('elements_lv0')
        t0 = time.time()
        update_local_problem(
            neumann_subds.neumann_subds,
            T_noCC,
            params['diagonal_term'],
            solution,
            Ft_internal_faces_adm,
            elements_lv0['remaped_internal_faces'],
            elements_lv0['volumes'],
            elements_lv0['neig_internal_faces'],
            all_coarse_intersect_faces,
            params,
            **kwargs
        )
        t1 = time.time()
        print('######################################')
        print(f'\nUpdate neumann local problem time: {t1 - t0} \n')
        print('######################################')

        # local_solution = master_neumann.run() ## paralell
        t0 = time.time()
        local_solution = master_neumann.run_serial()
        t1 = time.time()
        print('######################################')
        print(f'\nSolve inner local solution time: {t1 - t0} \n')
        print('######################################')

        return local_solution

    def get_local_solution_paralell(self, Ft_internal_faces_adm, all_coarse_intersect_faces, neumann_subds, T_noCC, diagonal_term, local_params, solution, master_neumann: MasterLocalSolver, **kwargs):

        elements_lv0 = kwargs.get('elements_lv0')
        v0 = elements_lv0['neig_internal_faces']
        map_internal_faces = elements_lv0['remaped_internal_faces']

        ft_internal_faces_for_prescription = np.zeros(Ft_internal_faces_adm.shape)
        ft_internal_faces_for_prescription[:, map_internal_faces[all_coarse_intersect_faces]] = Ft_internal_faces_adm[:, map_internal_faces[all_coarse_intersect_faces]]

        global_molar_flux_prescription = Gips.update_flux(
            ft_internal_faces_for_prescription,
            local_params['rho_j_internal_faces'],
            local_params['mobilities_internal_faces'],
            local_params['Pcap'],
            v0,
            local_params['z_centroids'],
            local_params['pretransmissibility_internal_faces'],
            local_params['g'],
            local_params['xkj_internal_faces'],
            local_params['Csi_j_internal_faces'],
            local_params['n_components'],
            local_params['n_volumes']
        )

        local_kwargs = {
            'global_molar_flux_prescription': global_molar_flux_prescription,
            'fine_scale_transmissibility_no_bc': T_noCC,
            'map_internal_faces': elements_lv0['remaped_internal_faces'],
            'adm_pressure': solution,
            'diagonal_term': diagonal_term
        }

        local_solution = master_neumann.run2(
            neumann_subds.neumann_subds, local_params, **local_kwargs
        )
        return local_solution

    def get_basis_functions(self, OP_AMS, mlo: MultilevelOperators, T_advec, D, dual_subdomains, global_vector_update, diagonal_to_op, master_local_operator: MasterLocalOperator):
        
        # OP_AMS, cfs  = mlo.run_paralel_2(
        #     T_advec,
        #     D,
        #     dual_subdomains,
        #     global_vector_update,
        #     diagonal_to_op,
        #     OP_AMS,
        #     1,
        #     master_local_operator
        # )

        # #################
        # ### op totally paralell
        # OP_AMS, cfs  = mlo.run_paralel_3(
        #     T_advec,
        #     D,
        #     dual_subdomains,
        #     global_vector_update,
        #     diagonal_to_op,
        #     OP_AMS,
        #     1,
        #     master_local_operator
        # )
        # #################

        #######################
        ## op serial
        if global_vector_update.sum() > 0:
            OP_AMS, cfs  = mlo.run_serial(
                T_advec,
                D,
                dual_subdomains,
                global_vector_update,
                diagonal_to_op,
                OP_AMS,
                1,
                master_local_operator
            )
            ########################

            return OP_AMS, cfs
        else:
            return OP_AMS, np.zeros(OP_AMS.shape[0])

    
        