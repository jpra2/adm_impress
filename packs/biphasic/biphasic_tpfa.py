from ..pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from ..directories import data_loaded
from ..utils import relative_permeability
from ..solvers.solvers_scipy.solver_sp import SolverSp
from ..solvers.solvers_trilinos.solvers_tril import solverTril
import os
from .. import directories as direc
import numpy as np
import scipy.sparse as sp
import time
from .biphasic_properties import biphasicProperties
import math
from ..errors.err import MaxLoopIterationError
from ..data_class.structured_mesh_properties import StructuredMeshProperties
from .tests_general import testsGeneral
from ..utils.utils_old import get_box

class BiphasicTpfa(FineScaleTpfaPressureSolver, biphasicProperties, testsGeneral):

    def __init__(self, M, data_impress, elements_lv0, wells, data_name: str='BiphasicTpfa.npz'):
        load = data_loaded['load_biphasic_data']
        self.load = load
        super().__init__(data_impress, elements_lv0, wells, data_name=data_name, load=load)
        self.biphasic_data = data_loaded['biphasic_data']
        self.relative_permeability = getattr(relative_permeability, self.biphasic_data['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.V_total = (data_impress['volume']*data_impress['poro']).sum()
        self.max_contador_vtk = len(self.biphasic_data['vpis_para_gravar_vtk'])
        self.delta_sat_max = 0.1
        self.lim_flux_w = 9e-8
        self.name_current_biphasic_results = os.path.join(direc.flying, 'current_biphasic_results.npy')
        self.name_all_biphasic_results = os.path.join(direc.flying, 'all_biphasic_results_')
        self.mesh_name = os.path.join(direc.flying, 'biphasic_')
        self.all_biphasic_results = self.get_empty_current_biphasic_results()
        self.mesh = M
        self.solver = SolverSp()
        self.meshset_volumes = M.core.mb.create_meshset()
        M.core.mb.add_entities(self.meshset_volumes, M.core.all_volumes)
        self.meshset_faces = M.core.mb.create_meshset()
        M.core.mb.add_entities(self.meshset_faces, M.core.all_faces)

        # self.solver = solverTril()

        if not load:
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.delta_t = 0.0
            self.contador_vtk = 0
            self.update_relative_permeability()
            self.update_mobilities()
            self.update_transmissibility_ini()
            self.export_to_npz()
            self.update_current_biphasic_results(0.0)
            self.save_infos()
        else:
            self.load_infos()

    def update_gama(self):

        if self.gravity:

            gama_w = self.biphasic_data['gama_w']
            gama_o = self.biphasic_data['gama_o']

            gama_w = np.repeat(gama_w, self.n_volumes)
            gama_o = np.repeat(gama_o, self.n_volumes)

            saturations = self.data_impress['saturation']

            gama = gama_w*saturations + gama_o*(1-saturations)

            self.data_impress['gama'] = gama

    def update_relative_permeability(self):

        krw, kro = self.relative_permeability(self.data_impress['saturation'])

        self.data_impress['krw'] = krw
        self.data_impress['kro'] = kro

    def update_mobilities(self):

        n = self.n_volumes
        # lambda_w = self.data_impress['krw']/self.biphasic_data['mi_w']
        # lambda_o = self.data_impress['kro']/self.biphasic_data['mi_o']
        # lambda_t = lambda_w + lambda_o
        # fw_vol = lambda_w/lambda_t
        lambda_w = self.lambda_w_volumes
        lambda_o = self.lambda_o_volumes
        lambda_t = self.lambda_t_volumes
        fw_vol = self.fw_volumes

        self.data_impress['lambda_w'] = lambda_w
        self.data_impress['lambda_o'] = lambda_o
        self.data_impress['lambda_t'] = lambda_t
        self.data_impress['fw_vol'] = fw_vol

    def update_transmissibility_ini(self):

        pretransmissibility = self.data_impress['pretransmissibility'].copy()
        internal_faces = self.elements_lv0['internal_faces']
        b_faces = self.elements_lv0['boundary_faces']
        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        vols_viz_boundary_faces = self.elements_lv0['neig_boundary_faces']
        self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        saturations = self.data_impress['saturation']
        delta_sat = saturations[vols_viz_internal_faces[:,0]] - saturations[vols_viz_internal_faces[:,1]]
        pos = delta_sat >= 0
        self._data['upwind_identificate'][pos, 0] = np.full(pos.sum(), True, dtype=bool)
        pos = ~pos
        self._data['upwind_identificate'][pos, 1] = np.full(pos.sum(), True, dtype=bool)
        self._data['upwind_identificate_o'] = ~self._data['upwind_identificate'].copy()
        self.upwind_wells()

        self.visualize_upwind_vec()

        # ws_inj = self.wells['ws_inj']
        #
        # v0 = vols_viz_internal_faces[:, 0]
        # v1 = vols_viz_internal_faces[:, 1]
        # ids = np.arange(len(v0))
        #
        # idv0 = np.array([], dtype=np.int64)
        # idv1 = idv0.copy()
        # for w in ws_inj:
        #     vv0 = ids[v0 == w]
        #     vv1 = ids[v1 == w]
        #     idv0 = np.append(idv0, vv0)
        #     idv1 = np.append(idv1, vv1)
        #
        # idv0 = np.array(idv0).flatten()
        # idv1 = np.array(idv1).flatten()
        # idsv = np.union1d(idv0, idv1)
        # ids_fora = np.setdiff1d(ids, idsv)
        #
        # self._data['upwind_identificate'][idv0, 0] = np.full(len(idv0), True, dtype=bool)
        # self._data['upwind_identificate'][idv1, 1] = np.full(len(idv1), True, dtype=bool)
        # self._data['upwind_identificate'][ids_fora, 0] = np.full(len(ids_fora), True, dtype=bool)

        # total_mobility_internal_faces = np.zeros(len(internal_faces))
        # fw_internal_faces = total_mobility_internal_faces.copy()

        # total_mobility_internal_faces[idv0] = lambda_t[v0[idv0]]
        # total_mobility_internal_faces[idv1] = lambda_t[v1[idv1]]
        # total_mobility_internal_faces[ids_fora] = (lambda_t[v0[ids_fora]] + lambda_t[v1[ids_fora]])/2

        total_mobility_internal_faces = self.lambda_t_internal_faces

        # fw_internal_faces[idv0] = fw_vol[v0[idv0]]
        # fw_internal_faces[idv1] = fw_vol[v1[idv1]]

        fw_internal_faces = self.fw_internal_faces
        # gama_internal_faces[idv0] = gama[v0[idv0]]
        # gama_internal_faces[idv1] = gama[v1[idv1]]
        # fw_internal_faces[ids_fora] = (fw_vol[v0[ids_fora]] + fw_vol[v1[ids_fora]])/2

        total_mobility_faces = np.zeros(len(pretransmissibility))
        fw_faces = total_mobility_faces.copy()

        total_mobility_faces[internal_faces] = total_mobility_internal_faces
        total_mobility_faces[b_faces] = self.lambda_t_volumes[vols_viz_boundary_faces]

        fw_faces[internal_faces] = fw_internal_faces
        fw_faces[b_faces] = self.fw_volumes[vols_viz_boundary_faces]
        # gama_faces[internal_faces] = gama_internal_faces
        # gama_faces[b_faces] = gama[vols_viz_boundary_faces]

        transmissibility = pretransmissibility*total_mobility_faces

        self.data_impress['transmissibility'] = transmissibility.copy()
        self.data_impress['lambda_t_faces'] = total_mobility_faces.copy()
        self.data_impress['fw_faces'] = fw_faces.copy()
        # self.data_impress['gama_faces'] = gama_faces.copy()

    def get_empty_current_biphasic_results(self):

        return [np.array(['loop', 'delta_t [s]', 'simulation_time [s]',
            'oil_production [m3/s]', 'water_production [m3/s]', 't [s]', 'wor', 'vpi', 'contador_vtk'])]

    def load_infos(self):
        self.current_biphasic_results = np.load(self.name_current_biphasic_results)
        self.loop = int(self.current_biphasic_results[0])
        self.t = self.current_biphasic_results[5]
        self.vpi = self.current_biphasic_results[7]
        self.contador_vtk = self.current_biphasic_results[8]

    def update_flux_w_and_o_volumes(self) -> None:

        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        v0 = vols_viz_internal_faces
        internal_faces = self.elements_lv0['internal_faces']
        u_normal_internal_faces = self.data_impress['u_normal'][internal_faces]
        total_flux_faces = self.data_impress['flux_faces']
        fw_faces = self.data_impress['fw_faces']
        fw_vol = self.data_impress['fw_vol']
        flux_volumes = self.data_impress['flux_volumes']
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']
        x = self.data_impress['pressure']
        grav_source_term_water_volumes = self._data['grav_source_term_water_volumes']

        areas_internal_faces = self.data_impress['area'][internal_faces]
        k_harm_internal_faces = self.data_impress['k_harm'][internal_faces]
        dh_internal_faces = self.data_impress['dist_cent'][internal_faces]

        # ps0 = x[v0[:, 0]]
        # ps1 = x[v0[:, 1]]

        # flux_w_faces = fw_faces*total_flux_faces
        # flux_w_internal_faces2 = flux_w_faces[internal_faces]
        # flux_w_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*self.lambda_w_internal_faces/dh_internal_faces - self._data['grav_source_term_water_faces'][internal_faces])
        # flux_w_internal_faces = self.flux_w_internal_faces

        # self.data_impress['flux_press_w_faces_vec'][internal_faces] = (-((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*self.lambda_w_internal_faces/dh_internal_faces)).reshape(len(internal_faces), 1)*u_normal_internal_faces
        self.data_impress['flux_press_w_faces_vec'][internal_faces] = self.flux_press_w_internal_faces.reshape(len(internal_faces), 1)*u_normal_internal_faces

        # lambda_o_internal_faces = self.data_impress['lambda_o'][v0[self._data['upwind_identificate']]]
        # flux_o_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*self.lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_faces'][internal_faces]))
        # flux_o_internal_faces = self.flux_o_internal_faces

        self.data_impress['flux_press_o_faces_vec'][internal_faces] = self.flux_press_o_internal_faces.reshape(len(internal_faces), 1)*u_normal_internal_faces

        # flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        # soma = flux_w_internal_faces + flux_o_internal_faces
        # vv = abs(soma - total_flux_faces[internal_faces])


        # verif = np.allclose(flux_w_internal_faces, flux_w_internal_faces2)
        # import pdb; pdb.set_trace()

        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        # flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        #
        # # import pdb; pdb.set_trace()
        #
        # flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*self.fw_volumes[ws_prod]
        # flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*self.fw_volumes[ws_inj]
        # flux_w_volumes = self.flux_w_volumes

        # flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        # lambda_o_internal_faces = self.data_impress['lambda_t'][v0[self._data['upwind_identificate']]] - self.data_impress['lambda_w'][v0[self._data['upwind_identificate']]]
        # flux_o_internal_faces2 = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_internal_faces']))

        # lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        # cols = np.repeat(0, len(lines))
        # data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        # flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        # flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - self.fw_volumes[ws_prod])

        # flux_o_volumes = self.flux_o_volumes
        #
        # flux_w_faces = np.zeros(len(self.data_impress['flux_w_faces']))
        # flux_w_faces[internal_faces] = flux_w_internal_faces

        self.data_impress['flux_w_faces'] = self.flux_w_faces
        self.data_impress['flux_o_faces'] = self.flux_o_faces
        self.data_impress['flux_w_volumes'] = self.flux_w_volumes
        self.data_impress['flux_o_volumes'] = self.flux_o_volumes

        # u_normal = self.data_impress['u_normal']
        # flux_w_vec_internal_faces = u_normal[internal_faces]*self.data_impress['flux_w_faces'][internal_faces].reshape([len(internal_faces), 1])
        # flux_o_vec_internal_faces = u_normal[internal_faces]*flux_o_internal_faces.reshape([len(internal_faces), 1])
        # self.data_impress['flux_w_faces_vec'] = np.zeros(self.data_impress['flux_w_faces_vec'].shape)
        # self.data_impress['flux_o_faces_vec'] = np.zeros(self.data_impress['flux_w_faces_vec'].shape)
        # self.data_impress['flux_w_faces_vec'][internal_faces] = flux_w_vec_internal_faces
        # self.data_impress['flux_o_faces_vec'][internal_faces] = flux_o_vec_internal_faces
        self.data_impress['flux_w_faces_vec'] = self.flux_w_faces_vec
        self.data_impress['flux_o_faces_vec'] = self.flux_o_faces_vec

    def update_delta_t_dep0(self):
        ###
        ## de acordo com o fluxo nos volumes
        ###

        flux_volumes = np.absolute(self.data_impress['flux_volumes'])
        phis = self.data_impress['poro']
        volume = self.data_impress['volume']
        delta_t = (self.biphasic_data['cfl']*(volume*phis)/flux_volumes).min()
        return delta_t

    def update_delta_t(self):
        ###
        ## de acordo com o fluxo de agua nos volumes
        ###

        deltas_t = []
        deltas_t.append(self.update_delta_t_for_delta_sat_max())
        deltas_t.append(self.update_delta_t_dep0())
        deltas_t.append(self.update_delta_t_new())

        # import pdb; pdb.set_trace()

        self.delta_t = min(deltas_t)

        # self.delta_t = flux_w_volumes/(volumes*phis)
        # self.delta_t = (self.biphasic_data['cfl']*(volume*phis)/flux_volumes).min()

    def update_delta_t_for_delta_sat_max(self):
        flux_w_volumes = self.data_impress['flux_w_volumes']
        phis = self.data_impress['poro']
        volume = self.data_impress['volume']
        test = flux_w_volumes/(volume*phis)
        test = np.absolute(test[np.absolute(test)>0]).max()
        delta_t = self.delta_sat_max/test
        return delta_t

    def update_delta_t_new(self):

        ####
        # de acordo com as faces
        ####

        lim_ds = 1e-10

        phis = self.data_impress['poro']
        volume = self.data_impress['volume']
        velocity_faces = np.absolute(self.data_impress['velocity_faces'])
        internal_faces = self.elements_lv0['internal_faces']
        fw_vol = self.data_impress['fw_vol']
        saturation = self.data_impress['saturation']
        viz_int = self.elements_lv0['neig_internal_faces']
        dists = self.data_impress['dist_cent']

        dists_int = dists[internal_faces]
        vel_internal_faces = np.linalg.norm(velocity_faces[internal_faces], axis=1)

        v0 = viz_int[:, 0]
        v1 = viz_int[:, 1]
        df = np.absolute(fw_vol[v1] - fw_vol[v0])
        ds = np.absolute(saturation[v1] - saturation[v0])
        ids = np.arange(len(ds))
        ids_2 = ids[ds > lim_ds]

        dists_int = dists_int[ids_2]
        vel_internal_faces = vel_internal_faces[ids_2]
        df = df[ids_2]
        ds = ds[ids_2]
        dfds = df/ds
        delta_t = (self.biphasic_data['cfl']*(dists_int/(vel_internal_faces*dfds))).min()
        return delta_t

    def update_saturation(self):
        cont = 0
        max_loops = 100
        self.data_impress['saturation_last'] = self.data_impress['saturation'].copy()
        verif = -1
        while verif != 0:
            verif = self.update_sat()
            if verif == 1:
                self.reduce_delta_t()

            if cont == max_loops:
                raise MaxLoopIterationError('Loop maximo atingido')

            cont += 1

    def update_sat(self):

        # import pdb; pdb.set_trace()

        saturations0 = self.data_impress['saturation'].copy()
        saturations = saturations0.copy()
        ids = np.arange(len(saturations))

        fw_volumes = -self.data_impress['flux_w_volumes']
        volumes = self.data_impress['volume']
        phis = self.data_impress['poro']

        # import pdb; pdb.set_trace()

        # ids_2 = ids[fw_volumes < 0]
        # if len(ids_2) > 0:
        #     self.data_impress['test_fw'][ids_2] = np.repeat(1.0, len(ids_2))
        #     self.data_impress.update_variables_to_mesh()
        #     self.mesh.core.print(file='test', extension='.vtk', config_input="input_cards/print_settings0.yml")
        #     import pdb; pdb.set_trace()
        #     self.data_impress['test_fw'] = np.repeat(0.0, len(self.data_impress['test_fw']))
        #     self.data_impress.update_variables_to_mesh(['test_fw'])

        ###########################
        ## teste
        test = ids[(saturations < 0) | (saturations > 1)]
        if len(test) > 0:
            import pdb; pdb.set_trace()
            raise ValueError(f'valor errado da saturacao {saturations[test]}')
        del test
        ###########################

        #####
        ### correct flux
        # test = ids[(fw_volumes < 0) & (np.absolute(fw_volumes) < self.lim_flux_w)]
        # if len(test) > 0:
        #     import pdb; pdb.set_trace()
        #     fw_volumes[test] = np.zeros(len(test))

        #####

        # ##########################
        # # teste
        # test = ids[fw_volumes < -self.lim_flux_w]
        # if len(test) > 0:
        #     import pdb; pdb.set_trace()
        #     raise ValueError(f'fluxo negativo de agua {fw_volumes[test]}')
        # ##########################


        ids_var = ids[np.absolute(fw_volumes) > self.lim_flux_w]

        ###################
        ## teste variacao do fluxo de agua
        if len(ids_var) == 0:
            import pdb; pdb.set_trace()
        ###################

        fw_volumes = fw_volumes[ids_var]
        volumes = volumes[ids_var]
        phis = phis[ids_var]
        sats = saturations[ids_var]

        sats += (fw_volumes*self.delta_t)/(phis*volumes)

        delta_sat = sats - saturations[ids_var]
        ids2 = np.arange(len(delta_sat))

        #############
        ## teste variacao maxima de saturacao
        test = ids2[delta_sat > self.delta_sat_max+0.000001]
        if len(test) > 0:
            return 1
        del test
        ##############

        saturations[ids_var] = sats

        # if np.allclose(saturations, saturations0):
        #     import pdb; pdb.set_trace()

        ########################
        # import pdb; pdb.set_trace()
        test = ids[(saturations < 0) | (saturations > 1)]
        if len(test) > 0:
            self.data_impress['saturation'] = saturations
            self.data_impress.update_variables_to_mesh()
            self.mesh.core.print(file='results/test_', extension='.vtk', config_input="input_cards/print_settings0.yml")
            import pdb; pdb.set_trace()

            raise ValueError(f'valor errado da saturacao {saturations[test]}')
        del test
        #########################

        min_sat = saturations.min()
        max_sat = saturations.max()

        if min_sat < self.biphasic_data['Swc'] or max_sat > 1-self.biphasic_data['Sor']:
            return 1
            # raise ValueError(f'\nprint max_sat: {max_sat} ; min_sat: {min_sat}\n')

        self.data_impress['saturation'] = saturations

        return 0

    def reduce_delta_t(self):
        d = self.delta_t
        self.delta_t *= 1/2
        # print(f'\nreducing delta_t: {d} -> {self.delta_t} \n')

    def update_t(self):
        self.t += self.delta_t

    def update_vpi(self):

        flux_total_inj = np.absolute(self.data_impress['flux_volumes'][self.wells['ws_inj']])
        self.vpi += (flux_total_inj.sum()*self.delta_t)/self.V_total

    def update_upwind_phases_old0(self):

        internal_faces = self.elements_lv0['internal_faces']
        flux_o_internal_faces = self.data_impress['flux_o_faces'][internal_faces]
        flux_w_internal_faces = self.data_impress['flux_w_faces'][internal_faces]

        pos = flux_w_internal_faces >= 0
        self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        self._data['upwind_identificate'][pos, 0] = np.full(pos.sum(), True, dtype=bool)
        pos = ~pos
        self._data['upwind_identificate'][pos, 1] = np.full(pos.sum(), True, dtype=bool)

        pos = self.flux_o_internal_faces >= 0
        self._data['upwind_identificate_o'] = np.full((len(internal_faces), 2), False, dtype=bool)
        self._data['upwind_identificate_o'][pos, 0] = np.full(pos.sum(), True, dtype=bool)
        pos = ~pos
        self._data['upwind_identificate_o'][pos, 1] = np.full(pos.sum(), True, dtype=bool)

    def update_upwind_phases_old1(self):
        '''
            paper Starnoni
        '''
        d1 = 0

        internal_faces = self.elements_lv0['internal_faces']
        # flux_o_internal_faces = self.data_impress['flux_o_faces'][internal_faces]
        flux_w_internal_faces = self.data_impress['flux_w_faces'][internal_faces]
        q_sigma_internal_faces = self.flux_sigma_internal_faces
        # q_sigma_internal_faces = self.data_impress['flux_faces'][internal_faces]
        ak = self.data_impress['area'][internal_faces]*self.data_impress['k_harm'][internal_faces]

        pos_sigma = q_sigma_internal_faces > 0 - d1
        pos_w = flux_w_internal_faces > 0 - d1
        and_pos = pos_w & pos_sigma
        self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        self._data['upwind_identificate'][pos_w & pos_sigma, 0] = True
        self._data['upwind_identificate'][~(pos_w & pos_sigma), 1] = True
        self._data['upwind_identificate_o'] = self._data['upwind_identificate'].copy()

        flux_w_internal_faces = -(self.grad_p_internal_faces*ak*self.lambda_w_internal_faces - self.flux_grav_w_internal_faces)
        flux_o_internal_faces = -(self.grad_p_internal_faces*ak*self.lambda_o_internal_faces - self.flux_grav_o_internal_faces)

        flux_total = flux_w_internal_faces + flux_o_internal_faces
        pos_flux_total = flux_total > 0 - d1
        self._data['upwind_identificate_o'][~(pos_flux_total & pos_sigma)] = ~self._data['upwind_identificate_o'][~(pos_flux_total & pos_sigma)]

        self.visualize_upwind_vec()

    def upwind_wells(self):

        wells_inj = self.wells['ws_inj']
        if len(wells_inj) > 0:
            set_wells_inj = set(wells_inj)
            faces = np.unique(np.concatenate(self.elements_lv0['volumes_face_faces'][wells_inj]))
            faces = np.setdiff1d(faces, self.elements_lv0['boundary_faces'])

            ids_faces_internal = self.rmap_internal_faces(faces)
            self._data['upwind_identificate'][ids_faces_internal] = False
            self._data['upwind_identificate_o'][ids_faces_internal] = False

            v0 = self.elements_lv0['neig_internal_faces'][ids_faces_internal]

            for volumes, i in zip(v0, ids_faces_internal):
                if set_wells_inj & set([volumes[0]]):
                    self._data['upwind_identificate'][i, 0] = True
                    self._data['upwind_identificate_o'][i, 0] = True
                if set_wells_inj & set([volumes[1]]):
                    self._data['upwind_identificate'][i, 1] = True
                    self._data['upwind_identificate_o'][i, 1] = True

    def update_upwind_phases(self):
        '''
            paper Starnoni
        '''
        k0 = 9e-7

        internal_faces = self.elements_lv0['internal_faces']
        v0 = self.elements_lv0['neig_internal_faces']
        saturation = self.data_impress['saturation']
        q_sigma_internal_faces = self.data_impress['flux_faces'][internal_faces]
        q_sigma_internal_faces[np.absolute(q_sigma_internal_faces) < k0] = 0
        qposi = q_sigma_internal_faces >= 0
        qneg = ~qposi
        ak = self.data_impress['area'][internal_faces]*self.data_impress['k_harm'][internal_faces]
        G = ak*(self.g_w_internal_faces - self.g_o_internal_faces)
        gnegi = G <= 0
        gpos = ~gnegi
        test1 = qposi & gnegi
        test2 = qposi & gpos
        test3 = qneg

        self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        self._data['upwind_identificate_o'] = self._data['upwind_identificate'].copy()

        if test1.sum() > 0:
            test1_1 = (q_sigma_internal_faces + G*self.lambda_w_volumes[v0[:, 1]] >= 0) & test1
            test1_2 = (~test1_1) & test1
            sats = saturation[v0[test1_1]]
            sats2 = saturation[v0[test1_2]]

            if test1_1.sum() > 0:
                self._data['upwind_identificate'][test1_1, 0] = True
                self._data['upwind_identificate_o'][test1_1, 0] = True
            if test1_2.sum() > 0:
                self._data['upwind_identificate'][test1_2, 1] = True
                self._data['upwind_identificate_o'][test1_2, 0] = True

        if test2.sum() > 0:
            test2_1 = (q_sigma_internal_faces - G*self.lambda_o_volumes[v0[:,1]] >= 0) & test2
            test2_2 = (~test2_1) & test2

            if test2_1.sum() > 0:
                self._data['upwind_identificate'][test2_1, 0] = True
                self._data['upwind_identificate_o'][test2_1, 0] = True
            if test_2_2.sum() > 0:
                self._data['upwind_identificate'][test2_2, 0] = True
                self._data['upwind_identificate_o'][test2_2, 1] = True

        if test3.sum() > 0:
            q_sigma_2 = q_sigma_internal_faces[test3]
            flux_w_faces = self.data_impress['flux_w_faces'][internal_faces][test3]
            flux_o_faces = self.data_impress['flux_o_faces'][internal_faces][test3]

            pos_sigma = np.full(len(q_sigma_2), True, dtype=bool)
            pos_w = flux_w_faces < 0
            pos_o = flux_o_faces < 0
            and_pos_w = pos_w & pos_sigma
            and_pos_o = pos_o & pos_sigma

            upwind_w = np.full((len(q_sigma_2), 2), False, dtype=bool)
            upwind_o = np.full((len(q_sigma_2), 2), False, dtype=bool)

            upwind_w[and_pos_w, 1] = True
            upwind_o[and_pos_o, 1] = True

            upwind_o[and_pos_w, 1] = True
            upwind_w[and_pos_o, 1] = True

            self._data['upwind_identificate'][test3] = upwind_w
            self._data['upwind_identificate_o'][test3] = upwind_o

            flux_w_internal_faces = (-(self.grad_p_internal_faces*ak*self.lambda_w_internal_faces - self.flux_grav_w_internal_faces))[test3]
            flux_o_internal_faces = (-(self.grad_p_internal_faces*ak*self.lambda_o_internal_faces - self.flux_grav_o_internal_faces))[test3]
            flux_total = flux_w_internal_faces + flux_o_internal_faces

            pos_flux_total = fluxflux_total < 0
            and_pos_sigma = pos_flux_total & pos_sigma

            if ~and_pos_sigma.sum() > 1:
                t1 = ~and_pos_sigma & and_pos_w & ~and_pos_o
                t2 = ~and_pos_sigma & and_pos_o & ~and_pos_w

                upwind_o[t1] = False
                upwind_o[t1, 0] = True

                upwind_w[t2] = False
                upwind[t2, 0] = True

                self._data['upwind_identificate'][test3] = upwind_w
                self._data['upwind_identificate_o'][test3] = upwind_o

        self.upwind_wells()

        verif1 = self._data['upwind_identificate'][:,0] ^ self._data['upwind_identificate'][:,1]
        verif2 = self._data['upwind_identificate_o'][:,0] ^ self._data['upwind_identificate_o'][:,1]
        verif1 = ~verif1
        verif2 = ~verif2
        if verif1.sum() > 0 or verif2.sum() > 0:
            import pdb; pdb.set_trace()

        self.visualize_upwind_vec()

    def update_upwind_phases_new2(self):
        d1 = 0

        internal_faces = self.elements_lv0['internal_faces']
        # flux_o_internal_faces = self.data_impress['flux_o_faces'][internal_faces]
        flux_w_internal_faces = self.data_impress['flux_w_faces'][internal_faces]
        flux_o_internal_faces = self.data_impress['flux_o_faces'][internal_faces]

        # q_sigma_internal_faces = self.data_impress['flux_faces'][internal_faces]
        ak = self.data_impress['area'][internal_faces]*self.data_impress['k_harm'][internal_faces]

        pos_o = flux_o_internal_faces > 0 - d1
        pos_w = flux_w_internal_faces > 0 - d1
        self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        self._data['upwind_identificate_o'] = self._data['upwind_identificate'].copy()
        self._data['upwind_identificate'][pos_w , 0] = True
        self._data['upwind_identificate'][~pos_w , 1] = True
        self._data['upwind_identificate_o'][pos_o, 0] = True
        self._data['upwind_identificate_o'][~pos_o , 1] = True

        self.visualize_upwind_vec()

    def update_upwind_phases_new1(self):
        k0 = 1e-8
        d1 = k0
        d2 = -k0

        internal_faces = self.elements_lv0['internal_faces']
        # flux_o_internal_faces = self.data_impress['flux_o_faces'][internal_faces]
        flux_w_internal_faces = self.data_impress['flux_w_faces'][internal_faces]
        flux_o_internal_faces = self.data_impress['flux_o_faces'][internal_faces]
        # q_sigma_internal_faces = self.flux_sigma_internal_faces
        q_sigma_internal_faces = self.flux_internal_faces
        # q_sigma_internal_faces[np.absolute(q_sigma_internal_faces) < k0] = 0
        # q_sigma_internal_faces = self.data_impress['flux_faces'][internal_faces]
        ak = self.data_impress['area'][internal_faces]*self.data_impress['k_harm'][internal_faces]

        pos_sigma = q_sigma_internal_faces > 0 + d1
        pos_w = flux_w_internal_faces > 0 + d1
        pos_o = flux_o_internal_faces > 0 + d1
        and_pos_w = pos_w & pos_sigma
        and_pos_o = pos_o & pos_sigma

        pos_sigma_2 = q_sigma_internal_faces < 0 + d2
        pos_w_2 = flux_w_internal_faces < 0 + d2
        pos_o_2 = flux_o_internal_faces < 0 + d2
        and_pos_w_2 = pos_w_2 & pos_sigma_2
        and_pos_o_2 = pos_o_2 & pos_sigma_2

        # import pdb; pdb.set_trace()

        ids = np.arange(len(internal_faces))
        # import pdb; pdb.set_trace()
        self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        self._data['upwind_identificate_o'] = np.full((len(internal_faces), 2), False, dtype=bool)

        self._data['upwind_identificate'][and_pos_w, 0] = True
        self._data['upwind_identificate_o'][and_pos_o, 0] = True
        self._data['upwind_identificate'][and_pos_w_2, 1] = True
        self._data['upwind_identificate_o'][and_pos_o_2, 1] = True

        self._data['upwind_identificate_o'][and_pos_w, 0] = True
        self._data['upwind_identificate'][and_pos_o, 0] = True
        self._data['upwind_identificate_o'][and_pos_w_2, 1] = True
        self._data['upwind_identificate'][and_pos_o_2, 1] = True

        verif1 = self._data['upwind_identificate'][:,0] ^ self._data['upwind_identificate'][:,1]
        verif2 = self._data['upwind_identificate_o'][:,0] ^ self._data['upwind_identificate_o'][:,1]
        verif1 = ~verif1
        verif2 = ~verif2

        if verif1.sum() > 0 or verif2.sum() > 0:
            g11 = g1 & verif1
            self._data['upwind_identificate'][:,0]



        self._data['upwind_identificate'][verif1, 0] = True
        self._data['upwind_identificate_o'][verif2, 0] = True

        flux_w_internal_faces = -(self.grad_p_internal_faces*ak*self.lambda_w_internal_faces - self.flux_grav_w_internal_faces)
        flux_o_internal_faces = -(self.grad_p_internal_faces*ak*self.lambda_o_internal_faces - self.flux_grav_o_internal_faces)

        flux_total = flux_w_internal_faces + flux_o_internal_faces
        pos_flux_total = flux_total > 0 + d1
        pos_flux_total_2 = flux_total < 0 + d2
        and_flux_sigma = pos_flux_total & pos_sigma
        and_flux_sigma_2 = pos_sigma_2 & pos_flux_total_2

        t1 = and_pos_w & pos_flux_total
        if not np.allclose(t1, and_pos_w):
            t1 = t1 ^ and_pos_w
            self._data['upwind_identificate_o'][t1] = np.full((t1.sum(), 2), False, dtype=bool)
            self._data['upwind_identificate_o'][t1, 1] = True

        t2 = and_pos_o & pos_flux_total
        if not np.allclose(t1, and_pos_w):
            t2 = t2 ^ and_pos_w
            self._data['upwind_identificate'][t2] = np.full((t2.sum(), 2), False, dtype=bool)
            self._data['upwind_identificate'][t2, 1] = True

        t3 = and_pos_w_2 & pos_flux_total_2
        if not np.allclose(t3, and_pos_w_2):
            t3 = t3 ^ and_pos_w_2
            self._data['upwind_identificate_o'][t3] = np.full((t3.sum(), 2), False, dtype=bool)
            self._data['upwind_identificate_o'][t3, 0] = True

        t4 = and_pos_o_2 & pos_flux_total_2
        if not np.allclose(t4, and_pos_w_2):
            t4 = t4 ^ and_pos_w_2
            self._data['upwind_identificate'][t4] = np.full((t4.sum(), 2), False, dtype=bool)
            self._data['upwind_identificate'][t4, 0] = True

        self.visualize_upwind_vec()

    def update_upwind_phases_new0(self):

        internal_faces = self.elements_lv0['internal_faces']
        flux_w_internal_faces = self.data_impress['flux_w_faces'][internal_faces]
        q_sigma_internal_faces = self.flux_sigma_internal_faces
        ak = self.data_impress['area'][internal_faces]*self.data_impress['k_harm'][internal_faces]

        pos_sigma = q_sigma_internal_faces >= 0
        pos_w = flux_w_internal_faces >= 0
        self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        self._data['upwind_identificate'][pos_w & pos_sigma, 0] = True
        self._data['upwind_identificate'][~(pos_w & pos_sigma), 1] = True
        self._data['upwind_identificate_o'] = self._data['upwind_identificate'].copy()

        qw = self.fw_internal_faces*(q_sigma_internal_faces - self.lambda_o_internal_faces*ak*(self.g_w_internal_faces - self.g_o_internal_faces))
        qo = (1 - self.fw_internal_faces)*(q_sigma_internal_faces + self.lambda_w_internal_faces*ak*(self.g_w_internal_faces - self.g_o_internal_faces))

        q_total = qw + qo
        pos_flux_total = q_total >= 0
        self._data['upwind_identificate_o'][~(pos_flux_total & pos_sigma)] = ~self._data['upwind_identificate_o'][~(pos_flux_total & pos_sigma)]

    def update_transmissibility(self):

        pretransmissibility = self.data_impress['pretransmissibility'].copy()
        internal_faces = self.elements_lv0['internal_faces']
        # b_faces = self.elements_lv0['boundary_faces']
        # t_internal = pretransmissibility[internal_faces]
        # vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        # vols_viz_boundary_faces = self.elements_lv0['neig_boundary_faces'].flatten()
        # fw_vol = self.data_impress['fw_vol']
        # lambda_t = self.data_impress['lambda_t']
        # v0 = vols_viz_internal_faces[:, 0]
        # v1 = vols_viz_internal_faces[:, 1]
        # flux_faces = self.data_impress['flux_faces']
        # flux_internal_faces = flux_faces[internal_faces]



        # total_mobility_internal_faces = np.zeros(len(internal_faces))
        # fw_internal_faces = total_mobility_internal_faces.copy()

        # total_mobility_internal_faces[fluxo_positivo] = lambda_t[v1[fluxo_positivo]]
        # total_mobility_internal_faces[outros] = lambda_t[v0[outros]]

        # total_mobility_internal_faces = self.lambda_t_internal_faces

        # fw_internal_faces[fluxo_positivo] = fw_vol[v1[fluxo_positivo]]
        # fw_internal_faces[outros] = fw_vol[v0[outros]]

        # fw_internal_faces = self.fw_internal_faces

        # gama_faces[internal_faces[fluxo_positivo]] = gama[v1[fluxo_positivo]]
        # gama_faces[internal_faces[outros]] = gama[v0[outros]]

        total_mobility_faces = self.lambda_t_faces
        # fw_faces = self.fw_faces
        # gama_faces[b_faces] = gama[vols_viz_boundary_faces]

        transmissibility = pretransmissibility*total_mobility_faces

        self.data_impress['transmissibility'] = transmissibility
        self.data_impress['lambda_t_faces'] = total_mobility_faces
        self.data_impress['fw_faces'] = self.fw_faces
        # self.data_impress['gama_faces'] = gama_faces.copy()

    def update_loop(self):
        self.loop += 1

    def update_current_biphasic_results(self, simulation_time: float=0.0):

        ws_prod = self.wells['ws_prod']
        fw_vol = self.data_impress['fw_vol']
        water_production = (self.data_impress['flux_volumes'][ws_prod]*fw_vol[ws_prod]).sum()
        oil_production = (self.data_impress['flux_volumes'][ws_prod]).sum() - water_production

        wor = water_production/oil_production

        self.current_biphasic_results = np.array([self.loop, self.delta_t, simulation_time,
            -oil_production, -water_production, self.t, wor, self.vpi, self.contador_vtk])

        self.all_biphasic_results.append(self.current_biphasic_results)

    def save_infos(self):
        self.export_current_biphasic_results()
        self.export_all_biphasic_results()
        self.data_impress.update_variables_to_mesh()
        self.data_impress.export_all_datas_to_npz()
        # self.mesh.core.print(file=self.mesh_name, extension='.h5m', config_input="input_cards/print_settings.yml")

    def export_current_biphasic_results(self):
        np.save(self.name_current_biphasic_results, self.current_biphasic_results)

    def export_all_biphasic_results(self):
        np.save(self.name_all_biphasic_results + str(self.loop) + '.npy', np.array(self.all_biphasic_results))
        self.all_biphasic_results = self.get_empty_current_biphasic_results()

    def run(self, save=False):

        # T, b = super().run()
        # self.update_gama()
        T, b = self.get_T_and_b()
        p = self.solver.direct_solver(T, b)
        self.test_rhs_term(T, b, p)
        self.data_impress['pressure'] = p
        # self.get_total_flux_faces()
        self.get_flux_faces_and_volumes()
        self.test_flux_volumes(self['Tini'], p, self.data_impress['flux_grav_volumes'], -self.data_impress['flux_volumes'])
        self.run_2(save = save)
        # return T, b

    def run_2(self, save=False):
        ######
        ## run for adm_method
        ######

        t0 = time.time()
        self.update_flux_w_and_o_volumes()
        self.test_flux_faces()
        self.update_delta_t()
        self.update_saturation()
        self.update_relative_permeability()
        self.update_mobilities()
        self.update_upwind_phases()
        self.update_transmissibility()
        self.update_t()
        self.update_vpi()
        self.update_loop()
        t1 = time.time()
        dt = t1-t0
        self.update_current_biphasic_results(dt)

        if save:
            self.save_infos()

    def get_T_and_b(self):

        T, b = super().run()
        return T, b

    def print_test(self):
        self.data_impress.update_variables_to_mesh()
        name = 'results/test_volumes_'+str(self.loop)+'.vtk'
        self.mesh.core.mb.write_file(name, [self.meshset_volumes])
        # self.mesh.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")

    def print_test_faces(self):

        self.data_impress.update_variables_to_mesh()
        name = 'results/test_faces_'+str(self.loop)+'.vtk'
        self.mesh.core.mb.write_file(name, [self.meshset_faces])
        # self.mesh.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings1.yml")

    def test_flux_faces(self):

        if data_loaded['_debug'] == False:
            return 0
        fo = self.data_impress['flux_o_faces']
        fw = self.data_impress['flux_w_faces']
        ft = self.data_impress['flux_faces']
        u_normal = self.data_impress['u_normal']
        ft = ft.reshape(len(ft), 1)
        self.data_impress['flux_faces_vec'] = u_normal*ft

        fo2 = fo[abs(fo) > 0]
        fw2 = fw[abs(fw) > 0]

    def visualize_upwind_vec(self):

        if data_loaded['_debug'] == False:
            return 0
        internal_faces = self.elements_lv0['internal_faces']
        u_normal_internal_faces = self.data_impress['u_normal'][internal_faces]
        internal_faces = self.elements_lv0['internal_faces']
        centroid_volumes = self.data_impress['centroid_volumes']
        v0 = self.elements_lv0['neig_internal_faces']
        c0 = centroid_volumes[v0[:,0]]
        c1 = centroid_volumes[v0[:,1]]
        dc = c1-c0
        norm = np.linalg.norm(dc, axis=1).reshape(dc.shape[0], 1)
        dcu = dc/norm
        upwind = np.zeros((len(internal_faces), 3), dtype=int)
        upwind[self._data['upwind_identificate'][:,0]] = dcu[self._data['upwind_identificate'][:,0]].astype(int)
        upwind[self._data['upwind_identificate'][:,1]] = -dcu[self._data['upwind_identificate'][:,1]].astype(int)
        self.data_impress['upwind_w_faces_vec'][internal_faces] = upwind
        upwind = np.zeros((len(internal_faces), 3), dtype=int)
        upwind[self._data['upwind_identificate_o'][:,0]] = dcu[self._data['upwind_identificate_o'][:,0]].astype(int)
        upwind[self._data['upwind_identificate_o'][:,1]] = -dcu[self._data['upwind_identificate_o'][:,1]].astype(int)
        self.data_impress['upwind_o_faces_vec'][internal_faces] = upwind

        faces = self.elements_lv0['faces']
        ctfaces = self.data_impress['centroid_faces']

        f1 = get_box(ctfaces, np.array([[0.0, 0.0, 6.0 - 0.1], [9.0, 9.0, 6.0+0.1]]))

        flux = self.data_impress['flux_faces'][f1]
        upw = self.data_impress['upwind_w_faces_vec'][f1]

        # self.print_test_faces()
        # self.print_test()

        # import pdb; pdb.set_trace()
