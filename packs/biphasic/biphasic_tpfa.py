from ..pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from ..pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from ..pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
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
import pdb

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
        self.delta_sat_max = 0.4
        self.lim_flux_w = 0
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
        self.pare1 = False

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
        self._data['upwind_identificate'], self._data['upwind_identificate_o'] = self.upwind_wells(
            self.wells['ws_inj'],
            self.mesh.volumes.bridge_adjacencies(self.mesh.volumes.all[:], 3, 2),
            b_faces,
            self.elements_lv0['remaped_internal_faces'],
            vols_viz_internal_faces,
            self._data['upwind_identificate'],
            self._data['upwind_identificate_o']
        )
        
        self.test_upwind_dup()

        self.visualize_upwind_vec()

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

        transmissibility = pretransmissibility*total_mobility_faces

        self.data_impress['transmissibility'] = transmissibility.copy()
        self.data_impress['lambda_t_faces'] = total_mobility_faces.copy()
        self.data_impress['fw_faces'] = fw_faces.copy()
        # self.data_impress['gama_faces'] = gama_faces.copy()

    def upwind_wells(self, injector_volumes, faces_adj_volumes, boundary_faces, map_internal_faces, volumes_adj_internal_faces, upwind_w, upwind_o):

        wells_inj = injector_volumes
        if len(wells_inj) > 0:
            set_wells_inj = set(wells_inj)
            faces = np.unique(np.concatenate(faces_adj_volumes[wells_inj]))
            faces = np.setdiff1d(faces, boundary_faces)

            ids_faces_internal = map_internal_faces[faces]
            upwind_w[ids_faces_internal] = False
            upwind_o[ids_faces_internal] = False

            v0 = volumes_adj_internal_faces[ids_faces_internal]

            for volumes, i in zip(v0, ids_faces_internal):
                if set_wells_inj & set([volumes[0]]):
                    upwind_w[i, 0] = True
                    upwind_o[i, 0] = True
                elif set_wells_inj & set([volumes[1]]):
                    upwind_w[i, 1] = True
                    upwind_o[i, 1] = True

        return upwind_w, upwind_o

    def get_empty_current_biphasic_results(self):

        return [np.array(['loop', 'delta_t [s]', 'simulation_time [s]',
            'oil_production [m3/s]', 'water_production [m3/s]', 't [s]', 'wor', 'vpi', 'contador_vtk'])]

    def load_infos(self):
        self.current_biphasic_results = np.load(self.name_current_biphasic_results)
        self.loop = int(self.current_biphasic_results[0])
        self.t = self.current_biphasic_results[5]
        self.vpi = self.current_biphasic_results[7]
        self.contador_vtk = self.current_biphasic_results[8]

    def update_flux_w_and_o_volumes(self):
        u_normal = self.data_impress['u_normal']
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

        flux_internal_faces = self.data_impress['flux_faces'][internal_faces]

        ps0 = x[v0[:, 0]]
        ps1 = x[v0[:, 1]]
        # ps0 = x[v0[:, 0]]
        # ps1 = x[v0[:, 1]]

        # flux_w_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_w_internal_faces/dh_internal_faces - self._data['grav_source_term_water_faces'][internal_faces])
        flux_w_internal_faces = flux_internal_faces*fw_internal_faces
        # self.data_impress['flux_press_w_faces_vec'][internal_faces] = (-((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_w_internal_faces/dh_internal_faces)).reshape(len(internal_faces), 1)*u_normal_internal_faces
        self.data_impress['flux_press_w_faces_vec'][internal_faces] = flux_w_internal_faces.reshape((ni, 1))*u_normal_internal_faces

        # flux_o_internal_faces = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_faces'][internal_faces]))
        flux_o_internal_faces =  flux_internal_faces*(1 - fw_internal_faces)
        flux_int_faces2  = flux_o_internal_faces + flux_w_internal_faces
        test = np.allclose(flux_int_faces2, flux_internal_faces)

        self.data_impress['flux_press_o_faces_vec'][internal_faces] = self.flux_press_o_internal_faces.reshape(len(internal_faces), 1)*u_normal_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        flux_w_volumes2 = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_w_volumes=np.bincount(lines,weights=data)




        flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*fw_vol[ws_prod]
        flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*fw_vol[ws_inj]


        # flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        # lambda_o_internal_faces = self.data_impress['lambda_t'][v0[self._data['upwind_identificate']]] - self.data_impress['lambda_w'][v0[self._data['upwind_identificate']]]
        # flux_o_internal_faces2 = -((ps1 - ps0)*areas_internal_faces*k_harm_internal_faces*lambda_o_internal_faces/dh_internal_faces - (self.data_impress['flux_grav_faces'][internal_faces] - self._data['grav_source_term_water_internal_faces']))

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        # cols = np.repeat(0, len(lines))
        data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        # flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_o_volumes=np.bincount(lines,weights=data)

        flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - fw_vol[ws_prod])
        flux_o_volumes[ws_inj] -= flux_volumes[ws_inj]*(1 - fw_vol[ws_inj])

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
        non_null_flux=flux_volumes>0
        phis = self.data_impress['poro']
        volume = self.data_impress['volume']
        delta_t = ((self.biphasic_data['cfl']*(volume*phis))[non_null_flux]/flux_volumes[non_null_flux]).min()
        return delta_t

    def update_delta_t(self):
        ###
        ## de acordo com o fluxo de agua nos volumes
        ###
        deltas_t=[]
        deltas_t.append(self.update_delta_t_for_delta_sat_max())

        deltas_t.append(self.update_delta_t_dep0())

        # deltas_t.append(self.update_delta_t_new())


        self.delta_t = min(deltas_t)

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
        # vel_internal_faces = np.linalg.norm(velocity_faces[internal_faces], axis=1)
        vel_internal_faces = velocity_faces[internal_faces].max(axis=1)

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
        saturations0 = self.data_impress['saturation'].copy()
        # import pdb; pdb.set_trace()

        saturations = saturations0.copy()
        ids = np.arange(len(saturations))

        fw_volumes = -self.data_impress['flux_w_volumes'].copy()
        volumes = self.data_impress['volume'].copy()
        phis = self.data_impress['poro'].copy()

        test = ids[(saturations < 0) | (saturations > 1)]
        if len(test) > 0:

            raise ValueError(f'valor errado da saturacao {saturations[test]}')
        del test

        ids_var = ids[np.absolute(fw_volumes) > self.lim_flux_w]

        if len(ids_var) == 0:
            import pdb; pdb.set_trace()

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

        min_sat = saturations.min()
        max_sat = saturations.max()
        # saturations[saturations<0.2]=0.2
        # saturations[saturations>0.8]=0.8

        '''
        if min_sat < self.biphasic_data['Swc']-0.1 or max_sat > 1-self.biphasic_data['Sor']+0.1:
            return 1
            # raise ValueError(f'\nprint max_sat: {max_sat} ; min_sat: {min_sat}\n')
        '''
        self.data_impress['saturation'] = saturations

        return 0

    def reduce_delta_t(self):
        d = self.delta_t
        self.delta_t *= 1/2
        print(f'\nreducing delta_t: {d} -> {self.delta_t} \n')
        # self.pare1 = True

    def update_t(self):
        self.t += self.delta_t

    def update_vpi(self):

        flux_total_inj = np.absolute(self.data_impress['flux_volumes'][self.wells['ws_inj']])
        self.vpi += (flux_total_inj.sum()*self.delta_t)/self.V_total

    def update_transmissibility(self):
        t0=time.time()
        pretransmissibility = self.data_impress['pretransmissibility'].copy()
        internal_faces = self.elements_lv0['internal_faces']

        b_faces = self.elements_lv0['boundary_faces']
        t_internal = pretransmissibility[internal_faces]
        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        vols_viz_boundary_faces = self.elements_lv0['neig_boundary_faces'].flatten()
        fw_vol = self.data_impress['fw_vol']
        lambda_t = self.data_impress['lambda_t']
        v0 = vols_viz_internal_faces[:, 0]
        v1 = vols_viz_internal_faces[:, 1]
        flux_faces = self.data_impress['flux_faces']
        flux_internal_faces = flux_faces[internal_faces]
        # gama = self.data_impress['gama']
        # gama_faces = np.zeros(len(self.data_impress['gama_faces']))

        # ids = np.arange(len(internal_faces))
        #
        pos = flux_internal_faces >= 0
        # outros = np.setdiff1d(ids, fluxo_positivo)


        # self._data['upwind_identificate'] = np.full((len(internal_faces), 2), False, dtype=bool)
        # self._data['upwind_identificate'][pos, 0] = np.full(pos.sum(), True, dtype=bool)
        # pos = ~pos
        # self._data['upwind_identificate'][pos, 1] = np.full(pos.sum(), True, dtype=bool)
        self._data['upwind_identificate']=np.vstack([pos,~pos]).T

        # total_mobility_internal_faces = np.zeros(len(internal_faces))
        # fw_internal_faces = total_mobility_internal_faces.copy()

        # total_mobility_internal_faces[fluxo_positivo] = lambda_t[v1[fluxo_positivo]]
        # total_mobility_internal_faces[outros] = lambda_t[v0[outros]]
        viz_up=vols_viz_internal_faces[self._data['upwind_identificate']]
        total_mobility_internal_faces = lambda_t[viz_up]

        # fw_internal_faces[fluxo_positivo] = fw_vol[v1[fluxo_positivo]]
        # fw_internal_faces[outros] = fw_vol[v0[outros]]

        fw_internal_faces = fw_vol[viz_up]

        # gama_faces[internal_faces[fluxo_positivo]] = gama[v1[fluxo_positivo]]
        # gama_faces[internal_faces[outros]] = gama[v0[outros]]

        total_mobility_faces = np.zeros(len(pretransmissibility))
        fw_faces = np.zeros(len(pretransmissibility))

        total_mobility_faces[internal_faces] = total_mobility_internal_faces
        total_mobility_faces[b_faces] = lambda_t[vols_viz_boundary_faces]

        fw_faces[internal_faces] = fw_internal_faces
        fw_faces[b_faces] = fw_vol[vols_viz_boundary_faces]
        # gama_faces[b_faces] = gama[vols_viz_boundary_faces]

        transmissibility = pretransmissibility*total_mobility_faces

        self.data_impress['transmissibility'] = transmissibility
        self.data_impress['lambda_t_faces'] = total_mobility_faces
        self.data_impress['fw_faces'] = fw_faces
        # self.data_impress['gama_faces'] = gama_faces.copy()

    def update_loop(self):
        self.loop += 1

    def update_current_biphasic_results(self, simulation_time: float=0.0):

        ws_prod = self.wells['ws_prod']
        fw_vol = self.data_impress['fw_vol']

        # water_production = (self.data_impress['flux_volumes'][ws_prod]*fw_vol[ws_prod]).sum()
        # oil_production = (self.data_impress['flux_volumes'][ws_prod]).sum() - water_production
        # import pdb; pdb.set_trace()
        faces=np.hstack(self.elements_lv0['volumes_face_faces'][ws_prod])
        water_production=-abs((self.data_impress['fw_faces'][faces]*self.data_impress['flux_w_faces'][faces])).sum()
        oil_production = (self.data_impress['flux_volumes'][ws_prod]).sum() - water_production
        if oil_production!=0:
            self.data_impress['saturation'][ws_prod]=(self.relative_permeability.Swc*oil_production+(1-self.relative_permeability.Sor)*water_production)/(water_production+oil_production)
            # oil_production=self.data_impress['fo_faces'][np.hstack(self.elements_lv0['volumes_face_volumes'][ws_prod])].sum()
            wor = water_production/oil_production
        else:
            wor = 0

        self.wor=wor
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
        t0=time.time()
        # T, b = super().run()
        # self.update_gama()
        T, b = self.get_T_and_b()
        p = self.solver.direct_solver(T, b)
        self.test_rhs_term(T, b, p)
        self.data_impress['pressure'] = p
        # self.get_total_flux_faces()
        self.get_flux_faces_and_volumes()
        # self.test_flux_volumes(self['Tini'], p, self.data_impress['flux_grav_volumes'], -self.data_impress['flux_volumes'])
        self.run_2(save = save)
        t1=time.time()
        dt=t1-t0
        self.t_comp=dt

    def run_2(self, save=False):
        ######
        ## run for adm_method
        ######        
        t0=time.time()
        self.update_flux_w_and_o_volumes()
        self.test_flux_faces()
        self.update_delta_t()
        # self.update_saturation()
        self.loop_for_saturation()
        self.update_relative_permeability()
        self.update_mobilities()
        # self.data_impress['flux_w_volumes'][:] = self.update_flux_w_volumes_wells
        # self.data_impress['flux_o_volumes'][:] = self.update_flux_o_volumes_wells
        self.update_upwind_phases()
        self.update_transmissibility()
        self.update_t()
        self.update_vpi()
        self.update_loop()
        t1=time.time()
        dt=t1-t0
        self.update_current_biphasic_results(dt)
        if save:
            self.save_infos()

    def get_T_and_b(self):

        T, b = super().run()
        return T, b

    def print_test(self):
        self.data_impress.update_variables_to_mesh()
        name = 'results/test_volumes_'+str(self.loop)+'.vtk'
        # self.mesh.core.print(file=name, extension='.vtk', config_input="input_cards/print_settings0.yml")
        self.mesh.core.mb.write_file(name, [self.meshset_volumes])

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

        try:
            fo2 = fo[abs(fo) > 0]
            fw2 = fw[abs(fw) > 0]
        except:
            pdb.set_trace()

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

    def test_upwind_dup(self):

        if data_loaded['_debug'] == False:
            return None

        verif1 = self._data['upwind_identificate'][:,0] ^ self._data['upwind_identificate'][:,1]
        verif2 = self._data['upwind_identificate_o'][:,0] ^ self._data['upwind_identificate_o'][:,1]
        verif1 = ~verif1
        verif2 = ~verif2
        if verif1.sum() > 0 or verif2.sum() > 0:
            import pdb; pdb.set_trace()

    def loop_for_saturation(self, n_loops=4):

        self.update_saturation()
        sat1w = self.data_impress['saturation'][self.wells['all_wells']]
        s1 = self.data_impress['saturation'][self.wells['all_wells']]
        print('\ns1: ', s1, '\n')

        delta = 0.001
        verif = True
        loop_int_max = 10
        loop_int = 0

        # for i in range(n_loops):
        #     self.data_impress['saturation'][self.wells['all_wells']] = self.update_saturation_wells
        #     s2 = self.data_impress['saturation'][self.wells['all_wells']]
        #     dsmax = np.absolute(s2 - s1).max()
        #     if dsmax > delta:
        #         print('\nentrou\n')
        #         print('dsmax: ', dsmax)
        #         print()
        #         s1 = s2.copy()
        #         pass
        # print('\nsaiu\n')

        if self.loop >= 22:
            self.pare1 = True
            pdb.set_trace()

        while verif and loop_int <= loop_int_max:
            self.data_impress['saturation'][self.wells['all_wells']] = self.update_saturation_wells
            s2 = self.data_impress['saturation'][self.wells['all_wells']]
            dsmax = np.absolute(s2 - s1).max()
            if dsmax > delta:
                print('\nentrou\n')
                print('dsmax: ', dsmax)
                print()
                s1 = s2.copy()
            else:
                print('\nsaiu\n')
                verif = False
            loop_int += 1


        sat2w = self.data_impress['saturation'][self.wells['all_wells']]

        # if self.loop >= 21:
        #     print(sat1w)
        #     print(sat2w)
        #     pdb.set_trace()

        # if self.loop >= 20:
        #     self.print_test()
        #     pdb.set_trace()
