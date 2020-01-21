from ..pressure_solver.fine_scale_tpfa import FineScaleTpfaPressureSolver
from ..directories import data_loaded
from ..utils import relative_permeability
from ..solvers.solvers_scipy.solver_sp import SolverSp
import os
from .. import directories as direc
import numpy as np
import scipy.sparse as sp
import time

class BiphasicTpfa(FineScaleTpfaPressureSolver):

        def __init__(self, M, data_impress, elements_lv0, wells, data_name: str='BiphasicTpfa.npz'):
        load = data_loaded['load_biphasic_data']
        super().__init__(data_impress, elements_lv0, wells, data_name=data_name, load=load)
        self.biphasic_data = data_loaded['biphasic_data']
        self.relative_permeability = getattr(relative_permeability, self.biphasic_data['relative_permeability'])
        self.relative_permeability = self.relative_permeability()
        self.V_total = (data_impress['volume']*data_impress['poro']).sum()
        self.max_contador_vtk = len(self.biphasic_data['vpis_para_gravar_vtk'])
        self.delta_sat_max = 0.4
        self.lim_flux_w = 9e-8
        self.name_current_biphasic_results = os.path.join(direc.flying, 'current_biphasic_results.npy')
        self.name_all_biphasic_results = os.path.join(direc.flying, 'all_biphasic_results_')
        self.mesh_name = os.path.join(direc.flying, 'biphasic_')
        self.all_biphasic_results = self.get_empty_current_biphasic_results()
        self.mesh = M
        self.solver = SolverSp()

        if not load:
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.contador_vtk = 0
            self.update_relative_permeability()
            self.update_mobilities()
            self.update_transmissibility_ini()
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
        lambda_w = self.data_impress['krw']/self.biphasic_data['mi_w']
        lambda_o = self.data_impress['kro']/self.biphasic_data['mi_o']
        lambda_t = lambda_w + lambda_o
        fw_vol = lambda_w/lambda_t

        self.data_impress['lambda_w'] = lambda_w
        self.data_impress['lambda_o'] = lambda_o
        self.data_impress['lambda_t'] = lambda_t
        self.data_impress['fw_vol'] = fw_vol

    def update_transmissibility_ini(self):

        pretransmissibility = self.data_impress['pretransmissibility'].copy()
        internal_faces = self.elements_lv0['internal_faces']
        b_faces = self.elements_lv0['boundary_faces']
        t_internal = pretransmissibility[internal_faces]
        vols_viz_internal_faces = self.elements_lv0['neig_internal_faces']
        vols_viz_boundary_faces = self.elements_lv0['neig_boundary_faces']
        fw_vol = self.data_impress['fw_vol']
        lambda_t = self.data_impress['lambda_t']

        ws_inj = self.wells['ws_inj']

        v0 = vols_viz_internal_faces[:, 0]
        v1 = vols_viz_internal_faces[:, 1]
        ids = np.arange(len(v0))

        idv0 = np.array([], dtype=np.int64)
        idv1 = idv0.copy()
        for w in ws_inj:
            vv0 = ids[v0 == w]
            vv1 = ids[v1 == w]
            idv0 = np.append(idv0, vv0)
            idv1 = np.append(idv1, vv1)

        idv0 = np.array(idv0).flatten()
        idv1 = np.array(idv1).flatten()
        idsv = np.union1d(idv0, idv1)
        ids_fora = np.setdiff1d(ids, idsv)

        total_mobility_internal_faces = np.zeros(len(internal_faces))
        fw_internal_faces = total_mobility_internal_faces.copy()

        total_mobility_internal_faces[idv0] = lambda_t[v0[idv0]]
        total_mobility_internal_faces[idv1] = lambda_t[v1[idv1]]
        total_mobility_internal_faces[ids_fora] = (lambda_t[v0[ids_fora]] + lambda_t[v1[ids_fora]])/2

        fw_internal_faces[idv0] = fw_vol[v0[idv0]]
        fw_internal_faces[idv1] = fw_vol[v1[idv1]]
        fw_internal_faces[ids_fora] = (fw_vol[v0[ids_fora]] + fw_vol[v1[ids_fora]])/2

        total_mobility_faces = np.zeros(len(pretransmissibility))
        fw_faces = total_mobility_faces.copy()

        total_mobility_faces[internal_faces] = total_mobility_internal_faces
        total_mobility_faces[b_faces] = lambda_t[vols_viz_boundary_faces]

        fw_faces[internal_faces] = fw_internal_faces
        fw_faces[b_faces] = fw_vol[vols_viz_boundary_faces]

        transmissibility = pretransmissibility*total_mobility_faces

        self.data_impress['transmissibility'] = transmissibility.copy()
        self.data_impress['lambda_t_faces'] = total_mobility_faces.copy()
        self.data_impress['fw_faces'] = fw_faces.copy()

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
        total_flux_faces = self.data_impress['flux_faces']
        fw_faces = self.data_impress['fw_faces']
        fw_vol = self.data_impress['fw_vol']
        flux_volumes = self.data_impress['flux_volumes']
        ws_prod = self.wells['ws_prod']
        ws_inj = self.wells['ws_inj']

        flux_w_faces = fw_faces*total_flux_faces
        flux_w_internal_faces = flux_w_faces[internal_faces]

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*fw_vol[ws_prod]
        flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*fw_vol[ws_inj]

        flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - fw_vol[ws_prod])

        self.data_impress['flux_w_faces'] = flux_w_faces
        self.data_impress['flux_w_volumes'] = flux_w_volumes
        self.data_impress['flux_o_volumes'] = flux_o_volumes

    def update_delta_t(self):
        ###
        ## de acordo com o fluxo nos volumes
        ###

        flux_volumes = np.absolute(self.data_impress['flux_volumes'])
        phis = self.data_impress['poro']
        volume = self.data_impress['volume']
        self.delta_t = (self.biphasic_data['cfl']*(volume*phis)/flux_volumes).min()

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
        self.delta_t = (self.biphasic_data['cfl']*(dists_int/(vel_internal_faces*dfds))).min()

    def update_saturation(self):
        self.data_impress['saturation_last'] = self.data_impress['saturation'].copy()
        verif = -1
        while verif != 0:
            verif = self.update_sat()
            if verif == 1:
                self.reduce_delta_t()

    def update_sat(self):

        saturations0 = self.data_impress['saturation'].copy()
        saturations = saturations0.copy()
        ids = np.arange(len(saturations))

        fw_volumes = -self.data_impress['flux_w_volumes']
        volumes = self.data_impress['volume']
        phis = self.data_impress['poro']

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
        test = ids[(saturations < 0) & (saturations > 1)]
        if len(test) > 0:
            import pdb; pdb.set_trace()
            raise ValueError(f'valor errado da saturacao {saturations[test]}')
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


        ids_var = ids[fw_volumes > self.lim_flux_w]

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
        test = ids2[delta_sat > self.delta_sat_max]
        if len(test) > 0:
            return 1
        ##############

        saturations[ids_var] = sats

        if np.allclose(saturations, saturations0):
            import pdb; pdb.set_trace()

        self.data_impress['saturation'] = saturations

        return 0

    def reduce_delta_t(self):
        d = self.delta_t
        self.delta_t *= 1/2
        print(f'\nreducing delta_t: {d} -> {self.delta_t} \n')

    def update_t(self):
        self.t += self.delta_t

    def update_vpi(self):

        flux_total_inj = np.absolute(self.data_impress['flux_volumes'][self.wells['ws_inj']])
        self.vpi += (flux_total_inj.sum()*self.delta_t)/self.V_total

    def update_transmissibility(self):

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

        ids = np.arange(len(internal_faces))

        fluxo_positivo = ids[flux_internal_faces <= 0]
        outros = np.setdiff1d(ids, fluxo_positivo)

        total_mobility_internal_faces = np.zeros(len(internal_faces))
        fw_internal_faces = total_mobility_internal_faces.copy()

        total_mobility_internal_faces[fluxo_positivo] = lambda_t[v1[fluxo_positivo]]
        total_mobility_internal_faces[outros] = lambda_t[v0[outros]]

        fw_internal_faces[fluxo_positivo] = fw_vol[v1[fluxo_positivo]]
        fw_internal_faces[outros] = fw_vol[v0[outros]]

        total_mobility_faces = np.zeros(len(pretransmissibility))
        fw_faces = total_mobility_faces.copy()

        total_mobility_faces[internal_faces] = total_mobility_internal_faces
        total_mobility_faces[b_faces] = lambda_t[vols_viz_boundary_faces]

        fw_faces[internal_faces] = fw_internal_faces
        fw_faces[b_faces] = fw_vol[vols_viz_boundary_faces]

        transmissibility = pretransmissibility*total_mobility_faces

        self.data_impress['transmissibility'] = transmissibility.copy()
        self.data_impress['lambda_t_faces'] = total_mobility_faces.copy()
        self.data_impress['fw_faces'] = fw_faces.copy()

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

        T, b = super().run()
        p = self.solver.direct_solver(T, b)
        self.data_impress['pressure'] = p
        self.get_flux_faces_and_volumes()
        self.run_2(save = save)

    def run_2(self, save=False):
        ######
        ## run for adm_method
        ######

        t0 = time.time()
        self.update_flux_w_and_o_volumes()
        self.update_delta_t()
        self.update_saturation()
        self.update_t()
        self.update_vpi()
        self.update_relative_permeability()
        self.update_mobilities()
        self.update_transmissibility()
        self.update_loop()
        t1 = time.time()
        dt = t1-t0
        self.update_current_biphasic_results(dt)

        if save:
            self.save_infos()
