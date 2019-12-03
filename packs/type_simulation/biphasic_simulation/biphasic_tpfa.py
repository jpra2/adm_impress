from ..monophasic_simulation.monophasic_tpfa import monophasicTpfa
from ... import directories as direc
from ...utils import relative_permeability
from ...preprocess.preprocess1 import set_saturation_regions
import numpy as np
import scipy.sparse as sp
from ...convert_unit.conversion import convert1
import time
import os


class biphasicTpfa(monophasicTpfa):

    def __init__(self, M, data_name: str='biphasicTpfa.npz', load=False) -> None:
        super().__init__(M, data_name)
        name_relative_permeability = direc.data_loaded['biphasic_data']['relative_permeability']
        self.relative_permeability = getattr(relative_permeability, name_relative_permeability)
        self.mi_w = direc.data_loaded['biphasic_data']['mi_w']
        self.mi_o = direc.data_loaded['biphasic_data']['mi_o']
        self.cfl = direc.data_loaded['biphasic_data']['cfl']
        self.V_total = (M.data['volume']*M.data['poro']).sum()
        self.Sor = float(direc.data_loaded['biphasic_data']['Sor'])
        self.Swc = float(direc.data_loaded['biphasic_data']['Swc'])
        self.delta_sat_max = 0.6
        self.lim_flux_w = 1e-9
        self.name_hist = os.path.join(direc.flying, 'hist.npy')
        self.name_hist2 = os.path.join(direc.flying, 'hist2_')

        if not load:

            set_saturation_regions(M)
            convert1(M)
            self.update_gama()
            self.update_relative_permeability()
            self.update_mobilities()
            self.update_transmissibility_ini()
            M.data.update_variables_to_mesh()
            M.data.export_variables_to_npz()
            self.loop = 0
            self.vpi = 0.0
            self.t = 0.0
            self.hist2 = self.get_empty_hist()
        else:
            self.load_infos()

    def update_gama(self):

        if self.gravity:

            M = self.mesh

            gama_w = direc.data_loaded['biphasic_data']['gama_w']
            gama_o = direc.data_loaded['biphasic_data']['gama_o']

            gama_w = np.repeat(gama_w, self.n_volumes)
            gama_o = np.repeat(gama_o, self.n_volumes)

            saturations = M.data['saturation']

            gama = gama_w*saturations + gama_o*(1-saturations)
            M.data['gama'] = gama

    def update_flux_w_and_o_volumes(self) -> None:

        M = self.mesh

        vols_viz_internal_faces = M.data.elements_lv0[direc.entities_lv0_0[2]]
        v0 = vols_viz_internal_faces
        internal_faces = M.data.elements_lv0[direc.entities_lv0_0[0]]
        total_flux_faces = M.data['flux_faces']
        fw_faces = M.data['fw_faces']
        fw_vol = M.data['fw_vol']
        flux_volumes = M.data['flux_volumes']
        flux_w_faces = fw_faces*total_flux_faces
        flux_w_internal_faces = flux_w_faces[internal_faces]

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_w_internal_faces, -flux_w_internal_faces]).flatten()
        flux_w_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()

        ws_prod = M.contours.datas['ws_prod']
        ws_inj = M.contours.datas['ws_inj']
        flux_w_volumes[ws_prod] -= flux_volumes[ws_prod]*fw_vol[ws_prod]
        flux_w_volumes[ws_inj] -= flux_volumes[ws_inj]*fw_vol[ws_inj]

        flux_o_internal_faces = total_flux_faces[internal_faces] - flux_w_internal_faces

        lines = np.array([v0[:, 0], v0[:, 1]]).flatten()
        cols = np.repeat(0, len(lines))
        data = np.array([flux_o_internal_faces, -flux_o_internal_faces]).flatten()
        flux_o_volumes = sp.csc_matrix((data, (lines, cols)), shape=(self.n_volumes, 1)).toarray().flatten()
        flux_o_volumes[ws_prod] -= flux_volumes[ws_prod]*(1 - fw_vol[ws_prod])

        M.data['flux_w_faces'] = flux_w_faces
        M.data['flux_w_volumes'] = flux_w_volumes
        M.data['flux_o_volumes'] = flux_o_volumes

    def update_mobilities(self):
        M = self.mesh
        n = self.n_volumes
        lambda_w = M.data['krw']/self.mi_w
        lambda_o = M.data['kro']/self.mi_o
        lambda_t = lambda_w + lambda_o
        fw_vol = lambda_w/lambda_t

        M.data['lambda_w'] = lambda_w
        M.data['lambda_o'] = lambda_o
        M.data['lambda_t'] = lambda_t
        M.data['fw_vol'] = fw_vol

    def update_relative_permeability(self):

        M = self.mesh
        krw, kro = self.relative_permeability(M.data['saturation'])
        M.data['krw'] = krw
        M.data['kro'] = kro

    def update_transmissibility_ini(self):
        M = self.mesh

        pretransmissibility = M.data['pretransmissibility'].copy()
        internal_faces = M.data.elements_lv0['internal_faces']
        b_faces = M.data.elements_lv0['boundary_faces']
        t_internal = pretransmissibility[internal_faces]
        vols_viz_internal_faces = M.data.elements_lv0['vols_viz_internal_faces']
        vols_viz_boundary_faces = M.data.elements_lv0['vols_viz_boundary_faces'].flatten()
        fw_vol = M.data['fw_vol']
        lambda_t = M.data['lambda_t']

        ws_inj = M.contours.datas['ws_inj']

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

        M.data['transmissibility'] = transmissibility.copy()
        M.data['lambda_t_faces'] = total_mobility_faces.copy()
        M.data['fw_faces'] = fw_faces.copy()

    def update_transmissibility(self):
        M = self.mesh

        pretransmissibility = M.data['pretransmissibility'].copy()
        internal_faces = M.data.elements_lv0['internal_faces']
        b_faces = M.data.elements_lv0['boundary_faces']
        t_internal = pretransmissibility[internal_faces]
        vols_viz_internal_faces = M.data.elements_lv0['vols_viz_internal_faces']
        vols_viz_boundary_faces = M.data.elements_lv0['vols_viz_boundary_faces'].flatten()
        fw_vol = M.data['fw_vol']
        lambda_t = M.data['lambda_t']
        v0 = vols_viz_internal_faces[:, 0]
        v1 = vols_viz_internal_faces[:, 1]

        flux_faces = M.data['flux_faces']
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

        M.data['transmissibility'] = transmissibility.copy()
        M.data['lambda_t_faces'] = total_mobility_faces.copy()
        M.data['fw_faces'] = fw_faces.copy()

    def update_delta_t(self):
        M = self.mesh

        flux_volumes = np.absolute(M.data['flux_volumes'])
        phis = M.data['poro']
        volume = M.data['volume']
        self.delta_t = (self.cfl*(volume*phis)/flux_volumes).min()

    def reduce_delta_t(self):
        d = self.delta_t
        self.delta_t *= 1/2
        print(f'\nreducing delta_t: {d} -> {self.delta_t} \n')

    def update_t(self):
        self.t += self.delta_t

    def update_vpi(self):
        M = self.mesh
        flux_total_inj = np.absolute(M.data['flux_volumes'][M.contours.datas['ws_inj']])
        self.vpi += (flux_total_inj.sum()*self.delta_t)/self.V_total

    def update_loop(self):
        self.loop += 1

    def update_hist(self, simulation_time: float=0.0):
        M = self.mesh
        ws_prod = M.contours.datas['ws_prod']
        fw_vol = M.data['fw_vol']
        water_production = (M.data['flux_volumes'][ws_prod]*fw_vol[ws_prod]).sum()
        oil_production = (M.data['flux_volumes'][ws_prod]).sum() - water_production

        wor = water_production/oil_production

        self.hist = np.array([self.loop, self.delta_t, simulation_time,
            -oil_production, -water_production, self.t, wor, self.vpi])

        self.hist2.append(self.hist)

    def update_sat(self):
        M = self.mesh

        saturations0 = M.data['saturation'].copy()
        saturations = saturations0.copy()
        ids = np.arange(len(saturations))

        fw_volumes = -M.data['flux_w_volumes']
        volumes = M.data['volume']
        phis = M.data['poro']

        ###########################
        ## teste
        test = ids[(saturations < 0) & (saturations > 1)]
        if len(test) > 0:
            import pdb; pdb.set_trace()
            raise ValueError(f'valor errado da saturacao {saturations[test]}')
        ###########################

        ###########################
        ## teste
        test = ids[fw_volumes < -self.lim_flux_w]
        if len(test) > 0:
            import pdb; pdb.set_trace()
            raise ValueError(f'valor errado da saturacao {saturations[test]}')
        ###########################

        ids_var = ids[(fw_volumes > self.lim_flux_w)]

        fw_volumes = fw_volumes[ids_var]
        volumes = volumes[ids_var]
        phis = phis[ids_var]
        sats = saturations[ids_var]

        sats += (fw_volumes*self.delta_t)/(phis*volumes)

        delta_sat = sats - saturations[ids_var]
        ids2 = np.arange(len(delta_sat))

        #############
        ## teste
        test = ids2[delta_sat > self.delta_sat_max]
        if len(test) > 0:
            return 1
        ##############

        saturations[ids_var] = sats

        M.data['saturation'] = saturations

        return 0

    def update_saturation(self):
        self.mesh.data['saturation_last'] = self.mesh.data['saturation'].copy()
        verif = -1
        while verif != 0:
            verif = self.update_sat()
            if verif == 1:
                self.reduce_delta_t()

    def get_empty_hist(self):

        # self.dtype = np.dtype([
        # ('loop', int, (1,)),
        # ('delta_t (s)', float, (1,)),
        # ('simulation_time (s)', float, (1,)),
        # ('oil_production (m3/s)', float, (1,)),
        # ('water_production (m3/s)', float, (1,)),
        # ('t (s)', float, (1,)),
        # ('wor', float, (1,))])

        # return np.array([np.array(['loop', 'delta_t [s]', 'simulation_time [s]',
        #     'oil_production [m3/s]', 'water_production [m3/s]', 't [s]', 'wor'])])

        return [np.array(['loop', 'delta_t [s]', 'simulation_time [s]',
            'oil_production [m3/s]', 'water_production [m3/s]', 't [s]', 'wor', 'vpi'])]

    def export_hist(self):
        np.save(self.name_hist, self.hist)

    def export_hist2(self):
        np.save(self.name_hist2 + str(self.loop) + '.npy', np.array(self.hist2))
        self.hist2 = self.get_empty_hist()

    def load_infos(self):
        self.mesh.data.load_from_npz()
        self.mesh.data.update_variables_to_mesh()
        self.hist = np.load(self.name_hist)
        self.loop = int(self.hist[0])
        self.t = self.hist[5]
        self.vpi = self.hist[7]
        self.hist2 = self.get_empty_hist()

    def run(self, save=True):
        t0 = time.time()
        super().run()
        self.update_flux_w_and_o_volumes()
        self.update_delta_t()
        self.update_saturation()
        self.update_t()
        self.update_vpi()
        self.update_gama()
        self.update_relative_permeability()
        self.update_mobilities()
        self.update_transmissibility()
        self.update_loop()
        t1 = time.time()
        dt = t1-t0
        self.update_hist(dt)

        if save:
            self.export_hist()
            self.export_hist2()
            self.mesh.data.update_variables_to_mesh()
            self.mesh.data.export_variables_to_npz()
