from .biphasic_tpfa import BiphasicTpfa
import numpy as np
import sympy as sy
import scipy.sparse as sp

import time

class Gtime:
    ident = 0
    ts = [0, -1]

    def __call__(self):
        ident = self.__class__.ident
        self.__class__.ts[ident] = time.time()
        if self.__class__.ident == 0:
            self.__class__.ident = 1
        else:
            self.__class__.ident = 0

    def dt(self):
        ts = self.__class__.ts
        return ts[1] - ts[0]

class BiphasicSfi(BiphasicTpfa):

    def __init__(self, M, data_impress, elements_lv0, wells, data_name: str='BiphasicSfi.npz'):
        super().__init__(M, data_impress, elements_lv0, wells, data_name)

        self.create_symbolic_variables()
        self.tol_sat = 1e-3
        self.max_iter = 1000
        # self.delta_t = 0.0009
        self.initialize_iter_erro()

    def initialize_iter_erro(self):
        self.n_print = 0
        self.erro = 1
        self.iter = 1
        self.dt = 0
        self.rodar = True
        self.update_saturation_last = True
        self.ident_for_update_delta_t = True

    def create_symbolic_variables(self):

        self.symbolic_saturation = sy.symbols('Sw')
        self.symbolic_normalized_saturation = self.relative_permeability._stemp(self.symbolic_saturation)
        self.symbolic_krw = self.relative_permeability._krw(self.symbolic_normalized_saturation)
        self.symbolic_kro = self.relative_permeability._kro(self.symbolic_normalized_saturation)
        self.symbolic_lambda_w = self.symbolic_krw/self.biphasic_data['mi_w']
        self.symbolic_lambda_o = self.symbolic_kro/self.biphasic_data['mi_o']
        self.symbolic_lambda_t = self.symbolic_lambda_w + self.symbolic_lambda_o
        self.symbolic_fw = self.symbolic_lambda_w/self.symbolic_lambda_t
        self.symbolic_dfw = self.symbolic_fw.diff(self.symbolic_saturation)

    def identificate_volumes_for_faces_upwind(self):

        internal_faces = self.elements_lv0['internal_faces']
        flux_faces = self.data_impress['flux_faces']
        flux_internal_faces = flux_faces[internal_faces]

        self._data['sw_faces_identificate_vols'] = np.full((len(internal_faces), 2), False, dtype=bool)

        ids = np.arange(len(internal_faces))
        fluxo_positivo = ids[flux_internal_faces <= 0]
        outros = np.setdiff1d(ids, fluxo_positivo)

        self._data['sw_faces_identificate_vols'][fluxo_positivo, 1] = np.full(len(fluxo_positivo), True, dtype=bool)
        self._data['sw_faces_identificate_vols'][outros, 0] = np.full(len(outros), True, dtype=bool)

    def get_jacobian_for_saturation(self):

        internal_faces = self.elements_lv0['internal_faces']
        sw_faces_identificate_vols = self._data['sw_faces_identificate_vols']
        neig_internal_faces = self.elements_lv0['neig_internal_faces']
        flux_faces = self.data_impress['flux_faces']
        sw = self.data_impress['saturation']
        flux_volumes = self.data_impress['flux_volumes']
        vol_volumes = self.data_impress['volume']
        volumes = self.elements_lv0['volumes']
        phis = self.data_impress['poro']
        fws_vol = self.data_impress['fw_vol']
        flux_w_volumes = self.data_impress['flux_w_volumes']

        # import pdb; pdb.set_trace()
        #
        # wsp = self.wells['ws_p']
        # wsq = self.wells['ws_q']

        lines = []
        cols = []
        data = []
        arrl = np.zeros(2, dtype=int)
        arrc = arrl.copy()
        arrd = arrl.copy().astype(float)

        volumes_sw_face = neig_internal_faces[sw_faces_identificate_vols]
        sws_faces = sw[volumes_sw_face]
        flux_internal_faces = flux_faces[internal_faces]
        dfws = np.zeros(len(internal_faces))

        # for i, face in enumerate(internal_faces):
        #     vs = neig_internal_faces[i]
        #     ident = sw_faces_identificate_vols[i]
        #     volume_sw_face = vs[ident][0]
        #     sw_face = sw[volume_sw_face]
        #     flux_face = flux_faces[face]
        #     dfw = self.symbolic_dfw.subs(self.symbolic_saturation, sw_face)
        #
        #     arrl[:] = vs
        #     arrc[:] = [volume_sw_face, volume_sw_face]
        #     arrd[0] = dfw*flux_face
        #     arrd[1] = -dfw*flux_face
        #
        #     lines.append(arrl)
        #     cols.append(arrc)
        #     data.append(arrd)

        for i, face in enumerate(internal_faces):
            dfws[i] = (self.symbolic_dfw.subs(self.symbolic_saturation, sws_faces[i]))

        # k1 = -1
        k1 = 1

        lines.append(neig_internal_faces.flatten())
        cols.append(np.array([volumes_sw_face, volumes_sw_face]).T.flatten())
        data.append(np.array([k1*dfws*flux_internal_faces, -k1*dfws*flux_internal_faces]).T.flatten())

        qw = fws_vol*flux_volumes

        # arrl = np.array([], dtype=int)
        # arrc = arrl.copy()
        # arrd = arrl.copy().astype(float)

        # for vol in volumes:
        #     vol_vol = vol_volumes[vol]
        #     phi_vol = phis[vol]
        #     indep_term = (phi*vol_vol/self.delta_t)*(1-sw[vol])
        #
        #     arrl[0] = vol
        #     arrc[0] = vol
        #     arrd[0] = indep_term
        #
        #     lines.append(arrl)
        #     cols.append(arrc)
        #     data.append(arrd)

        # k0 = -1
        k0 = 1

        indep_term = k0*(1/self.delta_t)*(phis*vol_volumes)*(1 - sw)
        lines.append(volumes)
        cols.append(volumes)
        data.append(indep_term)

        n = len(volumes)

        lines = np.concatenate(lines)
        cols = np.concatenate(cols)
        data = np.concatenate(data)

        J = sp.csc_matrix((data, (lines, cols)), shape=(n, n))

        # if self.iter > 1:
        #     import pdb; pdb.set_trace()
        # import pdb; pdb.set_trace()

        fx = k0*(1/self.delta_t)*(phis*vol_volumes)*(self.data_impress['saturation'] - self.data_impress['saturation_last']) - flux_w_volumes + qw

        return J, fx

    def printar(self):
        self.n_print = 0
        self.data_impress.update_variables_to_mesh(['saturation'])
        self.mesh.core.print(folder='results', file='test'+ str(self.n_print), extension='.vtk', config_input='input_cards/print_settings0.yml')
        self.n_print += 1

    def run(self):
        if self.update_saturation_last:
            self.data_impress['saturation_last'] = self.data_impress['saturation'].copy()
            self.update_saturation_last = False

        gt = Gtime()
        gt()
        super().run()
        self.identificate_volumes_for_faces_upwind()
        gt()
        self.dt += gt.dt()

    def run_2(self, save=False):

        gt = Gtime()
        gt()

        if not self.rodar:
            self.update_flux_w_and_o_volumes()
            self.update_t()
            self.update_vpi()
            self.update_loop()
            gt()
            self.dt += gt.dt()
            self.update_current_biphasic_results(self.dt)
            if save:
                self.save_infos()
            self.initialize_iter_erro()
            return 0

        tol = self.tol_sat
        n_max_iter = 100
        verif = True
        self.update_flux_w_and_o_volumes()
        if self.ident_for_update_delta_t:
            self.update_delta_t()
            self.ident_for_update_delta_t = False

        # import pdb; pdb.set_trace()

        # self._data['sw_ni'] = self.data_impress['saturation'].copy()
        J, fx = self.get_jacobian_for_saturation()
        dx = self.solver.direct_solver(J, fx)
        self.data_impress['saturation'] += dx
        import pdb; pdb.set_trace()
        # erro = np.absolute(self.data_impress['saturation'] - self._data['sw_ni']).max()
        self.erro = np.absolute(dx).max()

        if self.iter > self.max_iter:
            print('Max iter')
            self.rodar = False

        if self.erro < tol:
            self.rodar = False

        self.iter += 1

        self.update_relative_permeability()
        self.update_mobilities()
        self.update_transmissibility()

        gt()
        self.dt += gt.dt()

        self.printar()
        import pdb; pdb.set_trace()

        return 1










        import pdb; pdb.set_trace()
