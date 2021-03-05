from ..flux_schemes.tpfa_scheme import TpfaScheme
from ..flux_calculation.flux_tpfa import TpfaFlux
from ..directories import data_loaded
from ..solvers.solvers_scipy.solver_sp import SolverSp
import numpy as np

class FineScaleTpfaPressureSolver(TpfaScheme, TpfaFlux):

    def __init__(self, data_impress, elements_lv0, wells, data_name: str='FineScaleTpfaPressureSolver.npz', load=False):
        super().__init__(data_impress, elements_lv0, data_name=data_name, load=load)
        self.gravity = data_loaded['gravity']
        self.wells = wells
        self.solver = SolverSp()

    def get_RHS_term(self):

        b = self.data_impress[self.data_impress.variables_impress['flux_grav_volumes']].copy()

        if self.gravity:

            self.wells.add_gravity()

        ws_p = self.wells['ws_p']  # pocos de pressao prescrita
        values_p = self.wells['values_p']  # valores de pressao prescrita
        ws_q = self.wells['ws_q']  # pocos de vazao prescrita
        values_q = self.wells['values_q']  # valor da vazao prescrita

        b[ws_q] += values_q
        b[ws_p] = values_p

        return b

    def set_boundary_conditions(self):
        '''
        transmissibility matrix with boundary conditions
        '''
        ws_p = self.wells['ws_p']

        T = self['Tini'].tolil().copy()
        T[ws_p] = np.zeros((len(ws_p), T.shape[0]))
        T[ws_p, ws_p] = np.ones(len(ws_p))

        return T.tocsc()

    def update_gama(self):

        gama = data_loaded['monophasic_data']['gama']
        self.data_impress['gama'] = np.repeat(gama, self.n_volumes)

    def run(self):
        self.update_gama()
        self.get_transmissibility_matrix_without_boundary_conditions()
        self.get_gravity_source_term()
        T = self.set_boundary_conditions()
        b = self.get_RHS_term()
        return T, b
        # p = self.solver.direct_solver(T, b)
        # self.data_impress['pressure'] = p
        # self.get_flux_faces_and_volumes()
