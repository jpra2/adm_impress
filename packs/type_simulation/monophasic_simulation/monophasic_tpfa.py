from ...schemes.tpfa_scheme import tpfaScheme
from ...solvers.solvers_scipy.solver_sp import SolverSp

class monophasicTpfa(tpfaScheme):

    def __init__(self, M, data_name: str='monophasicTpfa.npz') -> None:
        super().__init__(M, data_name)
        self.solver = SolverSp()
        self.update_gama()

    def update_gama(self):

        M = self.mesh
        gama = direc.data_loaded['monophasic_data']['gama']
        M.data['gama'] = np.repeat(gama, self.n_volumes)

    def run(self) -> None:

        self.get_transmissibility_matrix_without_boundary_conditions()
        self.get_transmissibility_matrix()
        self.get_RHS_term()

        self.update_initial_pressure_guess(self.solver.direct_solver(self.data['T'], self.data['b']))
        self.get_flux_faces_and_volumes()

        self.mesh.data.update_variables_to_mesh()
        # self.export_data()
