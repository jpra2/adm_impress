from packs.compositional.pressure_solver import TPFASolver
from packs.utils.test_functions import test_kwargs_keys


class AdmTpfaCompositionalSolver(TPFASolver):
    
    _kwargs_keys = {
        'get_pressure': set(['multilevel_data',
                             'multilevel_operators'])
    }
    
    def get_pressure(self, M, wells, fprop, delta_t, **kwargs):
        test_kwargs_keys(AdmTpfaCompositionalSolver._kwargs_keys['get_pressure'], kwargs.keys())
        T = self.update_transmissibility(M, wells, fprop, delta_t)
        D = self.update_independent_terms(M, fprop, wells, delta_t)
        self.P = self.update_pressure(T, D, fprop)
        Ft_internal_faces = self.update_total_flux_internal_faces(M, fprop)
        self.update_flux_wells(fprop, wells, delta_t)
        return self.P, Ft_internal_faces, self.q