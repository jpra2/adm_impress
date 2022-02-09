import numpy as np
from ..directories import data_loaded

def debug_verificate(func):

    def verificate(*args1, **kwargs):
        if data_loaded['_debug'] == True:
            func(*args1, **kwargs)
        else:
            return 0

    return verificate

class testsGeneral:

    @debug_verificate
    def test_rhs_term(self, T, b, p):
        '''
            T: transmissibility with boundary conditions
            b: rhs term
            p: pressure
        '''
        b2 = T*p
        assert np.allclose(b2, b) == True

    @debug_verificate
    def test_flux_volumes(self, T_without, p, gravity_source_term_volumes, flux_volumes):
        '''
            T_witout: transmissibility matrix without boundary conditions
            p: pressure
            gravity_source_term_volumes: autoexplicativo
            flux_volumes: fluxo nos volumes calculado a partir do fluxo nas faces
        '''
        assert np.allclose(T_without*p - gravity_source_term_volumes, flux_volumes) == True
