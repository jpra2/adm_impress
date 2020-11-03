import numpy as np
from cmath import acos
'''Class created to get cubic roots in the analytical way, so its possible to vectorize some stuff that involves the EOS cubic solution'''
'''Here we present the method/notations present in Chen's book, Computational Methods for Multiphase Flow...'''

class CubicRoots:
    def __init__(self):
        pass

    def run(self, coef):
        Q, delta = self.get_delta(coef)
        omegas = self.get_omegas()
        X = self.get_model_roots(omegas, Q, delta)
        real_roots = self.get_actual_roots(coef, X)
        return real_roots

    def get_delta(self, coef):
        Q = 2 * coef[1]**3 / 27 - coef[1] * coef[2] / 3 + coef[3]
        P = - coef[1]**2 / 3 + coef[2]
        delta = (Q / 2)**2 + (P / 3)**3
        return Q, delta

    def get_omegas(self):
        omega = (-1 + 1j*np.sqrt(3))/2
        omegas  = np.array([[1,1],[omega,omega**2],[omega**2,omega]])
        return omegas

    def get_model_roots(self, omegas, Q, delta):
        xs_args = np.empty([2,len(Q)], dtype=np.complex)

        aux = np.ones(len(Q), dtype=np.complex) #creating for errors that was having in sqrt numpy function
        aux[delta < 0] = 1j
        delta[delta < 0] = -delta[delta < 0]
        
        xs_args[0,:] = - Q/2 + (delta)**(1/2) * aux
        xs_args[1,:] = - Q/2 - (delta)**(1/2) * aux
        real_args = np.isreal(xs_args)
        xs_args[real_args] = np.cbrt(np.real(xs_args[real_args])) + 0j
        xs_args[~real_args] = (xs_args[~real_args])**(1/3)
        X = omegas@xs_args
        return X

    def get_actual_roots(self, coef, X):
        roots = X.T - coef[1][:,np.newaxis]/3
        return roots
