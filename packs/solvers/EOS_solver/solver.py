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
        coef1_square = coef[1]*coef[1]
        Q = 2 * coef1_square*coef[1] / 27 - coef[1] * coef[2] / 3 + coef[3]
        P = - coef1_square / 3 + coef[2]
        delta = (Q / 2)*(Q / 2) + (P / 3)*(P / 3)*(P / 3)
        return Q, delta

    def get_omegas(self):
        omega = (-1 + 1j*np.sqrt(3))/2
        omega_square = omega**2
        omegas  = np.array([[1,1],[omega,omega_square],[omega_square,omega]])
        return omegas

    def get_model_roots(self, omegas, Q, delta):
        xs_args = np.empty([2,len(Q)], dtype=np.complex)

        aux = np.ones(len(Q), dtype=np.complex) #creating for errors that was having in sqrt numpy function
        aux[delta < 0] = 1j
        delta[delta < 0] = -delta[delta < 0]

        delta_sqrtAUX = (delta)**(1/2) * aux
        Q_2 = Q/2
        xs_args[0,:] = - Q_2 + delta_sqrtAUX
        xs_args[1,:] = - Q_2 - delta_sqrtAUX
        real_args = np.isreal(xs_args)
        xs_args[real_args] = np.cbrt(np.real(xs_args[real_args])) + 0j
        xs_args[~real_args] = (xs_args[~real_args])**(1/3)
        #X = omegas@xs_args

        #if len(X.imag[(abs(X.imag)<1e-16) * (abs(X.imag)>0)])>0: import pdb; pdb.set_trace()
        #aux = X.imag
        #aux[abs(X.imag)<1e-16] = 0
        #X.imag = aux

        X_aux = omegas * xs_args.T[:,np.newaxis,:]
        X = X_aux.sum(axis=-1).T
        return X

    def get_actual_roots(self, coef, X):
        roots = X.T - coef[1][:,np.newaxis]/3
        return roots
