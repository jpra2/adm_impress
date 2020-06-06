import numpy as np
from cmath import acos
'''Class created to get cubic roots in the analytical way, so its possible to vectorize some stuff that involves the EOS cubic solution'''
'''Here we present the method/notations present in Chen's book, Computational Methods for Multiphase Flow...'''

class CubicRoots:
    def __init__(self):
        pass

    def run(self, coef):
        self.get_delta(coef)
        self.get_omegas()
        self.get_model_roots()
        real_roots = self.get_actual_roots(coef)
        return real_roots

    def get_delta(self, coef):
        self.Q = 2 * coef[1]**3 / 27 - coef[1] * coef[2] / 3 + coef[3]
        P = - coef[1]**2 / 3 + coef[2]
        self.delta = (self.Q / 2)**2 + (P / 3)**3

    def get_omegas(self):
        omega = (-1 + 1j*np.sqrt(3))/2
        self.omegas  = np.array([[1,1],[omega,omega**2],[omega**2,omega**3]])

    def get_model_roots(self):
        xs_args = np.empty([2,len(self.Q)])
        xs_args[0,:] = - self.Q/2 + self.delta**(1/2)
        xs_args[1,:] = - self.Q/2 - self.delta**(1/2)

        real_args = np.isreal(xs_args)
        real_index = np.where(real_args == True)[0]
        complex_index2 = np.where(real_args ==  False)[0]
        xs_args[real_index] = np.cbrt(xs_args[real_index])
        xs_args[complex_index2] = (xs_args[complex_index2])**(1/3)

        self.X = self.omegas@xs_args

    def get_actual_roots(self, coef):
        reais = np.isreal(self.X)
        r_pos = np.where(reais==True)
        Xreais = np.real(self.X[r_pos[:]])
        real_roots = Xreais - coef[1]/3
        return real_roots
