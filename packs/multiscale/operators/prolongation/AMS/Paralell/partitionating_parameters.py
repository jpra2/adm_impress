from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt
import sympy
import numpy as np
import time

class calibrate_partitioning_parameters:
    def __init__(self):
        As=np.array([10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500])
        bs=np.array([1, 8, 16, 24, 32, 40])
        Abt=self.generate_A_b_t(As, bs)

        # Abt=np.load("flying/A_b_t.npy")
        xs=Abt[:,0]
        ys=Abt[:,1]
        zs=Abt[:,2]

        coefs = self.calculate_coeficients(xs, ys, zs)
        np.save("flying/partitioning_coeffitients_bx_cy_dx2_exy_fy2_intercept.npy",coefs)
        fig = plt.figure()
        self.imprima3d2(xs, ys,zs, fig)
        lims=[xs.min(), xs.max(), ys.min(), ys.max()]
        self.reg3(coefs,lims, fig)
        plt.close('all')


    def calculate_coeficients(self, xs, ys, zs):
        xx=np.array([xs])
        dx2=(abs(xx.T-xx)**2).sum(axis=1)

        yy=np.array([ys])
        dy2=(abs(yy.T-yy)**2).sum(axis=1)
        d2=dx2+dy2
        d=d2**0.5
        pesos=d/d.sum()

        degree=2
        poly = PolynomialFeatures(degree=degree)
        X_t = poly.fit_transform(np.array([xs,ys]).T,pesos)

        clf = LinearRegression()
        clf.fit(X_t, zs)
        coefs = np.concatenate([np.array(clf.coef_)[1:],np.array([clf.intercept_])])
        return coefs

    def imprima3d2(self, xs, ys, zs, fig):
        ax = fig.add_subplot(111, projection='3d')
        m="o"
        ax.scatter(xs, ys, zs, marker=m)
        ax.set_xlabel('A')
        ax.set_ylabel('b')
        ax.set_zlabel('time_to_solve (s)')
        plt.savefig("results/teste_scatter3d.png")

    def reg3(self, coefsi, lims, fig):
        coefs=coefsi[:-1]
        intercept=coefsi[-1]
        x, y = sympy.symbols('x y')
        degree =2
        if degree==2:
            syms=np.array([x, y, x*x, x*y, y*y])
        if degree==1:
            syms=np.array([x, y])
        func=(coefs*syms).sum()+intercept

        ax = fig.add_subplot(111, projection='3d')

        ax.set_xlabel('A')
        ax.set_ylabel('b')
        ax.set_zlabel('time_to_solve (s)')
        plt.show()

        sympy.plotting.plot3d(func, (x, lims[0], lims[1]), (y, lims[2], lims[3]))
        plt.savefig("results/sym_surface2.png")

    def generate_A_b_t(self, As, bs):
        from scipy.sparse import csc_matrix
        from scipy.sparse import linalg
        A_b_t=[]
        for A in As:
            for b in bs:
                AA=np.random.rand(A,A)
                AA=csc_matrix(AA)
                bb=np.random.rand(A,b)
                bb=csc_matrix(bb)
                t0=time.time()
                linalg.spsolve(AA,bb)
                t=time.time()-t0

                A_b_t.append(np.array([A, b, t]))
        return(np.array(A_b_t))
