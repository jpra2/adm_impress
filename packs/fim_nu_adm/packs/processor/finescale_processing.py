import numpy as np
import scipy.sparse as sp
import time
from scipy.sparse import linalg
# from ..postprocessor.exporter import FieldVisualizer
from .assembler import Assembler
# visualize=FieldVisualizer()

class NewtonIterationFinescale():
    def __init__(self, wells, faces, volumes, case):
        self.Assembler = Assembler(wells, faces, volumes)
        self.time_solve=[]
        self.porosities=volumes['pore_volume']
        self.wells=wells
        self.viz=np.load('results/'+case+'/viz.npy')

    # @profile
    def newton_iteration_finescale(self, p, s, time_step, rel_tol=1e-3):
        pressure = p.copy()
        swns = s.copy()
        swn1s = s.copy()
        converged=False
        count=0
        dt=time_step
        while not converged:
            # swns[self.Assembler.wells['ws_inj']]=1
            self.Assembler.iteration=count
            J, q=self.Assembler.get_jacobian_matrix(swns, swn1s, pressure, time_step)
            t0=time.time()
            # print("resolvendo ...")
            sol=-linalg.spsolve(J, q)
            self.time_solve.append(time.time()-t0)
            # print("resolveu! ", time.time()-t0)
            # import pdb; pdb.set_trace()
            n=int(len(q)/2)
            pressure+=sol[0:n]
            sol[n+self.wells['ws_prod']]=0
            swns+=sol[n:]
            converged=max(abs(sol[n:]))<rel_tol and (swns.max()<1.00000001 and swns.min()>-0.000001)
            swns[swns>1]=1
            swns[swns<0]=0
            # print(count, max(abs(sol[n:])),'fs')
            print(count,f'{ max(abs(sol[n:])):.2e}')
            count+=1
            if count>20 or max(abs(sol[n:])>3):
                count=21
                print('excedded maximum number of iterations finescale')
                return False, count, pressure, swns
        # saturation[wells['ws_prod']]=saturation[wells['viz_prod']].sum()/len(wells['viz_prod'])
        swns[self.wells['ws_prod']]=swns[self.viz].sum()/len(self.viz)
        self.PVI=(swns*self.porosities).sum()/self.porosities.sum()
        return True, count, pressure, swns
