import numpy as np
import scipy.sparse
from scipy.sparse import linalg
from implicit_impress.jacobian.symbolic_jacobian import symbolic_J as s_J
from implicit_impress.jacobian.impress_assembly import assembly
import sympy as sym
# import pandas as pd

T, S_up, Sw, So, Swn, Son, Dt, k, phi, p_i, p_j, Dx, Dy=sym.symbols("T S Sw So Swn Son Dt k phi p_i p_j Dx Dy")
class newton():
    def __init__(self,M):
        self.set_properties(M)
        self.pvi_lim=1.0
        self.pvi_acum=0
        self.pvi=0.15
        self.i=0
        self.Dx=1.0
        self.Dy=1.0
        self.Dz=1.0
        self.porous_volume=0.3*self.Dx*self.Dy*self.Dz  #Volume*porosidade

        for j in range(8):
            self.assembly=assembly(M,0.00001)
            self.pvi_def=0.02+0.02*j
            ni=int(self.pvi_lim/self.pvi_def)

            self.qov=[]
            self.qwv=[]
            self.vpiv=[]
            self.pvi_acum=0.0
            self.pvi=0.5
            self.dt=self.porous_volume*self.pvi/self.assembly.vazao

            self.set_properties(M)
            self.iterac=0
            for i in range(ni):
                self.iteration(M)
                self.i+=1
            filepath = 'results/production_data.xlsx'
            try:
                old=numpy.savetxt(filepath,self.F_Jacobian,delimiter=",")
                data=np.hstack([old,self.data])
                print('Appended to an old file')
            except:
                print('First Data file')
                data=self.data
                self.iterac+=1
            # df = pd.DataFrame (data)
            # df.to_excel(filepath, index=False)

    def iteration(self,M):
        n=len(M.volumes.all)
        S0=M.swns[:].T[0]
        P0=M.pressure[:].T[0]
        SP=np.concatenate([S0,P0])
        gids=M.volumes.all
        max_iter=8
        if self.i==0:
            max_iter*=3
        i=0
        csc_matrix=scipy.sparse.csc_matrix
        while i<max_iter and (i==0 or max(abs(sol))>0.001):

            if i>4 or self.i>0:
                self.pvi=self.pvi_def
                self.dt=self.porous_volume*self.pvi/self.assembly.vazao
                self.dt*=36
            self.assembly=assembly(M,self.dt)
            J=self.assembly.J
            q=self.assembly.q
            Fp=-q[0:n]
            Fs=-q[n:]
            J=scipy.sparse.csc_matrix(J)

            Jpp1=J[0:n,0:n]
            Jps1=J[0:n,n:]
            Jsp1=J[n:,0:n]
            Jss1=J[n:,n:]

            Mp1=Jpp1-Jps1*linalg.spsolve(Jss1,Jsp1)

            qp1=(scipy.sparse.csc_matrix(Fp).T-scipy.sparse.csc_matrix(Jps1*linalg.spsolve(Jss1,Fs)).T).toarray().transpose()[0]

            dp1=linalg.spsolve(Mp1,qp1)
            ds1=linalg.spsolve(Jss1,Fs-Jsp1*dp1)

            sol=np.concatenate([dp1,ds1])

            S0+=sol[n:]
            S0[S0<0.2]=0.2
            S0[S0>0.8]=0.8
            P0+=sol[0:n]
            print(max(abs(sol[0:n])),max(abs(sol[n:])))
            M.pressure[:]=P0[gids]
            M.swns[:]=S0[gids]
            i+=1
        v0_facs=M.volumes.adjacencies(M.volumes.all_elements[0])
        v0_internal_faces=np.intersect1d(v0_facs, M.faces.internal_elements[:])
        pv=M.pressure[M.volumes.all_elements[0]]
        # pv=M.mb.tag_get_data(M.n_pressure_tag,M.all_volumes[0],flat=True)[0]

        qo=0.0
        qw=0.0

        for f in v0_internal_faces:
            adjs=M.faces.bridge_adjacencies(f,2,3)
            Sw_up=max(M.swns[adjs])
            # Sw_up=max(M.mb.tag_get_data(M.Swns_tag,np.array(adjs),flat=True))
            pj=max(M.pressure[adjs])
            # pj=max(M.mb.tag_get_data(M.n_pressure_tag,adjs,flat=True))
            qo+=float(self.assembly.F_Jacobian.F_o.subs({T:1.0, k:1.0, Sw:Sw_up, p_i:pv, p_j:pj}))/0.8
            qw+=float(self.assembly.F_Jacobian.F_w.subs({T:1.0, k:1.0, Sw:Sw_up, p_i:pv, p_j:pj}))
        self.pvi_acum+=self.pvi
        self.qov.append(qo)
        self.qwv.append(qw)
        self.vpiv.append(self.pvi_acum)

        self.data=np.array([self.vpiv,self.qov,self.qwv])

        print('oleo, agua',qo,qw,qo+qw,'pvi',self.pvi_acum)
        m4 = M.core.mb.create_meshset()
        # import pdb; pdb.set_trace()
        M.core.mb.add_entities(m4, M.core.all_volumes)
        M.core.mb.write_file('results/biphasic/fully_implicit'+str(self.i)+'.vtk',[m4])
        # import pdb; pdb.set_trace()
        M.swn1s[:]=S0[gids]
        # M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes,S0[gids])

    def set_properties(self,M):
        n=len(M.volumes.all)
        # M.mb.tag_set_data(M.k_eq,self.internal_faces,np.repeat(1.0,len(self.internal_faces)))
        M.k_eq[:]=M.k_harm[:][0]
        # M.mb.tag_set_data(M.phi,M.all_volumes,np.repeat(0.3,len(M.all_volumes)))
        M.phi[:]=np.repeat(0.3,n)
        # M.mb.tag_set_data(M.swns,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))

        M.swns[:]=np.repeat(0.2,n)
        # M.mb.tag_set_data(M.swns,M.all_volumes[-1],0.8)
        M.swns[-1]=0.8
        # M.mb.tag_set_data(M.swn1s,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))
        M.swn1s[:]=np.repeat(0.2,n)
        # M.mb.tag_set_data(M.swn1s,M.all_volumes[-1],0.8)
        M.swn1s[-1]=0.8
        # M.mb.tag_set_data(M.pressure,M.all_volumes,np.array(range(len(M.all_volumes)),dtype=np.float))
        M.pressure[:]=np.arange(n,dtype=np.float)

    def potential_ordering(self, M):
        #Reordena a matriz jacobiana
        d=2
