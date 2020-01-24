from fully_implicit.jacobian.symbolic_jacobian import symbolic_J as F_Jacobian
import sympy as sym
import numpy as np

# import pandas as pd
from sympy import Heaviside as H

T, S_up, Sw, So, Swn, Son, Dt, k, phi, p_i, p_j, Dx, Dy=sym.symbols("T S Sw So Swn Son Dt k phi p_i p_j Dx Dy")
class Ass:
    def __init__(self,M):
        self.internal_faces=[f for f in M.all_faces if len(M.mb.get_adjacencies(f,M.dimension))==2]

        self.boundary_faces=np.setdiff1d(M.all_faces,self.internal_faces)

        self.F_Jacobian=F_Jacobian().J
        # self.set_properties(M)

        self.vazao=1.0
        self.porous_volume=1*0.3  #Volume*porosidade
        self.Dx=1.0/6
        self.Dy=1.0/6
        self.pvi_lim=1.0

        for j in range(8):
            self.pvi_def=0.15-0.02*j
            ni=int(self.pvi_lim/self.pvi_def)
            self.qov=[]
            self.qwv=[]
            self.vpiv=[]
            self.pvi_acum=0.0
            self.pvi=0.02
            self.dt=self.porous_volume*self.pvi/self.vazao
            self.set_properties(M)
            self.i=0
            self.iterac=0
            for i in range(ni):
                self.Newton(M)
                self.i+=1
                filepath = 'output/production_data.xlsx'
            try:
                old=pd.read_excel(filepath)
                data=np.hstack([old,self.data])
                print('Appended to an old file')
            except:
                print('First Data file')
                data=self.data
                self.iterac+=1
            df = pd.DataFrame (data)
            df.to_excel(filepath, index=False)

    def Newton(self,M):
        n=len(M.all_volumes)
        S0=M.mb.tag_get_data(M.Swns_tag,M.all_volumes,flat=True)
        P0=M.mb.tag_get_data(M.n_pressure_tag,M.all_volumes,flat=True)
        SP=np.concatenate([S0,P0])
        gids=M.mb.tag_get_data(M.GLOBAL_ID_tag,M.all_volumes,flat=True)
        max_iter=8
        if self.i==0:
            max_iter*=3
        i=0

        while i<max_iter and (i==0 or max(abs(sol[n:]))>0.001):
            if i>4 or self.i>0:
                self.pvi=self.pvi_def
                self.dt=self.porous_volume*self.pvi/self.vazao
            J, q=self.get_jacobian_matrix(M)
            J[0]=0.0
            J[0][0]=1.0

            J[n,:]=0.0
            # J[:,n]=0.0
            J[n,n]=1.0
            q[n]=0.0

            q[0]=0.0
            q[-1]-=self.vazao
            # q[n-1]-=2.0
            # Ab=np.zeros((2*n,2*n+1),dtype=float)
            # Ab[:,:-1]+=J
            # Ab[:,-1]+=q
            # df = pd.DataFrame (Ab)
            # self.iterac+=1
            # filepath = 'output/Jacob_mod'+str(self.i)+str(i)+'.xlsx'
            # df.to_excel(filepath, index=False)

            Jpp=J[0:n,0:n]
            Jps=J[0:n,n:]
            Jsp=J[n:,0:n]
            Jss=J[n:,n:]

            Fp=-q[0:n]
            Fs=-q[n:]

            Mp=Jpp-np.dot(Jps,np.linalg.solve(Jss,Jsp))
            qp=Fp-np.dot(Jps,np.linalg.solve(Jss,Fs))
            dp=np.linalg.solve(Mp,qp)
            ds=np.linalg.solve(Jss,Fs-np.dot(Jsp,dp))

            sol=np.concatenate([dp,ds])
            # sol2=np.linalg.solve(J,-q)
            # import pdb; pdb.set_trace()
            # sol=np.linalg.solve(J,-q)
            S0+=sol[n:]
            pos_men=np.where(S0<0.2)[0]
            pos_mai=np.where(S0>0.8)[0]
            S0[pos_men]=0.2
            S0[pos_mai]=0.8
            P0+=sol[0:n]
            print(max(abs(sol[0:n])),max(abs(sol[n:])))

            M.mb.tag_set_data(M.n_pressure_tag,M.all_volumes,P0[gids])
            M.mb.tag_set_data(M.Swns_tag,M.all_volumes,S0[gids])
            i+=1

        # Traçando curva de produção
        v0_facs=M.mb.get_adjacencies(M.all_volumes[0], M.dimension-1)
        v0_internal_faces=np.setdiff1d(v0_facs, self.boundary_faces)
        pv=M.mb.tag_get_data(M.n_pressure_tag,M.all_volumes[0],flat=True)[0]

        qo=0.0
        qw=0.0

        for f in v0_internal_faces:

            adjs=M.mb.get_adjacencies(f,M.dimension)
            Sw_up=max(M.mb.tag_get_data(M.Swns_tag,np.array(adjs),flat=True))
            pj=max(M.mb.tag_get_data(M.n_pressure_tag,adjs,flat=True))
            qo+=float(F_Jacobian().F_o.subs({T:1.0, k:1.0, Sw:Sw_up, p_i:pv, p_j:pj}))/0.8
            qw+=float(F_Jacobian().F_w.subs({T:1.0, k:1.0, Sw:Sw_up, p_i:pv, p_j:pj}))
        self.pvi_acum+=self.pvi
        self.qov.append(qo)
        self.qwv.append(qw)
        self.vpiv.append(self.pvi_acum)

        self.data=np.array([self.vpiv,self.qov,self.qwv])

        print('oleo, agua',qo,qw,qo+qw,'pvi',self.pvi_acum)
        m4 = M.mb.create_meshset()
        M.mb.add_entities(m4, M.all_volumes)
        M.mb.write_file('results/biphasic/fully_implicit'+str(self.i)+'.vtk',[m4])
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes,S0[gids])

        # import pdb; pdb.set_trace()




    def get_jacobian_matrix(self,M):
        GID_volumes=M.mb.tag_get_data(M.GLOBAL_ID_tag,M.all_volumes, flat=True)
        count=0
        J=np.zeros((2*len(M.all_volumes),2*len(M.all_volumes)),dtype=float)
        q=np.zeros(2*len(M.all_volumes))
        Swns=M.mb.tag_get_data(M.Swns_tag,M.all_volumes,flat=True)
        Swn1s=M.mb.tag_get_data(M.Swn1s_tag,M.all_volumes,flat=True)
        n=len(M.all_volumes)
        for v in M.all_volumes:
            pv=M.mb.tag_get_data(M.n_pressure_tag,v,flat=True)[0]
            faces=M.mb.get_adjacencies(v, M.dimension-1)
            internal_faces=np.setdiff1d(faces, self.boundary_faces)
            Swn_v=Swns[count]
            ID_vol=GID_volumes[count]

            J[ID_vol][n+ID_vol]+=float(F_Jacobian().c_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))
            J[n+ID_vol][n+ID_vol]+=float(F_Jacobian().c_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt}))
            q[ID_vol]+=float(F_Jacobian().acum_o.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))
            q[n+ID_vol]+=float(F_Jacobian().acum_w.subs({Dx:self.Dx, Dy:self.Dy, phi:0.3, Dt:self.dt, Sw:Swns[count], Swn:Swn1s[count]}))


            adjs=np.array([M.mb.get_adjacencies(face, M.dimension) for face in internal_faces])
            
            adjs0=adjs[:,0]
            adjs1=adjs[:,1]

            ids0=M.mb.tag_get_data(M.GLOBAL_ID_tag,np.array(adjs0),flat=True)
            ids1=M.mb.tag_get_data(M.GLOBAL_ID_tag,np.array(adjs1),flat=True)

            Swns0=M.mb.tag_get_data(M.Swns_tag,np.array(adjs0),flat=True)
            Swns1=M.mb.tag_get_data(M.Swns_tag,np.array(adjs1),flat=True)

            press0=M.mb.tag_get_data(M.n_pressure_tag,np.array(adjs0),flat=True)
            press1=M.mb.tag_get_data(M.n_pressure_tag,np.array(adjs1),flat=True)

            count_fac=0

            for f in internal_faces:
                pf0=press0[count_fac]
                pf1=press1[count_fac]

                swf=Swns0[count_fac] if pf0>pf1 else Swns1[count_fac]
                id_up=ids0[count_fac] if pf0>pf1 else ids1[count_fac]

                pj=press0[count_fac] if press1[count_fac]==pv else press1[count_fac]
                id_j=ids1[count_fac] if ids0[count_fac]==ID_vol else ids0[count_fac]
                if id_j==ID_vol:
                    print("ERRO FATAL!!!")

                # id_up=ID_vol if pv>pj else id_j

                # if id_up2!=id_up and pv!=pj:
                    # print("id_up2!=id_up!")




                J00=float(self.F_Jacobian[0][0].subs({T:1, k:1, Sw:swf}))
                J01=float(self.F_Jacobian[0][1].subs({T:1, k:1, Sw:swf, p_i:pv, p_j:pj}))
                J10=float(self.F_Jacobian[1][0].subs({T:1, k:1, Sw:swf}))
                J11=float(self.F_Jacobian[1][1].subs({T:1, k:1, Sw:swf, p_i:pv, p_j:pj}))

                # if id_up<ID_vol:
                #     print("id_up<ID_vol! J01, J11, ID_vol, id_up: ", J01, J11, ID_vol, id_up,' self.i: ', self.i)
                    # import pdb; pdb.set_trace()

                q[ID_vol]-=float(F_Jacobian().F_o.subs({T:1.0, k:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))

                q[n+ID_vol]-=float(F_Jacobian().F_w.subs({T:1.0, k:1.0, Sw:Swns1[count_fac], p_i:pv, p_j:pj}))

                J[ID_vol][ID_vol]-=J00
                J[ID_vol][id_j]+=J00

                J[n+ID_vol][ID_vol]-=J10
                J[n+ID_vol][id_j]+=J10

                J[ID_vol][n+id_up]-=J01
                J[n+ID_vol][n+id_up]-=J11

                count_fac+=1

            count+=1
        # df = pd.DataFrame (J)
        self.iterac+=1
        # filepath = 'output/my_excel_file'+str(self.iterac)+'.xlsx'
        # df.to_excel(filepath, index=False)
        # print('salvou a jacobiana '+filepath)
        return(J, q)

    def set_properties(self,M):
        M.mb.tag_set_data(M.GLOBAL_ID_tag,M.all_volumes,range(len(M.all_volumes)))

        M.mb.tag_set_data(M.GLOBAL_ID_tag,M.all_volumes[5],2)
        M.mb.tag_set_data(M.GLOBAL_ID_tag,M.all_volumes[2],5)

        M.mb.tag_set_data(M.k_eq_tag,self.internal_faces,np.repeat(1.0,len(self.internal_faces)))
        M.mb.tag_set_data(M.phi_tag,M.all_volumes,np.repeat(0.3,len(M.all_volumes)))
        M.mb.tag_set_data(M.press_value_tag,np.uint64([M.all_volumes[0],M.all_volumes[-1]]),[0.0,15.0])
        M.mb.tag_set_data(M.Swns_tag,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))
        M.mb.tag_set_data(M.Swns_tag,M.all_volumes[-1],0.8)
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes[-1],0.8)
        M.mb.tag_set_data(M.n_pressure_tag,M.all_volumes,np.array(range(len(M.all_volumes)),dtype=np.float))
