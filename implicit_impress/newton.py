
class newton():
    def __init__(self):
        d=1
    def iteration(self,M):
        n=len(M.all_volumes)
        S0=M.mb.tag_get_data(M.Swns_tag,M.all_volumes,flat=True)
        P0=M.mb.tag_get_data(M.n_pressure_tag,M.all_volumes,flat=True)
        SP=np.concatenate([S0,P0])
        gids=M.mb.tag_get_data(M.GLOBAL_ID_tag,M.all_volumes,flat=True)
        max_iter=8
        if self.i==0:
            max_iter*=3
        i=0
        csc_matrix=scipy.sparse.csc_matrix
        while i<max_iter and (i==0 or max(abs(sol))>0.001):

            if i>4 or self.i>0:
                self.pvi=self.pvi_def
                self.dt=self.porous_volume*self.pvi/self.vazao
                self.dt*=40

            J, q=self.get_jacobian_matrix(M)
            J2=self.apply_dirichlet(J,[0,n])

            q[n]=0.0

            q[0]=0.0
            q[-1]-=self.vazao

            Fp=-q[0:n]
            Fs=-q[n:]
            J1=scipy.sparse.csc_matrix(J2)

            Jpp1=J1[0:n,0:n]
            Jps1=J1[0:n,n:]
            Jsp1=J1[n:,0:n]
            Jss1=J1[n:,n:]

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

            M.mb.tag_set_data(M.n_pressure_tag,M.all_volumes,P0[gids])
            M.mb.tag_set_data(M.Swns_tag,M.all_volumes,S0[gids])
            i+=1

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

    def set_properties(self,M):
        M.mb.tag_set_data(M.GLOBAL_ID_tag,M.all_volumes,range(len(M.all_volumes)))
        M.mb.tag_set_data(M.k_eq_tag,self.internal_faces,np.repeat(1.0,len(self.internal_faces)))
        M.mb.tag_set_data(M.phi_tag,M.all_volumes,np.repeat(0.3,len(M.all_volumes)))
        M.mb.tag_set_data(M.press_value_tag,np.uint64([M.all_volumes[0],M.all_volumes[-1]]),[0.0,15.0])
        M.mb.tag_set_data(M.Swns_tag,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))
        M.mb.tag_set_data(M.Swns_tag,M.all_volumes[-1],0.8)
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes,np.repeat(0.2,len(M.all_volumes)))
        M.mb.tag_set_data(M.Swn1s_tag,M.all_volumes[-1],0.8)
        M.mb.tag_set_data(M.n_pressure_tag,M.all_volumes,np.array(range(len(M.all_volumes)),dtype=np.float))

    def potential_ordering(self, M):
        #Reordena a matriz jacobiana
        d=2
