import numpy as np
from packs.stokes_brinkman.assembly import assembly
class stokes_solver:
    def __init__(self, M):
        dx, dy, dz = self.get_mesh_properties(M)
        self.nv=len(M.volumes.all)
        self.nfi=len(M.faces.internal)
        self.initiate_f_int_tag(M)
        self.M = assembly.global_assembly(M, dx, dy, dz).M
        self.rhs=np.zeros(self.nv+self.nfi)
        self.rhs[0]=0
        self.rhs[self.nv-1]=1
        np.savetxt("results/mat.csv",self.M,delimiter=",")
        My=self.M[0:self.nv,self.nv:self.nv+6].copy()
        Mx=self.M[0:self.nv,self.nv+6:self.nv+12].copy()
        self.M[0:self.nv,self.nv:self.nv+6]=Mx.copy()
        self.M[0:self.nv,self.nv+6:self.nv+12]=My.copy()



        My=self.M[self.nv:self.nv+6,self.nv:self.nv+6].copy()
        Mx=self.M[self.nv+6:self.nv+12,self.nv+6:self.nv+12].copy()
        # import pdb; pdb.set_trace()
        self.M[self.nv:self.nv+6,self.nv:self.nv+6]=Mx.copy()
        self.M[self.nv+6:self.nv+12,self.nv+6:self.nv+12]=My.copy()
        np.savetxt("results/mat2.csv",self.M,delimiter=",")
        import pdb; pdb.set_trace()
        self.solve(self.M,self.rhs)
    def get_mesh_properties(self,M):
        v0=M.volumes.all[0]
        vert_v0=M.volumes.bridge_adjacencies(v0,0,0)
        coords_vert=M.nodes.coords(vert_v0)
        dx, dy, dz=coords_vert.max(axis=0)-coords_vert.min(axis=0)
        return dx, dy, dz

    def initiate_f_int_tag(self, M):

        faces=M.faces.internal
        nodes=M.faces.bridge_adjacencies(faces,0,0).flatten()
        coords=M.nodes.coords[nodes].reshape(len(faces),4,3)
        deltas_x_fac=coords[:,:,0].max(axis=1)-coords[:,:,0].min(axis=1)
        horizontal=faces[deltas_x_fac>10**-8]
        vertical=faces[deltas_x_fac<10**-8]
        n_h=len(horizontal)
        n_v=len(vertical)
        M.id_fint[horizontal] = range(n_h)
        M.id_fint[vertical] = range(n_h,n_h+n_v)

        # int_faces=M.faces.internal
        # n_int_faces=len(int_faces)
        # M.id_fint[int_faces] = range(n_int_faces)#np.arange(n_int_faces)
    def solve(self,Mat, rhs):

        Mat[0]=np.zeros(self.nv+self.nfi)
        Mat[0][0]=1
        Mat[self.nv-1]=np.zeros(self.nv+self.nfi)
        Mat[self.nv-1][self.nv-1]=1


        sol=np.linalg.solve(self.M,self.rhs)
        import pdb; pdb.set_trace()
