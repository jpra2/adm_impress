import numpy as np
from packs.stokes_brinkman_3d.assembly import assembly
class stokes_solver:
    def __init__(self, M):
        dx, dy, dz = self.get_mesh_properties(M)
        self.initiate_f_int_tag(M)
        M.id_global[:]=M.volumes.all

        v=M.core.mb.create_meshset()
        M.core.mb.add_entities(v,M.core.all_volumes)
        M.core.mb.write_file("results/solution_volumes.vtk",[v])
        f=M.core.mb.create_meshset()
        M.core.mb.add_entities(f,M.core.all_faces[self.fx])
        M.core.mb.write_file("results/solution_faces_x.vtk",[f])

        self.nv=len(M.volumes.all)
        self.nfi=len(M.faces.internal)

        self.M = assembly.global_assembly(self,M, dx, dy, dz).M
        self.rhs=np.zeros(self.nv+self.nfi)
        self.rhs[0]=1
        self.rhs[self.nv-1]=0
        np.savetxt("results/mat3d.csv",self.M,delimiter=",")

        #import pdb; pdb.set_trace()
        self.sol=self.solve(self.M,self.rhs)
        self.write_output(M)

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

        cx=coords[:,:,0]
        cy=coords[:,:,1]
        cz=coords[:,:,2]
        fx=faces[(cx.max(axis=1)-cx.min(axis=1))<10**-2]
        fy=faces[(cy.max(axis=1)-cy.min(axis=1))<10**-2]
        fz=faces[(cz.max(axis=1)-cz.min(axis=1))<10**-2]

        self.fx=fx
        self.fy=fy
        self.fz=fz

        n_x=len(fx)
        n_y=len(fy)
        n_z=len(fz)

        M.id_fint[fx] = range(n_x)
        M.id_fint[fy] = range(n_x,n_x+n_y)
        if n_z>0:
            M.id_fint[fz] = range(n_x+n_y,n_x+n_y+n_z)
        f=M.core.mb.create_meshset()
        M.core.mb.add_entities(f,M.core.all_faces[:])
        M.core.mb.write_file("results/solution_faces_all.vtk",[f])

    def solve(self,Mat, rhs):
        Mat[0]=np.zeros(self.nv+self.nfi)
        Mat[0][0]=1
        
        Mat[self.nv-1]=np.zeros(self.nv+self.nfi)
        Mat[self.nv-1][self.nv-1]=1
        sol=np.linalg.solve(self.M,self.rhs)
        return sol

    def write_output(self,M):
        print(self.sol)
        N=len(M.volumes.all)+len(M.faces.internal)
        nx=len(self.fx)
        ny=len(self.fy)
        nz=len(self.fz)
        volumes=M.volumes.all
        M.pressure[volumes]=self.sol[volumes]

        M.velocity[self.fx]=np.array([self.sol[self.nv:self.nv+nx],np.zeros(nx),np.zeros(nx)]).T
        M.velocity[self.fy]=np.array([np.zeros(ny),self.sol[self.nv+nx:self.nv+nx+ny],np.zeros(ny)]).T
        M.velocity[self.fz]=np.array([np.zeros(nz),np.zeros(nz),self.sol[self.nv+nx+ny:self.nv+nx+ny+nz]]).T
        v=M.core.mb.create_meshset()
        M.core.mb.add_entities(v,M.core.all_volumes)
        M.core.mb.write_file("results/solution_volumes.vtk",[v])

        f=M.core.mb.create_meshset()
        M.core.mb.add_entities(f,M.core.all_faces)
        M.core.mb.write_file("results/solution_faces.vtk",[f])
