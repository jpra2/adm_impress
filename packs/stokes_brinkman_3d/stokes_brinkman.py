import numpy as np
from packs.stokes_brinkman_3d.assembly import assembly
from packs.preprocess.preprocess_stokes_brinkman import preprocess_stokes
class stokes_solver:
    def __init__(self,M):
        prep_sb=preprocess_stokes(M)
        assembled_matrices=assembly.global_assembly(prep_sb,M)
        self.M = assembled_matrices.M
        self.rhs= assembled_matrices.rhs
        self.sol=self.solve(self.M,self.rhs, prep_sb)
        self.write_output(prep_sb, M)

    def solve(self,Mat, rhs, prep_sb):
        Mat[0]=np.zeros(prep_sb.nv+prep_sb.nfi)
        Mat[0][0]=1
        Mat[prep_sb.nv-1]=np.zeros(prep_sb.nv+prep_sb.nfi)
        Mat[prep_sb.nv-1][prep_sb.nv-1]=1
        sol=np.linalg.solve(self.M,self.rhs)
        return sol

    def write_output(self,prep_sb,M):
        print(self.sol)
        N=len(M.volumes.all)+len(M.faces.internal)
        nx=len(prep_sb.fx)
        ny=len(prep_sb.fy)
        nz=len(prep_sb.fz)
        volumes=M.volumes.all
        M.pressure[volumes]=self.sol[volumes]
        M.velocity[prep_sb.fx]=np.array([self.sol[prep_sb.nv:prep_sb.nv+nx],np.zeros(nx),np.zeros(nx)]).T
        M.velocity[prep_sb.fy]=np.array([np.zeros(ny),self.sol[prep_sb.nv+nx:prep_sb.nv+nx+ny],np.zeros(ny)]).T
        M.velocity[prep_sb.fz]=np.array([np.zeros(nz),np.zeros(nz),self.sol[prep_sb.nv+nx+ny:prep_sb.nv+nx+ny+nz]]).T

        v=M.core.mb.create_meshset()
        M.core.mb.add_entities(v,M.core.all_volumes)
        M.core.mb.write_file("results/solution_volumes.vtk",[v])

        f=M.core.mb.create_meshset()
        M.core.mb.add_entities(f,M.core.all_faces)
        M.core.mb.write_file("results/solution_faces.vtk",[f])
