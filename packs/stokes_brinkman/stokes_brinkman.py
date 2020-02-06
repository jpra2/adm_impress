
from packs.stokes_brinkman.assembly import assembly
class stokes_solver:
    def __init__(self, M):
        dx, dy, dz = self.get_mesh_properties(M)
        self.M_cont = assembly.global_assembly(M, dx, dy, dz)

    def get_mesh_properties(self,M):
        v0=M.volumes.all[0]
        vert_v0=M.volumes.bridge_adjacencies(v0,0,0)
        coords_vert=M.nodes.coords(vert_v0)
        dx, dy, dz=coords_vert.max(axis=0)-coords_vert.min(axis=0)
        return dx, dy, dz
