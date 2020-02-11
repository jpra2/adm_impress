
class preprocess_stokes:
    def __init__(self,M):
        self.get_mesh_properties(M)
        self.nv=len(M.volumes.all)
        self.nfi=len(M.faces.internal)
        self.initiate_f_int_tag(M)

    def get_mesh_properties(self,M):
        v0=M.volumes.all[0]
        vert_v0=M.volumes.bridge_adjacencies(v0,0,0)
        coords_vert=M.nodes.coords(vert_v0)
        self.dx, self.dy, self.dz=coords_vert.max(axis=0)-coords_vert.min(axis=0)

    def initiate_f_int_tag(self, M):
        faces=M.faces.internal
        nodes=M.faces.bridge_adjacencies(faces,0,0).flatten()
        coords=M.nodes.coords[nodes].reshape(len(faces),4,3)

        cx=coords[:,:,0]
        cy=coords[:,:,1]
        cz=coords[:,:,2]
        fx=faces[(cx.max(axis=1)-cx.min(axis=1))<10**-5]
        fy=faces[(cy.max(axis=1)-cy.min(axis=1))<10**-5]
        fz=faces[(cz.max(axis=1)-cz.min(axis=1))<10**-5]

        self.fx=fx
        self.fy=fy
        self.fz=fz

        n_x=len(fx)
        n_y=len(fy)
        n_z=len(fz)
        if n_x>0:
            M.id_fint[fx] = range(n_x)
        if n_y>0:
            M.id_fint[fy] = range(n_x,n_x+n_y)
        if n_z>0:
            M.id_fint[fz] = range(n_x+n_y,n_x+n_y+n_z)
