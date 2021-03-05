import copy

class InfosForProcess:
    def __init__(self,
        T: 'global transmissibility matrix without boundary conditions',
        pms: 'global multiscale presure',
        g_flux_grav_faces,
        gids: 'global gids',
        g_faces: 'global_faces',
        g_neig_internal_faces: 'all neig internal faces',
        remaped_internal_faces,
        solver
    ):

        self.T = T
        self.pms = pms
        self.g_flux_grav_faces = g_flux_grav_faces
        self.gids = gids
        self.g_faces = g_faces
        self.g_neig_internal_faces = g_neig_internal_faces
        self.remaped_internal_faces = remaped_internal_faces
        self.solver = solver

    def copy(self):
        return copy.deepcopy(self)
