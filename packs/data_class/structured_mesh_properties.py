import numpy as np

class StructuredMeshProperties:

    def delta_h_internal_faces():
        doc = "The delta_h_internal_faces property."
        def fget(self):
            internal_faces = self.elements_lv0['internal_faces']
            return self.data_impress['dist_cent'][internal_faces]
        return locals()
    delta_h_internal_faces = property(**delta_h_internal_faces())

    def delta_z_internal_faces():
        doc = "The delta_z_internal_faces property."
        def fget(self):
            v0 = self.elements_lv0['neig_internal_faces']
            centroids = self.data_impress['centroid_volumes']
            zs = centroids[:, 2]
            return zs[v0[:, 1]] - zs[v0[:, 0]]
        return locals()
    delta_z_internal_faces = property(**delta_z_internal_faces())

    def delta_p_internal_faces():
        doc = "The delta_p_internal_faces property."
        def fget(self):
            x = self.data_impress['pressure']
            v0 = self.elements_lv0['neig_internal_faces']
            ps0 = x[v0[:, 0]]
            ps1 = x[v0[:, 1]]
            return (ps1 - ps0)
        return locals()
    delta_p_internal_faces = property(**delta_p_internal_faces())

    def grad_z_internal_faces():
        doc = "The grad_z_internal_faces property."
        def fget(self):
            return self.delta_z_internal_faces/self.delta_h_internal_faces
        return locals()
    grad_z_internal_faces = property(**grad_z_internal_faces())

    def grad_p_internal_faces():
        doc = "The delta_p_internal_faces property."
        def fget(self):
            return self.delta_p_internal_faces/self.delta_h_internal_faces
        return locals()
    grad_p_internal_faces = property(**grad_p_internal_faces())

    def rmap_internal_faces(self, value):

        resp = self.elements_lv0['remaped_internal_faces'][value]
        if isinstance(value, int):
            if set([-1]) & set([resp]):
                raise ValueError(f'{value} is not an internal face')
        elif set([-1]) & set(resp):
            raise ValueError(f'the entry may be an array of internal faces only')
        return resp

    def up_g():
        doc = "The up_g property."
        def fget(self):
            try:
                return self._up_g
            except AttributeError:
                internal_faces = self.elements_lv0['internal_faces']
                v0 = self.elements_lv0['neig_internal_faces']
                up_g = np.zeros(len(internal_faces), dtype=int)
                delta_z_internal_faces = self.delta_z_internal_faces
                up1 = delta_z_internal_faces >= 0
                up0 = delta_z_internal_faces < 0
                up_g[up1] = v0[up1, 1]
                up_g[up0] = v0[up0, 0]
                self._up_g = up_g
                return self._up_g
        return locals()
    up_g = property(**up_g())
