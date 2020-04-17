

class bipahsicMultiscaleProperties:

    def flux_internal_faces():
        doc = "The flux_internal_faces property."
        def fget(self):
            return self._flux_internal_faces
        def fset(self, value):
            self._flux_internal_faces = value
        def fdel(self):
            del self._flux_internal_faces
        return locals()
    flux_internal_faces = property(**flux_internal_faces())
