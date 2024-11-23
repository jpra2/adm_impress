import sympy as sym

class LocalAssembly:
    def __init__(self):
        self.continuity()
        self.momentum()
    def continuity(self, dx, dy):
        
