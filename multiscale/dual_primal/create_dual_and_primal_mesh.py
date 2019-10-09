

class DualPrimalMesh:

    def __init__(self, M):
        self.mesh = M
        M.dualprimal = self

    
