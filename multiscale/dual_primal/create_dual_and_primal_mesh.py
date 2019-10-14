import directories as direc
from utils import utils_old

class DualPrimalMesh:

    def __init__(self, M):
        self.__loaded = False
        M.dualprimal = self
        self.tags = dict()
        self.mesh = M

    def create_tags(self, M):
        if self.__loaded:
            return 0

        names = []
