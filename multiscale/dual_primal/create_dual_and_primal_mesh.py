import directories as direc
import utils_old

class DualPrimalMesh:

    def __init__(self, M):
        self.__loaded = False
        M.dualprimal = self
        self.tags = dict()
        self.info_tags = dict()
        self.entity_to_tag = dict()

    def create_tags(self, M):
        if self.__loaded:
            return 0





class DualPrimal:
    '''
    '''

    def __init__(self):
        self.Crs = direc.data_loaded[direc.names_data_loaded_lv0[3]]
