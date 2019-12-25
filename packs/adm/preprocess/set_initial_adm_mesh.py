

class InitialADMMesh:

    def __init__(self, wells, levels, M, data_impress):
        self.mesh = M
        self.wells = wells
        self.levels = levels
        self.ids_volumes = M.volumes.all
        self.data_impress = data_impress

    def set_coarse_ids_classic(self):
        pass
