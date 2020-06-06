import pdb
import numpy as np
from packs.data_class.data_manager import DataManager
from packs.directories.data_loaded import only_mesh_name


class MpfaOMonophasic(DataManager):

    def __init__(self, data_name: str='MpfaOMonophasic', load=False):
        data_name = data_name + '_' + only_mesh_name + '.npz'
        super().__init__(data_name=data_name, load=load)
