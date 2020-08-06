import os
from packs.directories import flying
from packs.directories import only_mesh_name
from packs.data_class.common_data_manager import CommonDataManager
import numpy as np
import weakref

class DataManager(CommonDataManager):
    all_datas = dict()

    def __init__(self, data_name: str='', load: bool=False) -> None:

        if data_name == '':
            data_name = type(self).__name__ + '_' + only_mesh_name + '.npz'

        self.name = os.path.join(flying, data_name)

        self._loaded = False
        if not isinstance(data_name, str):
            raise ValueError('data_name must be string')

        if data_name[-4:] != '.npz':
            raise NameError('data_name must end with ".npz"')

        if self.name in DataManager.all_datas.keys():
            raise ValueError('data_name cannot be repeated')

        self._data = dict()
        DataManager.all_datas[self.name] = weakref.proxy(self)
        if load:
            self.load_from_npz()

    def export_to_npz(self):

        np.savez(self.name, **self._data)

        print(f'\n{self.name} saved\n')

        # with open(self.name_info_data, 'rb') as f:
        #     pickle.dump

    def load_from_npz(self):

        arq = np.load(self.name, allow_pickle=True)

        for name, variable in arq.items():
            self._data[name] = variable

        print(f'\n{self.name} loaded\n')

    @classmethod
    def export_all_datas_to_npz(cls):
        for obj in DataManager.all_datas.values():
            obj.export_to_npz()

    @classmethod
    def load_all_datas_from_npz(cls):
        for obj in DataManager.all_datas.values():
            obj.load_from_npz()

    @classmethod
    def get_obj_by_name(cls, name):
        return DataManager.all_datas[name]

    def __str__(self):
        return str(type(self))

    def __del__(self):
        # if self.name == 'flying/CompositionalFVM.npz':
        #     import pdb; pdb.set_trace()

        try:
            DataManager.all_datas = self.removekey(DataManager.all_datas, self.name)
        except:
            pass
