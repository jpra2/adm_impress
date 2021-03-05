import os
from ..directories import flying
import numpy as np

class DataManager:
    all_datas = dict()

    def __init__(self, data_name: str, load: bool=False) -> None:

        self._loaded = False
        if not isinstance(data_name, str):
            raise ValueError('data_name must be string')

        if data_name[-4:] != '.npz':
            raise NameError('data_name must end with ".npz"')

        if data_name in self.__class__.all_datas.keys():
            raise ValueError('data_name cannot be repeated')

        self.name = os.path.join(flying, data_name)
        self._data = dict()
        DataManager.all_datas[self.name] = self
        if load:
            self.load_from_npz()

    def export_to_npz(self):

        np.savez(self.name, **self._data)

        # with open(self.name_info_data, 'rb') as f:
        #     pickle.dump

    def load_from_npz(self):

        arq = np.load(self.name, allow_pickle=True)

        for name, variable in arq.items():
            self._data[name] = variable

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def export_all_datas_to_npz(self):
        for obj in DataManager.all_datas.values():
            obj.export_to_npz()

    def load_all_datas_from_npz(self):
        for obj in DataManager.all_datas.values():
            obj.load_from_npz()

    def get_obj_by_name(self, name):
        return DataManager.all_datas[name]

    def __str__(self):
        return str(type(self))

    def __setitem__(self, key, value):
        self._data[key] = value

    def __getitem__(self, key):
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __hash__(self, key):
        return hash(self._data)

    def __contains__(self, key):
        return key in self._data

    def __del__(self):
        # if self.name == 'flying/CompositionalFVM.npz':
        #     import pdb; pdb.set_trace()
        try:
            del DataManager.all_datas[self.name]
        except:
            pass
