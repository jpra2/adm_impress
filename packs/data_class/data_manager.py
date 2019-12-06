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

        if data_name in DataManager.all_datas.keys():
            raise ValueError('data_name cannot be repeated')

        self.name = os.path.join(flying, data_name)
        self._data = dict()
        DataManager.all_datas[data_name] = self
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

        self._loaded = True

    def export_all_datas(self):
        for obj in DataManager.all_datas.values():
            obj.export_to_npz()

    def __str__(self):
        return str(list(self._data.keys()))

    def __setitem__(self, key, value):
        self._data[key] = value

    def __getitem__(self, key):
        return self._data[key]

    def __hash__(self, key):
        return hash(self._data[key])

# if __name__ == '__main__':
#
#     ff1 = dataManager('test.npz')
#     try:
#         ff = dataManager(1)
#     except Exception as e:
#         assert str(e) == 'data_name must be string'
#
#     try:
#         ff = dataManager('alguma_coisa')
#     except Exception as e:
#         assert str(e) == 'data_name must end with ".npz"'
#
#     try:
#         ff = dataManager('test.npz')
#     except Exception as e:
#         assert str(e) == 'data_name cannot be repeated'
#
#     print('\n Game over \n')
