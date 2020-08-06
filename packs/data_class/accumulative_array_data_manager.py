import h5py
from packs.directories import only_mesh_name
from packs.directories import flying
import os
import pdb
import numpy as np

class AccumulativeArrayDataManager:

    def __init__(self):
        self.file_name = type(self).__name__ + '_' + only_mesh_name + '_'
        self.ext = '.h5'

    def create(self, global_identifier=0) -> None:

        _resps = ['s', 'n']
        assert isinstance(global_identifier, int) == True or isinstance(global_identifier, np.int) or isinstance(global_identifier, np.int64)

        data_name = self.file_name + str(global_identifier) + self.ext

        test = self.searc_for_data(data_name)
        if test == True:
            resp = input(f'\nVoce deseja sobrescrever e iniciar a partir do arquivo {data_name}\nDigite "s" ou "n" ')
            while resp not in _resps:
                resp = input('Digite "s" ou "n"')
        else:
            resp = 's'

        if resp == 's':
            pass
        else:
            global_identifier = self.find_identifier(self.file_name)
            data_name = self.file_name + str(global_identifier) + self.ext

        self.global_identifier = global_identifier
        self._data = []

    def insert_data(self, data):
        self._data.append(data)

    def export(self, local_key_identifier):
        name_export = self.file_name + str(self.global_identifier) + self.ext
        with h5py.File(os.path.join(flying, name_export), 'w') as f:
            for data in self._data:
                local_key = data[local_key_identifier]
                f.create_dataset(str(local_key), data=data)

        self._data = []
        self.global_identifier += 1

        print(f'\n{name_export} exported\n')

    def searc_for_data(self, name):

        tt = os.path.exists(name)
        return tt

    def find_identifier(self, name):

        all_datas = os.listdir(flying)
        indentifiers = []
        for file in all_datas:
            if file.startswith(name):
                identifier = file.replace(name, '')
                identifier = identifier.replace('.h5', '')
                identifiers.append(int(identifier))

        identifier = max(identifiers) + 1
        return identifier

    def load_all_datas(self):

        arqs_name = os.listdir(flying)
        all_datas = []

        for file in arqs_name:
            if file.startswith(self.file_name):
                with h5py.File(os.path.join(flying, file), 'r') as f:
                    for name in f.keys():
                        all_datas.append(f[name].value)

        return all_datas
