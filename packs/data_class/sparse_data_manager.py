import os
from packs.directories import flying
from packs.directories import only_mesh_name
from packs.data_class.common_data_manager import CommonDataManager
import numpy as np
import h5sparse
import weakref
import pdb


class SparseDataManager(CommonDataManager):

    '''
        Dataset for sparse matrix manager
    '''
    all_datas = dict()

    def __init__(self, data_name='', load=False):

        if data_name == '':
            data_name = type(self).__name__ + '_' + only_mesh_name + '.h5'

        self.name = os.path.join(flying, data_name)

        if not isinstance(data_name, str):
            raise ValueError('data_name must be string')

        if data_name[-3:] != '.h5':
            raise NameError('data_name must end with ".h5"')

        if self.name in SparseDataManager.all_datas.keys():
            raise ValueError('data_name cannot be repeated')
        self._data = dict()

        SparseDataManager.all_datas[self.name] = weakref.proxy(self)
        if load == True:
            self.load()
        elif load == False:
            pass
        else:
            raise NameError('load must be True or False')

    def load(self):

        with h5sparse.File(self.name, 'r') as f:
            names = f.keys()
            for name in names:
                self._data[name] = f[name].value

        print(f'\n{self.name} loaded\n')

    def export(self):
        with h5sparse.File(self.name, 'w') as f:
            names = self._data.keys()
            for name in names:
                f.create_dataset(name, data=self[name])

        print(f'\n{self.name} saved\n')

    @classmethod
    def get_obj_by_name(cls, name):
        return SparseDataManager.all_datas[name]

    @classmethod
    def export_all_datas(cls):
        for obj in SparseDataManager.all_datas.values():
        # for obj in cls.all_datas.values():
            obj.export()

    @classmethod
    def load_all_datas(cls):
        # for obj in cls.all_datas.values():
        for obj in SparseDataManager.all_datas.values():
            obj.load()

    def __del__(self):
        # if self.name == 'flying/CompositionalFVM.npz':
        #     import pdb; pdb.set_trace()
        try:
            SparseDataManager.all_datas = self.removekey(SparseDataManager.all_datas, self.name)
        except:
            pass
