import numpy as np
import copy
from packs.errors import err as errors
import os
import packs.defpaths as defpaths
import pandas as pd

def test_array_instance(data):
    if isinstance(data, np.ndarray):
        pass
    else:
        raise TypeError

def test_str_instance(data):
    if isinstance(data, str):
        pass
    else:
        raise TypeError

class ArrayDataManager:
    
    def __init__(self, name=''):
        self.name = name
        self._data = dict()
    
    def insert_data(self, data):
        self._data.update(data)
    
    def export(self):
        np.savez(self.name, **self._data)
        print(f'\n{self.name} saved\n')
    
    def load(self):
        arq = np.load(self.name, allow_pickle=True)

        for name, variable in arq.items():
            self._data[name] = variable

        print(f'\n{self.name} loaded\n')
    
    def get_data_from_load(self):
        self.load()
        return copy.deepcopy(self._data)
    

class SuperArrayManager:
    
    def insert_name(self, name=''):
        self.__dict__['name'] = np.array([name])
    
    @classmethod
    def class_name(cls):
        return cls.__name__
    
    def insert_data(self, data: dict):
        """data is a dictionary with str keys and np.ndarray values

        Args:
            data (_type_): dict
        """
        names = list(data.keys())
        values = list(data.values())
        
        a = [test_array_instance(_) for _ in values]
        a = [test_str_instance(_) for _ in names]
        
        names_series = pd.DataFrame({
            'names': names
        })
        
        names_data_self = np.array(list(self.__dict__.keys()))
        names_data_self = names_data_self[names_data_self != 'mesh_name']
        test = names_series.isin(names_data_self)
        if test.any().values[0]:
            names_in = names_series[test.values].values.flatten()
            raise errors.NameExistsError(f'The names: - {names_in} - exists in mesh properties')
        
        
        self.__dict__.update(data)
    
    def __setattr__(self, name, value):
        raise Exception("It is read only! Use the 'insert_data' class function'")   
    
    @property
    def class_path(self):
        return os.path.join(defpaths.flying, self.class_name() + '_' +  self.name[0] + '.npz')
    
    def export_data(self):
        manager = ArrayDataManager(self.class_path)
        manager.insert_data(self.__dict__)
        manager.export()

    def load_data(self):
        if self.exists():
            pass
        else:
            raise FileNotFoundError
        
        manager = ArrayDataManager(self.class_path)
        self.insert_data(manager.get_data_from_load())

    def __getitem__(self, key):
        return self.__dict__[key]

    def get_all_data(self):
        return copy.deepcopy(self.__dict__)

    def rename_data(self, datas_to_rename: dict):
        """ Update the data name

        @param datas_to_rename: dict where key = old data name, value = new data name
        """

        new_data = dict()

        for name in list(datas_to_rename.keys()):
            data = copy.deepcopy(self[name])
            new_name = datas_to_rename[name]
            del self.__dict__[name]
            new_data.update({new_name: data})

        self.insert_data(new_data)
    
    def update_data(self, datas_to_update):
        
        my_datas_name = list(self.__dict__.keys())
        
        new_data = dict()
        for name in datas_to_update:
            if name not in my_datas_name:
                raise errors.NameExistsError
            new_data.update({
                name: datas_to_update[name]  
            })
            del self.__dict__[name]
        
        self.insert_data(new_data)
    
    def exists(self):
        return os.path.exists(self.class_path)
    
    def keys(self):
        return self.__dict__.keys()