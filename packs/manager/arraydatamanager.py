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
    
    def __init__(self):
        self.initialize_name()
    
    def initialize_name(self):
        """Modify initialize method if necessary for modify the name property
        """
        self.insert_name('')

    @classmethod
    def test_names(cls, names):
        """This function must be updated from the child class

        Args:
            names (_type_): _description_
        """
        pass
    
    def insert_name(self, name=''):
        self.__dict__['name'] = np.array([name])
    
    @classmethod
    def class_name(cls):
        return cls.__name__
    
    def insert_data(self, data: dict, load=False):
        """data is a dictionary with str keys and np.ndarray values

        Args:
            data (_type_): dict
        """
        
        names = list(data.keys())
        values = list(data.values())
        
        self.test_names(names)
        
        if load is True:
            indexes = np.arange(len(names))
            index_name = [i for i in indexes if names[i] == 'name']
            for i in index_name:
                names.pop(i)
                values.pop(i)
        
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
            raise errors.NameExistsError(f'The names: - {names_in} - exists in {self.class_name()}')
        
        
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
        self.verify_if_exists()
        
        manager = ArrayDataManager(self.class_path)
        self.insert_data(manager.get_data_from_load(), load=True)

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
    
    def remove_data(self, data_name: list):
        for name in data_name:
            del self.__dict__[name]

    def insert_or_update_data(self, datas: dict):
        
        data_names = list(datas.keys())
        names_in = self.verify_names_in_data_names(data_names)
        names_out = self.verify_names_out_data_names(data_names)
        
        self.insert_data({name: datas[name] for name in names_out})
        self.update_data({name: datas[name] for name in names_in})
    
    def exists(self):
        return os.path.exists(self.class_path)
    
    def verify_if_exists(self):
        if self.exists():
            pass
        else:
            raise FileExistsError
    
    def keys(self):
        return self.__dict__.keys()
    
    @property        
    def data_names(self):
        return list(self.keys())
    
    def verify_names_in_data_names(self, names: list):
        names_series = pd.DataFrame({
            'names': names
        })
        
        names_data_self = np.array(self.data_names)
        test = names_series.isin(names_data_self)
        test = test.values
        names_in = names_series[test].values.flatten()
        return names_in
    
    def verify_names_out_data_names(self, names:list):

        names_series = pd.DataFrame({
            'names': names
        })
        
        names_data_self = np.array(self.data_names)
        test = names_series.isin(names_data_self)
        test = ~test.values
        names_out = names_series[test].values.flatten()
        return names_out
    
    def verify_name_in_data_names(self, name: str):
        return name in self.data_names   

    def verify_name_in_data_names_or_raise_error(self, name: str):
        if self.verify_name_in_data_names(name):
            pass
        else:
            raise errors.NameExistsError(f'The name: - {name} - does not exists in mesh properties')
    
    def verify_name_not_in_data_names_or_raise_error(self, name: str):
        if self.verify_name_in_data_names(name):
            raise errors.NameExistsError(f'The name: - {name} - exists in mesh properties')


