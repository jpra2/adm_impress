import yaml
import os
from ..directories import flying

class InfoManager:
    all_files = dict()

    def __init__(self, from_file_name1: str, to_file_name2: str):
        # import pdb; pdb.set_trace()
        self.file_name = to_file_name2
        self.from_file_name1 = from_file_name1

        if set([self.file_name]) & set(InfoManager.all_files.keys()):
            raise NameError(f'\n O arquivo {self.file_name} ja foi carregado\n')
        self.load()
        InfoManager.all_files[self.file_name] = self

    def load(self):
        with open(self.from_file_name1, 'r') as f:
            self._data_loaded = yaml.safe_load(f)

    def save_obj(self):
        with open(self.file_name, 'w') as yaml_file:
            yaml.dump(self._data_loaded, yaml_file)

    @classmethod
    def get_obj_by_name(cls, name):
        return cls.all_files[name]

    def __setitem__(self, key, value):
        self._data_loaded[key] = value

    def __getitem__(self, key):
        return self._data_loaded[key]

    def __hash__(self, key):
        return hash(self._data_loaded[key])

    def __str__(self):
        return self.file_name

    def __contains__(self, key):
        return key in self._data_loaded
