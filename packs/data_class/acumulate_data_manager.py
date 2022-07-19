import h5py
from packs.directories import only_mesh_name
from packs.directories import flying
import os
import pdb
import numpy as np
from packs.data_class.accumulative_array_data_manager import AccumulativeArrayDataManager
import pandas as pd

class CumulativeDatamanager(AccumulativeArrayDataManager):
    
    def export(self):
        name_export = self.file_name + str(self.global_identifier) + self.ext
        i=0
        with h5py.File(os.path.join(flying, name_export), 'w') as f:
            for data in self._data:
                grp = f.create_group(str(i))
                for key, value in data.items():
                    grp.create_dataset(key, data=value)
                i += 1

        self._data = []
        self.global_identifier += 1

        print(f'\n{name_export} exported\n')

    def load_all_datas(self):
        arqs_name = self.get_files_with_name()
        all_datas = []

        # import pdb; pdb.set_trace()
        for file in arqs_name:
            if file.startswith(self.file_name):
                with h5py.File(os.path.join(flying, file), 'r') as f:
                    for group in f.keys():
                        grp = f[group]
                        data_set = dict()
                        for key in grp.keys():
                            data_set[key] = grp[key][:]
                        all_datas.append(data_set)
        
        return all_datas
    
    def load_all_datas_from_keys(self, keyword_list):
        arqs_name = self.get_files_with_name()
        all_datas = []

        # import pdb; pdb.set_trace()
        for file in arqs_name:
            if file.startswith(self.file_name):
                with h5py.File(os.path.join(flying, file), 'r') as f:
                    for group in f.keys():
                        grp = f[group]
                        data_set = dict()
                        for key in grp.keys():
                            if key in keyword_list:
                                data_set[key] = grp[key][:]
                        all_datas.append(data_set)
        
        return all_datas

    def get_files_with_name(self):
        arqs = pd.Series(os.listdir(flying))
        test = arqs.str.startswith(self.file_name)
        result = arqs.values[test.values]
        return result
        