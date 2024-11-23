import pdb
from packs.data_class import DataManager
import numpy as np
import os


def test1():
    print(DataManager.all_datas.keys())
    a = DataManager()
    a['p'] = np.repeat(0, 3)
    print(a['p'])
    a.export_to_npz()
    assert os.path.exists(a.name) == True
    del a

def test2():
    print(DataManager.all_datas.keys())
    b = DataManager(load=True)
    print(b['p'])

test1()
test2()
