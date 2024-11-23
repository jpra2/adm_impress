import pdb
from packs.data_class import SparseDataManager
import numpy as np
import os
from scipy import sparse as sp

class TestSparseData(SparseDataManager):
    pass

def test1():
    print(SparseDataManager.all_datas.keys())
    a = SparseDataManager()


    lines = np.array([1, 2, 2])
    cols = np.array([2, 3, 3])
    data = np.array([1, 1, 1])

    a['T'] = sp.csc_matrix((data, (lines, cols)), shape=(4, 4))
    print(a['T'])
    a.export()
    assert os.path.exists(a.name) == True
    del a

def test2():
    print(SparseDataManager.all_datas.keys())
    b = SparseDataManager(load=True)
    print(b['T'])

def test3():

    a = TestSparseData()
    lines = np.array([1, 2, 2, 3, 3])
    cols = np.array([2, 3, 3, 1, 1])
    data = np.array([1, 1, 1, 2, 2])

    a['T'] = sp.csc_matrix((data, (lines, cols)), shape=(5, 5))
    print(a['T'])

    b = SparseDataManager()

    lines = np.array([1, 2, 2, 3, 3])
    cols = np.array([2, 3, 3, 1, 1])
    data = np.array([1, 1, 1, 6, 6])

    b['T'] = sp.csc_matrix((data, (lines, cols)), shape=(5, 5))
    print(b['T'])

    b.export_all_datas()

# test1()
# test2()
test3()
