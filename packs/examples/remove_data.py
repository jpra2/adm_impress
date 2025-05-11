import os
from packs import defpaths
from typing import List
import shutil


def run_data_mesh():
    files: List[str] = os.listdir(defpaths.data_mesh)
    remove_exts = ['.npz']

    for arq in files:
        remove_file = False
        for rm_ext in remove_exts:
            if arq.endswith(rm_ext):
                remove_file = True
                break
        
        if remove_file is True:
            path = os.path.join(defpaths.data_mesh, arq)
            os.remove(path)
            print(f'{path} deleted \n')

def run_flying():
    files: List[str] = os.listdir(defpaths.flying)
    remove_exts = ['.npz', '.h5']

    for arq in files:
        remove_file = False
        for rm_ext in remove_exts:
            if arq.endswith(rm_ext):
                remove_file = True
                break
        
        if remove_file is True:
            path = os.path.join(defpaths.flying, arq)
            os.remove(path)
            print(f'{path} deleted \n')

def run_results():
    my_path = defpaths.results
    files: List[str] = os.listdir(my_path)
    # remove_exts = ['.npz', '.h5']
    remove_exts = ['.vtk']

    for arq in files:
        remove_file = False
        for rm_ext in remove_exts:
            if arq.endswith(rm_ext):
                remove_file = True
                break
        
        if remove_file is True:
            path = os.path.join(my_path, arq)
            os.remove(path)
            print(f'{path} deleted \n')
            
def delete_dual_infos():
    my_path = defpaths.results
    files: List[str] = os.listdir(my_path)
    dual_str = 'dual'

    for file in files:
        if dual_str in file:
            path = os.path.join(my_path, file)
            shutil.rmtree(path)

def run_delete_all():
    run_flying()
    run_data_mesh()
    run_results()
    delete_dual_infos()