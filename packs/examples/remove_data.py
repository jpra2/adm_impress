import os
from packs import defpaths
from typing import List
import shutil

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
            

                 
