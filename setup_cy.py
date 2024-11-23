import os

path = os.path.join('packs', 'cython_files', 'mpfa_transmissibility')

os.chdir(path)

os.system(f'python setup.py')
