import os
import shutil

path_to_execute = os.path.join('packs', 'utils', 'to_compile_funcs', 'scripts')
path_to_move = os.path.join('packs', 'utils', 'to_compile_funcs', 'funcs')

files = os.listdir(path_to_execute)

for file in files:
    juntar = os.path.join(path_to_execute, file)
    os.system(f'python {juntar}')

try:
    shutil.rmtree(path_to_move)
except FileNotFoundError:
    pass
os.makedirs(path_to_move)

all_files = os.listdir(path_to_execute)
ext = '.so'

my_files = [file for file in all_files if file.endswith(ext)]

for file in my_files:
    file_to_move = os.path.join(path_to_execute, file)
    shutil.move(file_to_move, path_to_move)