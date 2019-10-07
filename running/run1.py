
import os
import directories as direc
import numpy as np

__all__ = []

path_ant = os.getcwd()
from impress.preprocessor_load import init_mesh
last_file_name = direc.names_outfiles_steps[0]
M = init_mesh(last_file_name)
# M.data.init_datas()
import time

t0 = time.time()
M.data.init_dicts()
M.data.load_variables_from_npz(direc.names_outfiles_variables_steps[0])
M.data.update_variables_to_mesh()
t1 = time.time()
print(f'tempo para carregar variaveis: {t1 - t0}')
import pdb; pdb.set_trace()
