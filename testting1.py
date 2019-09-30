import os
path_ant = os.getcwd()
from adm import directories
os.chdir(directories.path_impress)
from impress.preprocessor0 import M
os.chdir(path_ant)
M.data.init_datas()


import pdb; pdb.set_trace()

# import pickle
#
# from tcc.load_save_initialize.load_infos import LoadInfos
# from tcc.dual_mesh.create_dual_mesh import DualMesh1
# import numpy as np
# from . import directories
# import pdb; pdb.set_trace()
#
# # file_name = os.path.join(path_flying, 'mesh_obj.txt')
# # with open(file_name, 'wb') as handle:
# #     pickle.dump(M, handle)
#
#
#
# LoadInfos(M)
# DualMesh1(M)
#
# import pdb; pdb.set_trace()
