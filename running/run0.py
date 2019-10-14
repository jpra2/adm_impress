import os
import directories

__all__ = []

path_ant = os.getcwd()
os.chdir(directories.path_impress)
from impress.preprocessor0 import M
M.data.init_datas()
M.data.init_dicts()
# from .impress.preprocessor0 import M
os.chdir(path_ant)
from preprocess.prep0 import Preprocess0
import pdb; pdb.set_trace()
Preprocess0(M)
