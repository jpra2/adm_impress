import os
import pdb
from .. import directories
from ..preprocess.data import Data


__all__ = []

path_ant = os.getcwd()
os.chdir(directories.path_impress)
from impress.preprocessor0 import M

n_nodes = len(M.nodes)
n_faces = len(M.faces)
n_edges = len(M.edges)
n_volumes = len(M.volumes)

data = Data(n_nodes, n_faces, n_edges, n_volumes, M)
verif_format = ['float', 'int']

name_variables_yml = 'variable_settings.yml'
data_variables = M.read_config(name_variables_yml)

entitties = data_variables.keys()

for entity in entitties:
    try:
        dic1 = data_variables[entity]
    except:
        continue

    name_variables = dic1.keys()
    for name in name_variables:
        infos = dic1[name]
        data.get_info_data(name, infos['data size'], infos['data format'], entity)

M.data.init_datas()
M.data.init_dicts()
# from .impress.preprocessor0 import M
os.chdir(path_ant)
from ..preprocess.prep0 import Preprocess0
Preprocess0(M)
