
import os
from .. import directories as direc
from ..preprocess.data import Data

__all__ = []
name0 = direc.names_outfiles_steps[0]

def load_mesh(name_mesh = name0):

    path_ant = os.getcwd()
    from ..load.preprocessor_load import init_mesh
    M = init_mesh(name_mesh)
    from ..preprocess.init_data_class import initDataClass
    initDataClass(M)

    M.data.init_dicts()
    M.data.load_variables_from_npz(direc.names_outfiles_variables_steps[0])
    M.data.update_variables_to_mesh()

    return M
