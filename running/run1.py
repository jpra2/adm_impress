
import os
import directories as direc
import numpy as np

__all__ = []
name0 = direc.names_outfiles_steps[0]

def load_mesh(name_mesh = name0):

    path_ant = os.getcwd()
    from impress.preprocessor_load import init_mesh
    M = init_mesh(name_mesh)
    M.data.init_dicts()
    M.data.load_variables_from_npz(direc.names_outfiles_variables_steps[0])
    M.data.update_variables_to_mesh()
    return M
