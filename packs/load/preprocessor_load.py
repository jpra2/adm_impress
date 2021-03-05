from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as msh
import pdb

def init_mesh(mesh_name):

    M = msh(mesh_name, dim = 2)
    # M = msh('mesh/malha03.msh', dim = 2)

    return M
